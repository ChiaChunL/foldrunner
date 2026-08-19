"""Adopt alignments produced by earlier runs into the cache.

Existing pipelines name their output by the position of the sequence in the
input FASTA (``0/``, ``1/``, ``2/`` ...). That mapping lives outside the files,
so it breaks as soon as the input is reordered or the directory is renamed, and
it breaks silently. Re-addressing by content removes the dependency entirely.

Nothing outside the alignment is needed to do it: the first record of an a3m is
the query itself, so the cache key can be recomputed from the file. An import is
therefore a pure rewrite of an existing tree, and can be re-run safely.
"""

from __future__ import annotations

import os
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from foldrunner.ir import sequence_key
from foldrunner.msa.a3m import parse_a3m, read_query, scan
from foldrunner.msa.cache import PAIRED, UNPAIRED, Entry, MsaCache, Recipe

# Layouts seen in the wild. The Protenix MSA tool splits its output into a
# pairing/non-pairing pair per sequence; the ColabFold service returns one file
# per database searched, all of them unpaired.
PROTENIX_PAIRED = "pairing.a3m"
PROTENIX_UNPAIRED = "non_pairing.a3m"
COLABFOLD_UNPAIRED = ("uniref.a3m", "bfd.mgnify30.metaeuk30.smag30.a3m")

LAYOUTS = ("auto", "protenix", "colabfold")


@dataclass
class Found:
    """One importable entry located on disk."""

    directory: Path
    paired: Path | None = None
    unpaired: list[Path] = field(default_factory=list)
    layout: str = ""


@dataclass
class ImportReport:
    imported: int = 0
    skipped: int = 0
    failed: list[tuple[str, str]] = field(default_factory=list)
    keys: set[str] = field(default_factory=set)
    collisions: list[tuple[str, str]] = field(default_factory=list)
    conflicts: list[tuple[str, str, str]] = field(default_factory=list)
    scanned: int = 0
    layouts: dict[str, int] = field(default_factory=dict)

    def summary(self) -> str:
        lines = [
            f"scanned {self.scanned} entries, imported {self.imported}, "
            f"already cached {self.skipped}, failed {len(self.failed)}",
            f"unique sequences {len(self.keys)}",
        ]
        if self.collisions:
            share = 100 * len(self.collisions) / max(self.scanned, 1)
            lines.append(
                f"{len(self.collisions)} entries ({share:.0f}%) repeated a sequence already "
                "seen, which the position-based layout could not reveal"
            )
        if self.conflicts:
            kinds = sorted({f"{a}/{b}" for _, a, b in self.conflicts})
            lines.append(
                f"WARNING: {len(self.conflicts)} of those came from a different source "
                f"layout ({', '.join(kinds)}). Alignments from different searches have "
                "different depths and must not share a recipe. Import each source "
                "separately with its own --backend, or restrict this run with --layout."
            )
        if len(self.layouts) > 1:
            found = ", ".join(f"{name}={count}" for name, count in sorted(self.layouts.items()))
            lines.append(f"source layouts present: {found}")
        if self.failed:
            lines.append(f"first failure: {self.failed[0][0]} — {self.failed[0][1]}")
        return "\n".join(lines)


def _merge_unpaired(paths: list[Path]) -> str:
    """Join several single-database alignments into one unpaired block.

    Each file repeats the query as its first record; only the first copy is
    kept so the result stays a valid a3m.
    """
    chunks: list[str] = []
    for position, path in enumerate(paths):
        alignment = parse_a3m(path.read_text())
        start = 0 if position == 0 else 1
        for header, row in zip(
            alignment.headers[start:], alignment.rows[start:], strict=True
        ):
            chunks.append(f">{header}\n{row}\n")
    return "".join(chunks)


def discover(root: str | Path, layout: str = "auto") -> Iterator[Found]:
    """Walk a tree and yield every directory that holds importable alignments."""
    if layout not in LAYOUTS:
        raise ValueError(f"unknown layout {layout!r}, expected one of {', '.join(LAYOUTS)}")
    root = Path(root)
    for directory in sorted(p for p in root.rglob("*") if p.is_dir()):
        paired = directory / PROTENIX_PAIRED
        unpaired_protenix = directory / PROTENIX_UNPAIRED
        unpaired_colabfold = [
            directory / n for n in COLABFOLD_UNPAIRED if (directory / n).is_file()
        ]

        if layout in ("auto", "protenix") and (paired.is_file() or unpaired_protenix.is_file()):
            yield Found(
                directory=directory,
                paired=paired if paired.is_file() else None,
                unpaired=[unpaired_protenix] if unpaired_protenix.is_file() else [],
                layout="protenix",
            )
        elif layout in ("auto", "colabfold") and unpaired_colabfold:
            yield Found(
                directory=directory, unpaired=unpaired_colabfold, layout="colabfold"
            )


def _link_or_copy(source: Path, target: Path) -> None:
    """Place a file in the cache without duplicating its bytes when possible.

    Proteome-scale trees run to hundreds of gigabytes, so copying them to
    re-address the entries would double the footprint for no gain. A hard link
    keeps one copy on disk; it falls back to a real copy across filesystems.
    """
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        target.unlink()
    try:
        os.link(source, target)
    except OSError:
        target.write_bytes(source.read_bytes())


def import_entry(
    found: Found,
    cache: MsaCache,
    recipe: Recipe,
    *,
    dry_run: bool = False,
    link: bool = False,
):
    """Import one directory. Returns (key, was_written)."""
    reference = found.paired or (found.unpaired[0] if found.unpaired else None)
    if reference is None:
        raise ValueError("no alignment files in directory")

    sequence = read_query(reference).replace("-", "").upper()
    if not sequence:
        raise ValueError("query record is empty")
    key = sequence_key(sequence)

    # Both halves must describe the same query, otherwise the directory holds a
    # mismatched pair and importing it would attach the wrong alignment.
    if found.paired and found.unpaired:
        other = read_query(found.unpaired[0])
        if other.replace("-", "").upper() != sequence:
            raise ValueError("paired and unpaired alignments have different query sequences")

    if cache.has(recipe, key) or (found.paired and cache.has(recipe, key, "paired.a3m")):
        return key, False
    if dry_run:
        return key, True

    # A single source file can be linked verbatim; several have to be merged
    # into one block first, which means writing new bytes either way.
    single_unpaired = found.unpaired[0] if len(found.unpaired) == 1 else None
    if link and (single_unpaired or not found.unpaired):
        directory = cache.entry_dir(recipe, key)
        entry = Entry(key=key, sequence=sequence, recipe_id=recipe.recipe_id,
                      source=str(found.directory))
        if single_unpaired:
            _link_or_copy(single_unpaired, directory / UNPAIRED)
            entry.depth_unpaired, _ = scan(single_unpaired)
        if found.paired:
            _link_or_copy(found.paired, directory / PAIRED)
            entry.depth_paired, entry.n_species = scan(found.paired)
        cache.write_meta(recipe, entry)
        return key, True

    cache.store(
        recipe,
        key,
        sequence,
        unpaired=_merge_unpaired(found.unpaired) if found.unpaired else None,
        paired=found.paired.read_text() if found.paired else None,
        source=str(found.directory),
    )
    return key, True


def import_tree(
    root: str | Path,
    cache: MsaCache,
    recipe: Recipe,
    *,
    layout: str = "auto",
    dry_run: bool = False,
    link: bool = False,
    limit: int | None = None,
) -> ImportReport:
    """Re-address a whole tree of alignments into the cache."""
    report = ImportReport()
    origin: dict[str, str] = {}
    for count, found in enumerate(discover(root, layout)):
        if limit is not None and count >= limit:
            break
        report.scanned += 1
        report.layouts[found.layout] = report.layouts.get(found.layout, 0) + 1
        try:
            key, written = import_entry(found, cache, recipe, dry_run=dry_run, link=link)
        except Exception as error:  # noqa: BLE001 - reported per entry, never fatal
            report.failed.append((str(found.directory), str(error)))
            continue

        if key in report.keys:
            report.collisions.append((key, str(found.directory)))
            # A repeat from a different tool is not a harmless duplicate: the two
            # searches produce different depths, and silently keeping whichever
            # was walked first would put an arbitrary mixture in one recipe.
            if origin.get(key, found.layout) != found.layout:
                report.conflicts.append((key, origin[key], found.layout))
            report.skipped += 1
            continue

        report.keys.add(key)
        origin[key] = found.layout
        if written:
            report.imported += 1
        else:
            report.skipped += 1
    return report
