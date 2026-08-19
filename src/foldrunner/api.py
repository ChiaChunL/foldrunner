"""Programmatic interface.

The command line is a thin layer over these functions; anything it can do is
available here, in the same order: work out the cost, get the alignments, write
the inputs, then run them.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from foldrunner.config import Config
from foldrunner.config import load as load_config
from foldrunner.engines import ENGINES, get_engine
from foldrunner.engines.base import MsaView, WrittenJob
from foldrunner.enumerate import (
    dedupe,
    enumerate_jobs,
    load_fasta,
    parse_fasta,
    unique_sequences,
)
from foldrunner.ir import Job, Protein
from foldrunner.manifest import PANEL, rows_for, write_manifest, write_panel
from foldrunner.msa.backends.base import Backend
from foldrunner.msa.cache import MsaCache, Recipe
from foldrunner.msa.importer import ImportReport, import_tree
from foldrunner.run import Invocation, plan, write_scripts

__all__ = [
    "PanelCost",
    "cost",
    "import_alignments",
    "library",
    "panel",
    "run_plan",
    "search",
    "write",
]


def library(source: str | Path | Iterable[Protein]) -> list[Protein]:
    """Load a protein library, collapsing entries that share a sequence."""
    if isinstance(source, (str, Path)):
        proteins = load_fasta(source) if Path(source).exists() else parse_fasta(str(source))
    else:
        proteins = list(source)
    kept, _ = dedupe(proteins)
    return kept


def panel(
    source: str | Path | Iterable[Protein],
    mode: str = "all",
    *,
    preys: str | Path | Iterable[Protein] | None = None,
    seeds: tuple[int, ...] = (),
    max_template_date: date | None = None,
) -> list[Job]:
    """Enumerate the complexes for a library."""
    return enumerate_jobs(
        library(source),
        mode,
        preys=library(preys) if preys is not None else None,
        seeds=seeds,
        max_template_date=max_template_date,
    )


@dataclass
class PanelCost:
    """What a panel will cost before anything is submitted."""

    complexes: int
    sequences: int
    engines: int

    @property
    def searches_saved(self) -> int:
        """Searches avoided by addressing alignments by sequence, not by pair."""
        return max(0, self.complexes - self.sequences)

    @property
    def files(self) -> int:
        return self.complexes * self.engines


def cost(jobs: Iterable[Job], engines: Iterable[str] | None = None) -> PanelCost:
    jobs = list(jobs)
    return PanelCost(
        complexes=len(jobs),
        sequences=len(unique_sequences(jobs)),
        engines=len(list(engines)) if engines else len(ENGINES),
    )


def search(
    jobs: Iterable[Job],
    cache: str | Path | MsaCache,
    backend: Backend,
    *,
    paired: bool = True,
    batch_size: int = 20,
) -> int:
    """Compute the alignments a panel needs, one per unique sequence.

    Returns how many sequences were searched; those already cached under the
    backend's recipe are skipped, so this is safe to call again after an
    interruption.
    """
    store = cache if isinstance(cache, MsaCache) else MsaCache(cache)
    recipe = backend.recipe
    # Recorded before the first search, so an interrupted run still leaves a
    # cache that can say what produced it. A recipe id hashes every field, and
    # retyping a few flags does not reproduce one.
    store.write_recipe(recipe)
    sequences = unique_sequences(list(jobs))

    searched = 0
    for kind, want in (("unpaired.a3m", False), ("paired.a3m", True)):
        if want and not paired:
            continue
        todo = [p for key, p in sequences.items() if not store.has(recipe, key, kind)]
        for start in range(0, len(todo), batch_size):
            batch = todo[start : start + batch_size]
            for protein, result in zip(
                batch, backend.fetch([p.sequence for p in batch], paired=want), strict=True
            ):
                store.store(
                    recipe,
                    protein.key,
                    protein.sequence,
                    unpaired=result.unpaired,
                    paired=result.paired,
                    source=type(backend).__name__,
                )
                searched += 1
    return searched


def import_alignments(
    source: str | Path,
    cache: str | Path | MsaCache,
    recipe: Recipe,
    *,
    layout: str = "auto",
    link: bool = False,
    dry_run: bool = False,
) -> ImportReport:
    """Adopt alignments produced by an earlier run, re-addressed by content."""
    store = cache if isinstance(cache, MsaCache) else MsaCache(cache)
    return import_tree(source, store, recipe, layout=layout, link=link, dry_run=dry_run)


def write(
    jobs: Iterable[Job],
    outdir: str | Path,
    *,
    engines: Iterable[str] | None = None,
    cache: str | Path | MsaCache | None = None,
    recipe: Recipe | None = None,
    mode: str = "all",
) -> list[WrittenJob]:
    """Generate the input files, the manifest and the panel snapshot."""
    jobs = list(jobs)
    outdir = Path(outdir)
    names = list(engines) if engines else sorted(ENGINES)

    view = None
    if cache is not None:
        if recipe is None:
            raise ValueError("a cache needs the recipe its alignments were produced under")
        store = cache if isinstance(cache, MsaCache) else MsaCache(cache)
        view = MsaView(cache=store, recipe=recipe)

    written: list[WrittenJob] = []
    for name in names:
        writer = get_engine(name)
        written.extend(writer(job, outdir / name, view) for job in jobs)

    rows = [row for item in written for row in rows_for(item)]
    write_manifest(outdir / "manifest.tsv", rows)
    write_panel(
        outdir / PANEL,
        written,
        seeds=jobs[0].seeds if jobs else (),
        mode=mode,
        recipe=recipe.recipe_id if recipe else "",
        cache=str(cache) if cache else "",
    )
    return written


def run_plan(
    written: Iterable[WrittenJob],
    out_root: str | Path,
    *,
    config: Config | str | Path | None = None,
    engines: Iterable[str] | None = None,
    scripts: str | Path | None = None,
) -> list[Invocation]:
    """Build the invocations for a written panel, optionally writing scripts."""
    conf = config if isinstance(config, Config) else load_config(config)
    invocations = plan(written, conf, Path(out_root), engines)
    if scripts is not None:
        write_scripts(invocations, Path(scripts))
    return invocations
