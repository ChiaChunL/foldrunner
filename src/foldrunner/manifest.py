"""The table that makes cross-engine comparison possible.

Three of the eight engines discard the chain names they are given: Protenix and
AlphaFold Server relabel everything to A/B/C, and AlphaFold2-Multimer starts at
B. Without a record of what each writer actually emitted, results coming back
from different engines cannot be lined up against the same entities.

The manifest is also what makes a run resumable, and it carries the paired MSA
depth per job, which is needed to read the scores: cross-species pairs have
almost no pairing, so their interface confidence is not on the same footing as a
same-species pair's.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from foldrunner.engines.base import WrittenJob

COLUMNS = (
    "job",
    "engine",
    "component",
    "entity",
    "chain",
    "count",
    "n_chains",
    "n_residues",
    "msa_used",
    "n_paired",
    "path",
)


@dataclass
class Row:
    job: str
    engine: str
    component: int
    entity: str
    chain: str
    count: int
    n_chains: int
    n_residues: int
    msa_used: bool
    n_paired: int
    path: str

    def as_tuple(self) -> tuple:
        return (
            self.job,
            self.engine,
            self.component,
            self.entity,
            self.chain,
            self.count,
            self.n_chains,
            self.n_residues,
            int(self.msa_used),
            self.n_paired,
            self.path,
        )


def rows_for(written: WrittenJob, n_paired: int = 0) -> list[Row]:
    out: list[Row] = []
    pairs = zip(written.job.components, written.chains, strict=True)
    for index, (component, labels) in enumerate(pairs):
        for label in labels:
            out.append(
                Row(
                    job=written.job.name,
                    engine=written.engine,
                    component=index,
                    entity=component.entity.name,
                    chain=label,
                    count=component.count,
                    n_chains=written.job.n_chains,
                    n_residues=written.job.n_residues,
                    msa_used=written.msa_used,
                    n_paired=n_paired,
                    path=str(written.path),
                )
            )
    return out


def write_manifest(path: str | Path, rows: Iterable[Row]) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow(row.as_tuple())
    return target


def read_manifest(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def completed_jobs(path: str | Path, engine: str) -> set[str]:
    """Job names already written for an engine, for resuming a run."""
    target = Path(path)
    if not target.is_file():
        return set()
    return {r["job"] for r in read_manifest(target) if r["engine"] == engine}


PANEL = "panel.json"


def _relative(value: str | Path, root: Path) -> str:
    """Express a path under the panel directory relative to it.

    Anything outside the panel — an alignment cache shared between runs, for
    instance — is left as it is, because it is not the panel's to move.
    """
    try:
        return str(Path(value).resolve().relative_to(root.resolve()))
    except (ValueError, OSError):
        return str(value)


def rebase(value: str, root: Path) -> str:
    """Turn a stored path back into one on this machine."""
    candidate = Path(value)
    return str(candidate if candidate.is_absolute() else (root / candidate))


def write_panel(
    path: str | Path,
    written: Iterable[WrittenJob],
    *,
    seeds: tuple[int, ...] = (),
    mode: str = "all",
    recipe: str = "",
    cache: str = "",
) -> Path:
    """Record everything `run` needs that the manifest cannot carry.

    The manifest is a table for people and for downstream scoring, so it holds
    one row per chain and only the columns that mean something there. Running
    the panel needs things that do not fit that shape — the seeds the jobs were
    written with, and the per-engine extras such as the path of a generated
    driver script. Reconstructing an invocation from the table alone silently
    drops both.

    Paths are stored relative to the panel directory. An absolute path names the
    machine the panel was written on, and running it elsewhere — which is the
    whole point of being able to configure a host — would produce commands
    pointing at a directory that does not exist there.
    """
    target = Path(path)
    root = target.parent
    root.mkdir(parents=True, exist_ok=True)
    payload = {
        "seeds": list(seeds),
        "mode": mode,
        "recipe": recipe,
        "cache": cache,
        "written": [
            {
                "engine": item.engine,
                "job": item.job.name,
                "path": _relative(item.path, root),
                "chains": item.chains,
                "msa_used": item.msa_used,
                "extra": {k: _relative(v, root) for k, v in item.extra.items()},
                "components": [
                    {"entity": c.entity.name, "count": c.count, "protein": c.is_protein}
                    for c in item.job.components
                ],
            }
            for item in written
        ],
    }
    target.write_text(json.dumps(payload, indent=2))
    return target


def read_panel(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())
