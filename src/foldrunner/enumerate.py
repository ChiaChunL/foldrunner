"""Turn a protein library into the set of complexes to predict.

Self-pairing means a homodimer: one entity with ``count=2``, not two components
that happen to share a sequence. Monomers are a separate mode, because a panel
usually wants them as controls rather than mixed in with the pairs.
"""

from __future__ import annotations

import itertools
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from foldrunner.ir import Component, IRError, Job, Protein

MODES = ("all", "hetero", "monomer", "bipartite")

JOIN = "__"


@dataclass(frozen=True)
class Duplicate:
    """Two labels in the library that carry the same sequence."""

    kept: str
    dropped: str
    key: str


def parse_fasta(text: str) -> list[Protein]:
    """Read a FASTA string into proteins, taking the first whitespace-delimited
    field of each header as the label."""
    proteins: list[Protein] = []
    label: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        if label is not None:
            proteins.append(Protein(label, "".join(chunks)))

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            label = header.split()[0].split("|")[0] if header else ""
            chunks = []
        elif label is not None:
            chunks.append(line)
    flush()
    if not proteins:
        raise IRError("no sequences found in FASTA input")
    return proteins


def load_fasta(path: str | Path) -> list[Protein]:
    return parse_fasta(Path(path).read_text())


def dedupe(proteins: Iterable[Protein]) -> tuple[list[Protein], list[Duplicate]]:
    """Collapse entries that share a sequence, keeping the first label seen.

    Identical sequences under different labels would otherwise be searched twice
    and produce two job names for one prediction.
    """
    by_key: dict[str, Protein] = {}
    duplicates: list[Duplicate] = []
    for protein in proteins:
        existing = by_key.get(protein.key)
        if existing is None:
            by_key[protein.key] = protein
        elif existing.name != protein.name:
            duplicates.append(Duplicate(kept=existing.name, dropped=protein.name, key=protein.key))
    return list(by_key.values()), duplicates


def pair_name(a: Protein, b: Protein) -> str:
    """Canonical job name for a pair, order-independent."""
    first, second = sorted((a.name, b.name))
    return f"{first}{JOIN}{second}"


def _pair_job(a: Protein, b: Protein, **kwargs: object) -> Job:
    if a.key == b.key:
        components = (Component(a, count=2),)
    else:
        first, second = sorted((a, b), key=lambda p: p.name)
        components = (Component(first), Component(second))
    return Job(name=pair_name(a, b), components=components, **kwargs)  # type: ignore[arg-type]


def enumerate_jobs(
    proteins: Iterable[Protein],
    mode: str = "all",
    *,
    preys: Iterable[Protein] | None = None,
    seeds: tuple[int, ...] = (),
    max_template_date: date | None = None,
) -> list[Job]:
    """Build the job list for a library.

    ``all`` gives every unordered pair including self-pairs, N(N+1)/2 jobs.
    ``hetero`` drops the self-pairs, N(N-1)/2. ``monomer`` gives one job per
    sequence. ``bipartite`` crosses ``proteins`` (baits) with ``preys``.
    """
    if mode not in MODES:
        raise IRError(f"unknown mode {mode!r}, expected one of {', '.join(MODES)}")

    baits = list(proteins)
    common = {"seeds": seeds, "max_template_date": max_template_date}

    if mode == "monomer":
        return [Job(name=p.name, components=(Component(p),), **common) for p in baits]

    if mode == "bipartite":
        if preys is None:
            raise IRError("bipartite mode needs a second library (preys)")
        jobs: dict[str, Job] = {}
        for bait, prey in itertools.product(baits, list(preys)):
            job = _pair_job(bait, prey, **common)
            jobs.setdefault(job.name, job)
        return list(jobs.values())

    combine = (
        itertools.combinations_with_replacement if mode == "all" else itertools.combinations
    )
    return [_pair_job(a, b, **common) for a, b in combine(baits, 2)]


def expected_job_count(n: int, mode: str = "all", n_preys: int = 0) -> int:
    """Job count for a library of ``n`` unique sequences, without building them."""
    if mode == "all":
        return n * (n + 1) // 2
    if mode == "hetero":
        return n * (n - 1) // 2
    if mode == "monomer":
        return n
    if mode == "bipartite":
        return n * n_preys
    raise IRError(f"unknown mode {mode!r}")


def unique_sequences(jobs: Iterable[Job]) -> dict[str, Protein]:
    """Every distinct protein across a job list, keyed by content address.

    This is the number that drives MSA cost: it grows with the library, not with
    the number of pairs, which is the whole point of computing alignments once.
    """
    found: dict[str, Protein] = {}
    for job in jobs:
        for protein in job.proteins:
            found.setdefault(protein.key, protein)
    return found


def iter_chunks(jobs: list[Job], size: int) -> Iterator[list[Job]]:
    """Split a job list into fixed-size batches.

    AlphaFold Server and SeedFold accept a list of jobs per upload but cap daily
    submissions, so their writers emit batches rather than one large file.
    """
    if size < 1:
        raise IRError(f"chunk size must be >= 1, got {size}")
    for start in range(0, len(jobs), size):
        yield jobs[start : start + size]
