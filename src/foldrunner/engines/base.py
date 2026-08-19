"""Shared machinery for engine writers.

A writer is a pure function of a job plus an optional alignment view, so every
format can be tested offline without a network call or a GPU.

Chain identifiers are assigned here rather than carried in the job, because the
engines disagree about them: Protenix and AlphaFold Server relabel everything to
A/B/C regardless of what the input said, and AlphaFold2-Multimer begins at B.
Whatever a writer assigns is recorded on the returned :class:`WrittenJob`, and
from there into the manifest, which is what lets results from different engines
be lined up afterwards.
"""

from __future__ import annotations

import string
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

from foldrunner.ir import Job, Protein
from foldrunner.msa.cache import PAIRED, UNPAIRED, MsaCache, Recipe


class EngineDependencyError(RuntimeError):
    """An engine needs an optional dependency that is not installed."""


def chain_labels(start: str = "A") -> Iterator[str]:
    """A, B, ... Z, AA, AB, ... beginning at ``start``."""
    letters = string.ascii_uppercase
    offset = letters.index(start.upper())
    width = 1
    while True:
        if width == 1:
            for index in range(offset, len(letters)):
                yield letters[index]
        else:
            for first in letters:
                for second in letters:
                    yield first + second
        width += 1


def assign_chains(job: Job, start: str = "A") -> list[list[str]]:
    """Chain identifiers per component, one list per component in job order."""
    labels = chain_labels(start)
    return [[next(labels) for _ in range(component.count)] for component in job.components]


@dataclass
class MsaView:
    """Answers alignment questions for one cache and recipe."""

    cache: MsaCache
    recipe: Recipe

    def unpaired(self, protein: Protein) -> Path | None:
        return self._path(protein, UNPAIRED)

    def paired(self, protein: Protein) -> Path | None:
        return self._path(protein, PAIRED)

    def _path(self, protein: Protein, kind: str) -> Path | None:
        candidate = self.cache.entry_dir(self.recipe, protein.key) / kind
        return candidate if candidate.is_file() else None

    def directory(self) -> Path:
        return self.cache.root / self.recipe.recipe_id


@dataclass
class WrittenJob:
    """What a writer produced, and how it named things."""

    job: Job
    engine: str
    path: Path
    chains: list[list[str]] = field(default_factory=list)
    msa_used: bool = False
    extra: dict[str, str] = field(default_factory=dict)

    def chain_rows(self) -> list[tuple[str, str, int, str, str]]:
        """Manifest rows: job, engine, component index, entity name, chain id."""
        rows = []
        pairs = zip(self.job.components, self.chains, strict=True)
        for index, (component, labels) in enumerate(pairs):
            for label in labels:
                rows.append((self.job.name, self.engine, index, component.entity.name, label))
        return rows


Writer = Callable[..., WrittenJob]
ENGINES: dict[str, Writer] = {}


def register(name: str) -> Callable[[Writer], Writer]:
    def wrap(writer: Writer) -> Writer:
        ENGINES[name] = writer
        return writer

    return wrap


def get_engine(name: str) -> Writer:
    try:
        return ENGINES[name]
    except KeyError:
        known = ", ".join(sorted(ENGINES)) or "none registered"
        raise KeyError(f"unknown engine {name!r}; available: {known}") from None
