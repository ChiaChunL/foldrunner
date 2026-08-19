"""What every alignment source has to provide."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from foldrunner.msa.cache import Recipe


class SearchError(RuntimeError):
    """The search could not be completed."""


class RateLimited(SearchError):
    """The service declined the request because of its usage limits.

    Distinguished from a general failure because it is worth waiting out: the
    same request usually succeeds a few minutes later.
    """


@dataclass
class MsaResult:
    """Alignments for one query sequence."""

    sequence: str
    unpaired: str | None = None
    paired: str | None = None


@runtime_checkable
class Backend(Protocol):
    """A source of alignments.

    Implementations take a batch of sequences rather than one at a time. Both
    the web services and a local mmseqs run amortise their cost over a batch —
    the databases are scanned once for the whole set — so searching one sequence
    per call is dramatically slower for no benefit.
    """

    @property
    def recipe(self) -> Recipe: ...

    def fetch(self, sequences: list[str], *, paired: bool = False) -> list[MsaResult]: ...
