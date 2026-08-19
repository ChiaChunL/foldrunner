"""Content-addressed alignment store.

Entries are keyed by ``sha256(sequence.upper())``, which is also the name Chai-1
expects for ``<hash>.aligned.pqt``, so a cache directory can be handed to it as
``msa_directory`` unchanged. Addressing by content rather than by input position
is what makes the store survive a reordered library, and lets an alignment
computed for one panel be reused by the next.

Every entry sits under a ``recipe`` — the backend, databases, tool version and
search parameters that produced it. Alignments from different recipes have
different depths, and depth moves ipTM, so a panel that silently mixes them
produces differences that cannot be attributed after the fact.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

from foldrunner.msa.a3m import Alignment, parse_a3m

UNPAIRED = "unpaired.a3m"
PAIRED = "paired.a3m"
META = "meta.json"
RECIPE = "recipe.json"

KINDS = (UNPAIRED, PAIRED)

GZIP_MAGIC = b"\x1f\x8b"
# A tar header carries this at offset 257 whether or not the stream is
# compressed on top.
TAR_MAGIC_OFFSET = 257
TAR_MAGICS = (b"ustar\x0000", b"ustar  \x00")


class CacheError(RuntimeError):
    pass


@dataclass(frozen=True)
class Recipe:
    """How a set of alignments was produced."""

    backend: str
    databases: tuple[str, ...] = ()
    tool_version: str = ""
    params: tuple[tuple[str, str], ...] = ()

    @property
    def recipe_id(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True, default=list)
        digest = hashlib.sha256(payload.encode()).hexdigest()[:12]
        return f"{self.backend}-{digest}"

    def describe(self) -> str:
        dbs = "+".join(self.databases) if self.databases else "unspecified"
        version = self.tool_version or "unspecified"
        return f"{self.backend} (databases: {dbs}, version: {version})"


@dataclass
class Entry:
    """One cached sequence."""

    key: str
    sequence: str
    recipe_id: str
    depth_unpaired: int = 0
    depth_paired: int = 0
    n_species: int = 0
    source: str = ""
    created: float = field(default_factory=time.time)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, sort_keys=True)


def looks_like_archive(blob: bytes) -> bool:
    """Whether a payload is really an archive rather than an error page.

    The services answer a throttled request with an HTML page under a normal
    status code; unpacking it raises ``tarfile.ReadError``, which reads like a
    damaged download rather than the rate limit it is.

    The check accepts both a gzip stream and a bare tar, because whether the
    caller sees compression depends on transport: the services send the tarball
    with ``Content-Encoding: gzip``, and an HTTP client that honours that header
    hands back an already-decompressed tar. Testing only for the gzip magic
    would reject every successful response.
    """
    if blob[:2] == GZIP_MAGIC:
        return True
    return blob[TAR_MAGIC_OFFSET : TAR_MAGIC_OFFSET + 8] in TAR_MAGICS


def looks_like_gzip(blob: bytes) -> bool:
    """Deprecated name kept for callers; prefer :func:`looks_like_archive`."""
    return looks_like_archive(blob)


class MsaCache:
    """Filesystem store for alignments."""

    def __init__(self, root: str | Path) -> None:
        self.root = Path(root)

    def entry_dir(self, recipe: Recipe | str, key: str) -> Path:
        recipe_id = recipe if isinstance(recipe, str) else recipe.recipe_id
        return self.root / recipe_id / key

    def has(self, recipe: Recipe | str, key: str, kind: str = UNPAIRED) -> bool:
        return (self.entry_dir(recipe, key) / kind).is_file()

    def read(self, recipe: Recipe | str, key: str, kind: str = UNPAIRED) -> Alignment:
        path = self.entry_dir(recipe, key) / kind
        if not path.is_file():
            raise CacheError(f"no {kind} cached for {key[:12]} under {self.root}")
        return parse_a3m(path.read_text())

    def path(self, recipe: Recipe | str, key: str, kind: str = UNPAIRED) -> Path:
        path = self.entry_dir(recipe, key) / kind
        if not path.is_file():
            raise CacheError(f"no {kind} cached for {key[:12]} under {self.root}")
        return path

    def write(self, recipe: Recipe | str, key: str, kind: str, text: str) -> Path:
        """Store one alignment atomically.

        The write goes to a temporary name and is renamed into place, so an
        interrupted download leaves no half-written entry that a later run would
        mistake for a complete one.
        """
        if kind not in KINDS:
            raise CacheError(f"unknown alignment kind {kind!r}, expected one of {KINDS}")
        directory = self.entry_dir(recipe, key)
        directory.mkdir(parents=True, exist_ok=True)
        target = directory / kind
        tmp = directory / f".{kind}.{os.getpid()}.tmp"
        tmp.write_text(text)
        os.replace(tmp, target)
        return target

    def write_meta(self, recipe: Recipe | str, entry: Entry) -> Path:
        directory = self.entry_dir(recipe, entry.key)
        directory.mkdir(parents=True, exist_ok=True)
        target = directory / META
        tmp = directory / f".{META}.{os.getpid()}.tmp"
        tmp.write_text(entry.to_json())
        os.replace(tmp, target)
        return target

    def read_meta(self, recipe: Recipe | str, key: str) -> Entry | None:
        path = self.entry_dir(recipe, key) / META
        if not path.is_file():
            return None
        return Entry(**json.loads(path.read_text()))

    def store(
        self,
        recipe: Recipe,
        key: str,
        sequence: str,
        *,
        unpaired: str | None = None,
        paired: str | None = None,
        source: str = "",
    ) -> Entry:
        """Store alignments plus the metadata describing them."""
        # Merge rather than replace: the two halves are fetched in separate
        # passes, and a second call that rewrote the record from scratch would
        # report the first half as empty even though its file is right there.
        entry = self.read_meta(recipe, key) or Entry(
            key=key, sequence=sequence, recipe_id=recipe.recipe_id, source=source
        )
        if source:
            entry.source = source
        if unpaired is not None:
            self.write(recipe, key, UNPAIRED, unpaired)
            entry.depth_unpaired = parse_a3m(unpaired).depth
        if paired is not None:
            self.write(recipe, key, PAIRED, paired)
            alignment = parse_a3m(paired)
            entry.depth_paired = alignment.depth
            entry.n_species = len({s for s in alignment.species if s})
        self.write_meta(recipe, entry)
        return entry

    def write_recipe(self, recipe: Recipe) -> Path:
        """Record what produced these alignments, next to them.

        A recipe id is a hash of the whole recipe, so it cannot be rebuilt from
        a few command-line flags without getting every field right — including
        ones the caller never typed. Storing it means a later step can ask the
        cache what it holds instead of asking the user to remember.
        """
        directory = self.root / recipe.recipe_id
        directory.mkdir(parents=True, exist_ok=True)
        target = directory / RECIPE
        tmp = directory / f".{RECIPE}.{os.getpid()}.tmp"
        tmp.write_text(json.dumps(asdict(recipe), indent=2, sort_keys=True, default=list))
        os.replace(tmp, target)
        return target

    def read_recipe(self, recipe_id: str) -> Recipe | None:
        path = self.root / recipe_id / RECIPE
        if not path.is_file():
            return None
        raw = json.loads(path.read_text())
        return Recipe(
            backend=raw["backend"],
            databases=tuple(raw.get("databases", ())),
            tool_version=raw.get("tool_version", ""),
            params=tuple(tuple(p) for p in raw.get("params", ())),
        )

    def stored_recipes(self) -> list[Recipe]:
        """Every recipe this cache can describe, newest layout first."""
        found = []
        for recipe_id in self.recipes():
            recipe = self.read_recipe(recipe_id)
            if recipe is not None:
                found.append(recipe)
        return found

    def missing(self, recipe: Recipe | str, keys: list[str], kind: str = UNPAIRED) -> list[str]:
        """Which of these sequences still need work — the basis for resuming."""
        return [k for k in keys if not self.has(recipe, k, kind)]

    def recipes(self) -> list[str]:
        if not self.root.is_dir():
            return []
        return sorted(p.name for p in self.root.iterdir() if p.is_dir())

    def keys(self, recipe: Recipe | str) -> list[str]:
        recipe_id = recipe if isinstance(recipe, str) else recipe.recipe_id
        directory = self.root / recipe_id
        if not directory.is_dir():
            return []
        return sorted(p.name for p in directory.iterdir() if p.is_dir())
