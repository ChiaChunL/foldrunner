"""Engine-neutral representation of a prediction job.

The model is *entity plus copy count*, never a flat list of chains. Six of the
eight supported engines are natively count-based: AlphaFold 3 and Boltz-2 accept
a list of ``id`` values for one entity, Protenix, AlphaFold Server and SeedFold
take ``count``/``copies``, and ColabFold encodes cardinality in the a3m header.
Only Chai-1 and AlphaFold2-Multimer need the entity expanded into repeated
records, which is a concern for those writers alone.

Chain identifiers deliberately do not appear here. Protenix and AlphaFold Server
relabel every chain to A/B/C, and AlphaFold2-Multimer starts numbering at B, so a
chain name chosen up front would be fiction for three of the eight engines.
Writers assign identifiers and record what they assigned in the manifest.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from datetime import date

# The 20 standard residues, plus X for an unknown one. AlphaFold Server rejects
# anything outside the standard 20; other engines vary, so the check lives in the
# writers rather than here.
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
ALLOWED_AA = STANDARD_AA | {"X"}

_UNSAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


class IRError(ValueError):
    """Raised when a job or entity cannot be represented."""


def sanitize(name: str) -> str:
    """Reduce a label to characters that survive every engine and filesystem.

    Job names become file names for six of the eight engines, and ColabFold
    additionally uses them as FASTA headers, so the safe set is narrow.
    """
    cleaned = _UNSAFE_NAME.sub("_", name.strip()).strip("_")
    if not cleaned:
        raise IRError(f"name {name!r} is empty once reduced to safe characters")
    return cleaned


def sequence_key(sequence: str) -> str:
    """Content address for a sequence.

    SHA-256 of the uppercased sequence. This matches the naming Chai-1 expects
    for ``<hash>.aligned.pqt``, so one MSA cache can serve as its
    ``msa_directory`` without a rename step.
    """
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


@dataclass(frozen=True)
class Protein:
    """One unique protein sequence with a human-readable label."""

    name: str
    sequence: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", sanitize(self.name))
        seq = "".join(self.sequence.split()).upper()
        if not seq:
            raise IRError(f"protein {self.name!r} has an empty sequence")
        bad = sorted(set(seq) - ALLOWED_AA)
        if bad:
            raise IRError(
                f"protein {self.name!r} contains non-amino-acid characters: {''.join(bad)}"
            )
        object.__setattr__(self, "sequence", seq)

    @property
    def key(self) -> str:
        """SHA-256 content address, also the MSA cache key."""
        return sequence_key(self.sequence)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def has_unknown_residues(self) -> bool:
        return "X" in self.sequence


@dataclass(frozen=True)
class Ligand:
    """A small molecule, given either as SMILES or as a CCD code.

    AlphaFold Server accepts only a fixed CCD allow-list and rejects SMILES
    outright; Protenix and Boltz-2 accept either. Both fields are kept so a
    writer can pick whichever its engine supports and fail loudly otherwise.
    """

    name: str
    smiles: str | None = None
    ccd: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", sanitize(self.name))
        if not self.smiles and not self.ccd:
            raise IRError(f"ligand {self.name!r} needs either smiles or ccd")
        if self.ccd:
            object.__setattr__(self, "ccd", self.ccd.strip().upper().removeprefix("CCD_"))

    @property
    def key(self) -> str:
        return sequence_key(self.ccd or self.smiles or "")


Entity = Protein | Ligand


@dataclass(frozen=True)
class Component:
    """One entity together with how many copies of it the complex contains.

    A homodimer is ``Component(protein, count=2)`` — a single entity with two
    copies, not two components sharing a sequence. Every engine either expresses
    that directly or expands it in its own writer.
    """

    entity: Entity
    count: int = 1

    def __post_init__(self) -> None:
        if self.count < 1:
            raise IRError(f"component {self.entity.name!r} has count {self.count}, must be >= 1")

    @property
    def is_protein(self) -> bool:
        return isinstance(self.entity, Protein)


@dataclass(frozen=True)
class Job:
    """One complex to predict, in engine-neutral form."""

    name: str
    components: tuple[Component, ...]
    seeds: tuple[int, ...] = ()
    max_template_date: date | None = None
    affinity_binder: str | None = None
    notes: str = field(default="", compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", sanitize(self.name))
        if not self.components:
            raise IRError(f"job {self.name!r} has no components")
        if self.affinity_binder is not None:
            ligand_names = {c.entity.name for c in self.components if not c.is_protein}
            if self.affinity_binder not in ligand_names:
                raise IRError(
                    f"job {self.name!r} names affinity binder {self.affinity_binder!r}, "
                    f"which is not one of its ligands ({sorted(ligand_names) or 'none'})"
                )

    @property
    def proteins(self) -> tuple[Protein, ...]:
        return tuple(c.entity for c in self.components if isinstance(c.entity, Protein))

    @property
    def ligands(self) -> tuple[Ligand, ...]:
        return tuple(c.entity for c in self.components if isinstance(c.entity, Ligand))

    @property
    def unique_protein_keys(self) -> tuple[str, ...]:
        """Cache keys for every distinct protein sequence, in component order."""
        seen: dict[str, None] = {}
        for protein in self.proteins:
            seen.setdefault(protein.key, None)
        return tuple(seen)

    @property
    def n_chains(self) -> int:
        """Total chain count after copies are expanded."""
        return sum(c.count for c in self.components)

    @property
    def n_protein_chains(self) -> int:
        return sum(c.count for c in self.components if c.is_protein)

    @property
    def n_residues(self) -> int:
        """Residue total across all protein copies.

        Ligands are excluded: they contribute one PAE token per heavy atom rather
        than per residue, so any residue-based sizing has to treat them apart.
        """
        return sum(c.entity.length * c.count for c in self.components if c.is_protein)

    @property
    def is_homomer(self) -> bool:
        """True when the complex is copies of a single protein entity."""
        return len(self.components) == 1 and self.components[0].is_protein and self.n_chains > 1

    @property
    def is_monomer(self) -> bool:
        return self.n_chains == 1

    def with_seeds(self, seeds: tuple[int, ...]) -> Job:
        from dataclasses import replace

        return replace(self, seeds=tuple(seeds))
