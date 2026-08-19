"""a3m parsing and taxonomy-based pairing.

Pairing exists to give the model cross-chain coevolution: a row is only
informative if both halves come from the same organism. That makes a
taxonomy-labelled source mandatory. Hits returned by the pairing endpoints carry
the taxon in the header (``UniRef100_<accession>_<taxid>/``); metagenomic hits
used for unpaired alignments carry none, which is why the two are fetched and
stored separately.

The grouping rules follow AlphaFold2-Multimer (``alphafold/data/msa_pairing.py``)
and Protenix (``protenix/data/msa/msa_utils.py``): group by taxon, drop taxa that
appear in only one chain, and rank within a taxon by similarity to that chain's
own query.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

# UniRef hits arrive as "UniRef100_<accession>_<taxid>/", UniProt hits as
# "tr|<accession>|<name>_<SPECIES>". Both forms are matched so a cache can hold
# alignments from either source.
_UNIREF_RE = re.compile(r"^UniRef\d*_[^_\s]+_([^_/\s]+)")
_UNIPROT_RE = re.compile(r"^(?:tr|sp)\|[A-Z0-9]{6,10}(?:_\d+)?\|[A-Z0-9]{1,10}_([A-Z0-9]{1,5})")

GAP = "-"


def species_id(header: str) -> str | None:
    """Taxon for one a3m header, or None when the source carries no taxonomy."""
    header = header.lstrip(">").strip()
    match = _UNIREF_RE.match(header) or _UNIPROT_RE.match(header)
    return match.group(1) if match else None


def ungapped_length(row: str) -> int:
    """Length of a row once a3m insertions (lowercase) are removed.

    Every row in a valid a3m must reduce to the query length; engines reject
    alignments where it does not.
    """
    return sum(1 for ch in row if not ch.islower())


@dataclass
class Alignment:
    """One query sequence and its hits, in a3m order with the query first."""

    headers: list[str]
    rows: list[str]

    def __post_init__(self) -> None:
        if not self.rows:
            raise ValueError("alignment has no rows")
        if len(self.headers) != len(self.rows):
            raise ValueError(
                f"alignment has {len(self.headers)} headers but {len(self.rows)} rows"
            )

    @property
    def query(self) -> str:
        return self.rows[0]

    @property
    def depth(self) -> int:
        return len(self.rows)

    @property
    def species(self) -> list[str | None]:
        return [species_id(h) for h in self.headers]

    @property
    def n_with_species(self) -> int:
        return sum(1 for s in self.species if s)

    def gap_row(self) -> str:
        return GAP * len(self.query)


def read_query(path: Path) -> str:
    """Read only the first record of an a3m file.

    The query is all that is needed to work out an entry's content address, and
    alignments run to tens of megabytes each. At proteome scale, parsing whole
    files to recover a hash means reading hundreds of gigabytes to produce a few
    kilobytes of keys.
    """
    with Path(path).open() as handle:
        header = handle.readline()
        if not header.startswith(">"):
            raise ValueError(f"{path} does not start with an a3m header")
        rows: list[str] = []
        for line in handle:
            if line.startswith(">"):
                break
            rows.append(line.replace("\x00", "").strip())
    query = "".join(rows)
    if not query:
        raise ValueError(f"{path} has an empty query record")
    return query


def scan(path: Path) -> tuple[int, int]:
    """Count records and distinct taxa without holding the file in memory."""
    depth = 0
    taxa: set[str] = set()
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith(">"):
                depth += 1
                taxon = species_id(line)
                if taxon:
                    taxa.add(taxon)
    return depth, len(taxa)


def parse_a3m(text: str) -> Alignment:
    """Read an a3m string. Null bytes, used by the MSA services to separate
    per-query blocks, are stripped."""
    headers: list[str] = []
    rows: list[str] = []
    current: list[str] = []
    for line in text.replace("\x00", "").splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if headers:
                rows.append("".join(current))
            headers.append(line[1:])
            current = []
        elif headers:
            current.append(line.strip())
    if headers:
        rows.append("".join(current))
    if not headers:
        raise ValueError("no records found in a3m input")
    return Alignment(headers=headers, rows=rows)


def _group_by_species(alignment: Alignment) -> dict[str, list[int]]:
    """Row indices per taxon, keeping the a3m order.

    Hits arrive sorted by alignment score, so a taxon's first entry is its best
    match to that chain's query. Preserving the order is what makes "rank within
    a taxon" work without re-scoring anything.
    """
    groups: dict[str, list[int]] = {}
    for index, taxon in enumerate(alignment.species):
        if index == 0 or not taxon:
            continue
        groups.setdefault(taxon, []).append(index)
    return groups


def pair_alignments(
    alignments: list[Alignment],
    *,
    max_per_species: int = 1,
    require_all_chains: bool = False,
    max_rows: int | None = None,
) -> list[list[int | None]]:
    """Match rows across chains by taxon.

    Returns one entry per paired row, holding the row index to take from each
    chain, or None where that chain has no hit for the taxon. The first entry is
    always the queries themselves.

    ``max_per_species`` caps how many rows a single taxon may contribute, taking
    the best-scoring ones from each chain. ``require_all_chains`` switches from
    the greedy behaviour, which keeps a taxon as soon as two chains have it, to
    the stricter one that keeps only taxa present in every chain.
    """
    if len(alignments) < 2:
        raise ValueError("pairing needs at least two chains")

    per_chain = [_group_by_species(a) for a in alignments]
    paired: list[list[int | None]] = [[0] * len(alignments)]

    counts: dict[str, int] = {}
    for groups in per_chain:
        for taxon in groups:
            counts[taxon] = counts.get(taxon, 0) + 1

    threshold = len(alignments) if require_all_chains else 2
    candidates = [t for t, n in counts.items() if n >= threshold]
    # Taxa found in more chains carry more cross-chain signal, so they are
    # emitted first and survive the max_rows cut.
    candidates.sort(key=lambda t: (-counts[t], t))

    for taxon in candidates:
        available = [len(groups.get(taxon, ())) for groups in per_chain]
        depth = min(max_per_species, max(n for n in available if n) if any(available) else 0)
        for rank in range(depth):
            row: list[int | None] = []
            for groups in per_chain:
                rows = groups.get(taxon, ())
                row.append(rows[rank] if rank < len(rows) else None)
            paired.append(row)
            if max_rows is not None and len(paired) >= max_rows:
                return paired
    return paired


def render_paired(alignments: list[Alignment], pairing: list[list[int | None]]) -> list[str]:
    """Build concatenated rows from a pairing table, gapping absent chains."""
    out: list[str] = []
    for entry in pairing:
        parts = [
            alignment.rows[index] if index is not None else alignment.gap_row()
            for alignment, index in zip(alignments, entry, strict=True)
        ]
        out.append("".join(parts))
    return out


def render_block_diagonal(alignments: list[Alignment], skip_query: bool = True) -> list[str]:
    """Build the unpaired block: each chain's hits, gapped across every other."""
    out: list[str] = []
    for position, alignment in enumerate(alignments):
        start = 1 if skip_query else 0
        for index in range(start, alignment.depth):
            parts = [
                alignment.rows[index] if other == position else alignments[other].gap_row()
                for other in range(len(alignments))
            ]
            out.append("".join(parts))
    return out
