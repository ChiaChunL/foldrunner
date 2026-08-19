"""ColabFold writer.

ColabFold joins the chains of a complex with ``:`` inside one sequence, the
opposite of the multi-record layout AlphaFold2-Multimer expects from the same
file extension.

Precomputed alignments go in through the CSV form of the input, which carries an
``a3mpath`` column. They are written to a subdirectory rather than beside the
FASTA files: ``colabfold_batch`` accepts a directory and reads every ``.fasta``
*and* ``.a3m`` in it, so a flat layout makes it enumerate each complex twice.

The alignment itself uses the concatenated format ColabFold expects: a header
line giving each unique chain's length and copy count, then the paired rows,
then each chain's own hits padded across the others.
"""

from __future__ import annotations

import csv
from pathlib import Path

from foldrunner.engines.base import MsaView, WrittenJob, assign_chains, register
from foldrunner.ir import Job
from foldrunner.msa.a3m import pair_alignments, parse_a3m, render_block_diagonal, render_paired

NAME = "colabfold"

# Consumes a directory, so one invocation covers the whole panel;
# emitting one per complex would run every complex once per complex.
RUN_SCOPE = "directory"

COMMAND = "{executable} {input_dir}/queries.csv {out} {extra_args}"
DEFAULT_EXECUTABLE = "colabfold_batch"
# The CSV form is what carries the precomputed alignments.

CHAIN_JOIN = ":"

# Alignments live here rather than beside the FASTA files, so that pointing
# colabfold_batch at the engine directory yields one query per complex.
MSA_SUBDIR = "msas"


def query_sequence(job: Job) -> str:
    """The chains of the complex joined the way ColabFold reads them."""
    parts = []
    for component in job.components:
        parts.extend([component.entity.sequence] * component.count)
    return CHAIN_JOIN.join(parts)


def build_a3m(job: Job, msa: MsaView) -> str | None:
    """Concatenated alignment for one complex.

    The first line is ``#`` followed by the unique chain lengths and their copy
    counts, separated by a tab. ColabFold reads the query back by slicing the
    third line at those lengths, so the numbers have to describe the unique
    entities rather than the expanded chains.
    """
    proteins = [c.entity for c in job.components if c.is_protein]
    counts = [c.count for c in job.components if c.is_protein]
    if not proteins:
        return None

    alignments = []
    for protein in proteins:
        path = msa.paired(protein) or msa.unpaired(protein)
        if path is None:
            return None
        alignments.append(parse_a3m(path.read_text()))

    lengths = ",".join(str(len(p.sequence)) for p in proteins)
    cardinality = ",".join(str(c) for c in counts)
    header = f"#{lengths}\t{cardinality}\n"

    query = "".join(a.query for a in alignments)
    rows = [f">101\n{query}\n"]

    if len(alignments) > 1:
        paired = render_paired(alignments, pair_alignments(alignments))
        rows += [f">{101 + index}\n{row}\n" for index, row in enumerate(paired[1:], start=1)]
        for index, row in enumerate(render_block_diagonal(alignments)):
            rows.append(f">unpaired_{index}\n{row}\n")
    else:
        alignment = alignments[0]
        for index in range(1, alignment.depth):
            rows.append(f">{alignment.headers[index]}\n{alignment.rows[index]}\n")
    return header + "".join(rows)


@register(NAME)
def write(job: Job, outdir: Path, msa: MsaView | None = None) -> WrittenJob:
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    chains = assign_chains(job)

    a3m_path = None
    if msa is not None:
        text = build_a3m(job, msa)
        if text is not None:
            msa_dir = directory / MSA_SUBDIR
            msa_dir.mkdir(parents=True, exist_ok=True)
            a3m_path = msa_dir / f"{job.name}.a3m"
            a3m_path.write_text(text)

    path = directory / f"{job.name}.fasta"
    path.write_text(f">{job.name}\n{query_sequence(job)}\n")
    extra = {"a3m": str(a3m_path)} if a3m_path else {}
    return WrittenJob(
        job=job, engine=NAME, path=path, chains=chains, msa_used=a3m_path is not None, extra=extra
    )


def write_queries_csv(written: list[WrittenJob], outdir: Path, stem: str = "queries") -> Path:
    """Index the panel as a CSV, which is how ColabFold takes precomputed a3m files."""
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{stem}.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "sequence", "a3mpath"])
        for item in written:
            writer.writerow(
                [item.job.name, query_sequence(item.job), item.extra.get("a3m", "")]
            )
    return path
