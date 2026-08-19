"""AlphaFold2-Multimer FASTA writer.

One file per complex with one record per chain, which is the opposite of what
ColabFold reads from the same extension.

Chain identifiers start at B rather than A, matching what the pipeline writes
into its own output, so results can be matched back without an off-by-one.
There is no seed setting; repeated sampling comes from
``--num-multimer-predictions-per-model``, which draws several samples from one
alignment rather than re-running with a different seed.

Alignment reuse here is **partial**, and deliberately so. ``--use-precomputed-msas``
is decided per file: for each search the pipeline checks whether its output
already exists and only runs the search when it does not. Of the files it looks
for, exactly one is a3m — ``bfd_uniref_hits.a3m`` — while ``uniref90_hits.sto``,
``mgnify_hits.sto`` and the ``uniprot_hits.sto`` used for pairing are Stockholm.
Only the a3m is written; the rest are left for AlphaFold2 to search, because
inventing Stockholm records whose identifiers no longer support species pairing
would quietly change the alignment rather than reuse it.
"""

from __future__ import annotations

from pathlib import Path

from foldrunner.engines.base import (
    MsaView,
    WrittenJob,
    assign_chains,
    chain_labels,
    register,
)
from foldrunner.ir import Job

NAME = "af2_multimer"

# One invocation per complex.
RUN_SCOPE = "job"

COMMAND = (
    "mkdir -p {out}/{job} && "
    "if [ -d {msa_dir} ]; then cp -rn {msa_dir} {out}/{job}/ 2>/dev/null || true; fi && "
    "{executable} --fasta {input} --outdir {out} {extra_args}"
)
DEFAULT_EXECUTABLE = "run_af2.sh"
# AlphaFold2 is normally driven by a site wrapper script; point executable at
# it in the config.
#
# The multimer preset has to be requested explicitly. Both AlphaFold2 itself
# and the usual wrappers default to the monomer pipeline, which rejects a
# multi-record FASTA with "More than one input sequence found" — after the
# job has already started. The flag spelling differs between entry points
# (--model_preset=multimer for run_alphafold.py, --model-preset multimer for
# a typical wrapper), so it belongs in extra_args rather than here.
#
# The alignments have to be staged into the *output* tree before the run:
# --use-precomputed-msas looks for them under <outdir>/<name>/msas/<chain>/,
# not next to the input. Copying rather than moving keeps the panel
# reproducible, and -n leaves anything the pipeline already wrote alone.

# The pipeline labels the first chain B, not A.
FIRST_CHAIN = "B"

# The one file in AlphaFold2's precomputed set that is already a3m.
BFD_A3M = "bfd_uniref_hits.a3m"

# Files the pipeline will still search for itself; naming them keeps the
# limitation visible rather than implied.
STOCKHOLM_SLOTS = ("uniref90_hits.sto", "mgnify_hits.sto", "uniprot_hits.sto")


def msa_directory(outdir: Path, job: Job) -> Path:
    """Where AlphaFold2 looks for a target's alignments."""
    return Path(outdir) / job.name / "msas"


def write_precomputed_msas(job: Job, outdir: Path, msa: MsaView) -> list[Path]:
    """Place the a3m AlphaFold2 can reuse, one directory per unique sequence.

    AlphaFold2 uses two different chain namings in the same output, and they are
    easy to conflate. The directories under ``msas/`` are indexed from **A** over
    the *distinct* sequences, because the pipeline searches once per sequence
    rather than once per copy — a homodimer has one directory, not two. The
    chains written into the predicted structure start at **B** instead, and
    ``msas/chain_id_map.json`` is what ties the two together.
    """
    written: list[Path] = []
    root = msa_directory(outdir, job)
    seen: dict[str, str] = {}
    labels = chain_labels()
    for component in job.components:
        if not component.is_protein:
            continue
        key = component.entity.key
        if key in seen:
            continue
        seen[key] = next(labels)
        source = msa.unpaired(component.entity)
        if source is None:
            continue
        target = root / seen[key] / BFD_A3M
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(source.read_text())
        written.append(target)
    return written


@register(NAME)
def write(job: Job, outdir: Path, msa: MsaView | None = None) -> WrittenJob:
    chains = assign_chains(job, start=FIRST_CHAIN)
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)

    records = []
    for component, ids in zip(job.components, chains, strict=True):
        if not component.is_protein:
            continue
        for label in ids:
            records.append(f">{job.name}_{label}\n{component.entity.sequence}\n")

    path = directory / f"{job.name}.fasta"
    path.write_text("".join(records))

    extra: dict[str, str] = {}
    used = False
    if msa is not None and write_precomputed_msas(job, directory, msa):
        used = True
        extra["msa_directory"] = str(msa_directory(directory, job))
        extra["msa_partial"] = ",".join(STOCKHOLM_SLOTS)
    return WrittenJob(job=job, engine=NAME, path=path, chains=chains, msa_used=used, extra=extra)
