"""AlphaFold 3 JSON writer.

Copies are expressed as a list of identifiers on one entity, and the identifiers
are kept in the output, so this is one of the two engines whose chain names can
be chosen up front.

Alignments are referenced by path rather than embedded. Both forms are accepted,
but a panel repeats the same alignment across every complex a sequence appears
in, and an alignment for a few-hundred-residue protein runs to tens of megabytes;
embedding turns a hundred jobs into gigabytes of duplicated text. Supplying them
also matters for a second reason: left to itself AlphaFold 3 runs its own data
pipeline once per job, searching the same sequence once per pair it appears in.
"""

from __future__ import annotations

import json
from pathlib import Path

from foldrunner.engines.base import MsaView, WrittenJob, assign_chains, register
from foldrunner.ir import Job, Protein

NAME = "af3"

# One invocation per complex.
RUN_SCOPE = "job"

COMMAND = (
    "python /app/alphafold/run_alphafold.py "
    "--json_path={input} --output_dir={out} {extra_args}"
)
# AlphaFold 3 needs --model_dir and --db_dir, which are site paths; set them
# through extra_args in the config rather than guessing them here.
DIALECT = "alphafold3"
DEFAULT_VERSION = 3


def _protein_entry(protein: Protein, ids: list[str], msa: MsaView | None) -> dict:
    entry: dict[str, object] = {
        "id": ids if len(ids) > 1 else ids[0],
        "sequence": protein.sequence,
    }
    if msa is not None:
        unpaired = msa.unpaired(protein)
        paired = msa.paired(protein)
        if unpaired is not None:
            entry["unpairedMsaPath"] = str(unpaired)
        if paired is not None:
            entry["pairedMsaPath"] = str(paired)
        if "unpairedMsaPath" in entry or "pairedMsaPath" in entry:
            # Supplying alignments means the data pipeline is skipped, and a
            # skipped pipeline cannot fill templates in either. Left unset the
            # run fails with "Protein chain N is missing Templates" only once
            # the model is already loaded. An empty list means "no templates",
            # which is what a panel wants: searching them per job would reopen
            # the per-job cost the alignments were cached to avoid.
            entry.setdefault("templates", [])
    return {"protein": entry}


def build(
    job: Job, msa: MsaView | None = None, version: int = DEFAULT_VERSION
) -> tuple[dict, list]:
    chains = assign_chains(job)
    sequences = []
    for component, ids in zip(job.components, chains, strict=True):
        if isinstance(component.entity, Protein):
            sequences.append(_protein_entry(component.entity, ids, msa))
        else:
            ligand: dict[str, object] = {"id": ids if len(ids) > 1 else ids[0]}
            if component.entity.smiles:
                ligand["smiles"] = component.entity.smiles
            else:
                ligand["ccdCodes"] = [component.entity.ccd]
            sequences.append({"ligand": ligand})

    document: dict[str, object] = {
        "name": job.name,
        "sequences": sequences,
        "modelSeeds": list(job.seeds) or [1],
        "dialect": DIALECT,
        "version": version,
    }
    return document, chains


@register(NAME)
def write(job: Job, outdir: Path, msa: MsaView | None = None) -> WrittenJob:
    document, chains = build(job, msa)
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{job.name}.json"
    path.write_text(json.dumps(document, indent=2))
    used = any("unpairedMsaPath" in e.get("protein", {}) for e in document["sequences"])
    return WrittenJob(job=job, engine=NAME, path=path, chains=chains, msa_used=used)
