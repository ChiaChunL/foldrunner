"""Protenix JSON writer.

The input is always a list, even for a single complex, and copies are given as a
count rather than as named chains. Chain names in the output are assigned by
Protenix itself, always starting at A, so whatever is written here is advisory
only; the manifest records the mapping that will actually come back.

Alignments go in as paths on each protein chain. The older
``msa: {precomputed_msa_dir, pairing_db}`` form is deprecated and is not written.
"""

from __future__ import annotations

import json
from pathlib import Path

from foldrunner.engines.base import MsaView, WrittenJob, assign_chains, register
from foldrunner.ir import Job, Protein

NAME = "protenix"

# Takes the whole panel in a single file, in one invocation.
RUN_SCOPE = "file"

COMMAND = (
    "protenix pred -i {input} -o {out} -s {seeds} "
    "--need_atom_confidence True {extra_args}"
)
# Protenix accepts the whole panel and every seed in one invocation.
#
# --need_atom_confidence is not optional for a panel, whatever its name
# suggests. Without it Protenix writes only the summary confidences and the
# structure; the full_data file carrying token_pair_pae is never produced, so
# every interface metric that reads PAE — ipSAE, pDockQ2, LIS — is unavailable
# for every complex, and cannot be recovered without re-running the model.


def _entity(component, msa: MsaView | None) -> dict:
    entity = component.entity
    if isinstance(entity, Protein):
        chain: dict[str, object] = {"sequence": entity.sequence, "count": component.count}
        if msa is not None:
            unpaired = msa.unpaired(entity)
            paired = msa.paired(entity)
            if unpaired is not None:
                chain["unpairedMsaPath"] = str(unpaired)
            if paired is not None:
                chain["pairedMsaPath"] = str(paired)
        return {"proteinChain": chain}
    if entity.smiles:
        return {"ligand": {"ligand": entity.smiles, "count": component.count}}
    return {"ligand": {"ligand": f"CCD_{entity.ccd}", "count": component.count}}


def build(job: Job, msa: MsaView | None = None) -> tuple[dict, list]:
    chains = assign_chains(job)
    return (
        {"name": job.name, "sequences": [_entity(c, msa) for c in job.components]},
        chains,
    )


def write_batch(
    jobs: list[Job], outdir: Path, msa: MsaView | None = None, stem: str = "panel"
) -> list[WrittenJob]:
    """Write many complexes into one file.

    ``protenix pred`` takes a whole panel in a single input and several seeds in
    one invocation, so batching avoids re-loading the model per complex.
    """
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{stem}.json"
    entries = []
    written = []
    for job in jobs:
        document, chains = build(job, msa)
        entries.append(document)
        used = any("unpairedMsaPath" in e.get("proteinChain", {}) for e in document["sequences"])
        written.append(
            WrittenJob(job=job, engine=NAME, path=path, chains=chains, msa_used=used)
        )
    path.write_text(json.dumps(entries, indent=2))
    return written


@register(NAME)
def write(job: Job, outdir: Path, msa: MsaView | None = None) -> WrittenJob:
    document, chains = build(job, msa)
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{job.name}.json"
    path.write_text(json.dumps([document], indent=2))
    used = any("unpairedMsaPath" in e.get("proteinChain", {}) for e in document["sequences"])
    return WrittenJob(job=job, engine=NAME, path=path, chains=chains, msa_used=used)
