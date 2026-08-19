"""SeedFold job writer.

Entities carry a ``copies`` count, and every kind of content — sequence, CCD
code or SMILES — travels in the same ``sequence`` field, distinguished only by
the entity label.

This is the one engine that cannot be given a precomputed alignment or a seed:
the web service builds its own alignment and its sampling is not reproducible.
Any comparison that includes it has to say so, or a difference in how the
alignment was built will be read as a difference between models.
"""

from __future__ import annotations

import json
from pathlib import Path

from foldrunner.engines.base import MsaView, WrittenJob, assign_chains, register
from foldrunner.ir import Job, Protein

NAME = "seedfold"

# Uploaded by hand; nothing to invoke.
RUN_SCOPE = "manual"

# Submitted through a web form, so there is no command to run. The runner
# reports this rather than inventing one; SeedFold jobs are uploaded by hand.
COMMAND = None

DEFAULT_MODEL = "SeedFold_v1.0.0"
LINEAR_MODEL = "SeedFold-Linear_v1.0.0"

PROTEIN = "Protein"
LIGAND_CCD = "Ligand/Ion-CCD"
LIGAND_SMILES = "Ligand-Smiles"


def _entity(component) -> dict:
    entity = component.entity
    if isinstance(entity, Protein):
        return {"entity": PROTEIN, "copies": component.count, "sequence": entity.sequence}
    if entity.smiles:
        return {"entity": LIGAND_SMILES, "copies": component.count, "sequence": entity.smiles}
    return {"entity": LIGAND_CCD, "copies": component.count, "sequence": entity.ccd}


def build(job: Job, model: str = DEFAULT_MODEL) -> dict:
    """One batch entry. The single-prediction form is the entity list alone."""
    return {
        "job_name": job.name,
        "model": model,
        "entities": [_entity(c) for c in job.components],
    }


def write_batch(
    jobs: list[Job],
    outdir: Path,
    msa: MsaView | None = None,
    model: str = DEFAULT_MODEL,
    chunk: int = 100,
) -> list[WrittenJob]:
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    written: list[WrittenJob] = []
    for index in range(0, len(jobs), chunk):
        batch = jobs[index : index + chunk]
        path = directory / f"batch_{index // chunk:03d}.json"
        path.write_text(json.dumps([build(job, model) for job in batch], indent=2))
        for job in batch:
            written.append(
                WrittenJob(
                    job=job,
                    engine=NAME,
                    path=path,
                    chains=assign_chains(job),
                    msa_used=False,
                )
            )
    return written


@register(NAME)
def write(job: Job, outdir: Path, msa: MsaView | None = None) -> WrittenJob:
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{job.name}.json"
    # The single-prediction form is a bare list of entities, with no job wrapper.
    path.write_text(json.dumps([_entity(c) for c in job.components], indent=2))
    return WrittenJob(
        job=job, engine=NAME, path=path, chains=assign_chains(job), msa_used=False
    )
