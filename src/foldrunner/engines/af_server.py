"""AlphaFold Server job writer.

The web form differs from local AlphaFold 3 in ways that are easy to miss: the
top level is a list, the entity key is ``proteinChain`` rather than ``protein``,
seeds are strings, and copies are a count instead of named identifiers.

Two limits shape the output. Ligands may only be given as one of nineteen CCD
codes — custom SMILES are rejected outright — and an account may submit only a
few dozen jobs a day, so a panel is written as batches sized to a day's quota
rather than as one file nobody can upload.
"""

from __future__ import annotations

import json
from pathlib import Path

from foldrunner.engines.base import MsaView, WrittenJob, assign_chains, register
from foldrunner.ir import IRError, Job, Protein

NAME = "af_server"

# Uploaded by hand; nothing to invoke.
RUN_SCOPE = "manual"

# Submitted through a web form, so there is no command to run. The runner
# reports this rather than inventing one; AlphaFold Server jobs are uploaded by hand.
COMMAND = None
DIALECT = "alphafoldserver"

# The complete accepted ligand set. Anything else has to be dropped or the job
# is refused on upload.
LIGAND_ALLOWLIST = frozenset(
    """ADP ATP AMP GTP GDP FAD NAD NAP NDP HEM HEC PLM OLA MYR CIT CLA CHL BCL BCB""".split()
)
ION_ALLOWLIST = frozenset("MG ZN CL CA NA MN K FE CU CO".split())

# Academic accounts get a few dozen jobs a day; batches are sized so one file is
# one day's submission.
DEFAULT_DAILY_QUOTA = 20


def _protein(component, msa: MsaView | None, max_template_date) -> dict:
    entity = component.entity
    chain: dict[str, object] = {"sequence": entity.sequence, "count": component.count}
    if max_template_date is not None:
        chain["maxTemplateDate"] = max_template_date.isoformat()
    if msa is not None:
        path = msa.unpaired(entity)
        if path is not None:
            # The web form has no filesystem to point at, so the alignment has
            # to travel inside the job description.
            chain["unpairedMsa"] = path.read_text()
    return {"proteinChain": chain}


def _ligand(component) -> dict:
    entity = component.entity
    code = (entity.ccd or "").upper()
    if code in ION_ALLOWLIST:
        return {"ion": {"ion": code, "count": component.count}}
    if code in LIGAND_ALLOWLIST:
        return {"ligand": {"ligand": f"CCD_{code}", "count": component.count}}
    raise IRError(
        f"AlphaFold Server does not accept ligand {entity.name!r}: it takes only the "
        f"nineteen allowed CCD codes and ten ions, never custom SMILES"
    )


def build(job: Job, msa: MsaView | None = None, seed: int | None = None) -> tuple[dict, list]:
    chains = assign_chains(job)
    sequences = []
    for component in job.components:
        if isinstance(component.entity, Protein):
            sequences.append(_protein(component, msa, job.max_template_date))
        else:
            sequences.append(_ligand(component))

    name = job.name if seed is None else f"{job.name}_s{seed}"
    return (
        {
            "name": name,
            "modelSeeds": [str(seed)] if seed is not None else [],
            "sequences": sequences,
            "dialect": DIALECT,
            "version": 1,
        },
        chains,
    )


def write_batch(
    jobs: list[Job],
    outdir: Path,
    msa: MsaView | None = None,
    quota: int = DEFAULT_DAILY_QUOTA,
) -> list[WrittenJob]:
    """Write the panel as quota-sized upload batches.

    One job per seed: the published schema allows several, but the upload form
    has been seen to reject them, and one-per-job is accepted either way.
    """
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)

    entries: list[tuple[Job, dict, list]] = []
    for job in jobs:
        for seed in job.seeds or (None,):
            document, chains = build(job, msa, seed)
            entries.append((job, document, chains))

    written: list[WrittenJob] = []
    for index in range(0, len(entries), quota):
        chunk = entries[index : index + quota]
        path = directory / f"batch_{index // quota:03d}.json"
        path.write_text(json.dumps([document for _, document, _ in chunk], indent=2))
        for job, document, chains in chunk:
            used = any(
                "unpairedMsa" in e.get("proteinChain", {}) for e in document["sequences"]
            )
            written.append(
                WrittenJob(job=job, engine=NAME, path=path, chains=chains, msa_used=used)
            )
    return written


@register(NAME)
def write(job: Job, outdir: Path, msa: MsaView | None = None) -> WrittenJob:
    return write_batch([job], outdir, msa, quota=DEFAULT_DAILY_QUOTA)[0]
