"""Boltz-2 YAML writer.

Boltz-2 keeps the chain identifiers it is given and expresses copies as a list of
ids on one entity, so a homodimer is written once with ``id: [A, B]`` rather than
as two entries sharing a sequence.

Jobs carrying an affinity property are written to a separate directory. Boltz-2
runs the affinity stage over every entry in an input directory, so an entry
without the property fails there looking for a ``pre_affinity_*.npz`` that was
never produced; keeping the two kinds apart is the only way both run.
"""

from __future__ import annotations

from pathlib import Path

import yaml

from foldrunner.engines.base import MsaView, WrittenJob, assign_chains, register
from foldrunner.ir import Job, Protein

NAME = "boltz2"

# Consumes a directory, so one invocation covers the whole panel;
# emitting one per complex would run every complex once per complex.
RUN_SCOPE = "directory"

COMMAND = "boltz predict {input_dir} --out_dir {out} --seed {seed} {extra_args}"
# Boltz-2 takes a directory and one seed per run; the file name carries no
# seed, so several seeds mean several output directories.

# Boltz-2 reads an unset "msa" as "search for me" and the literal "empty" as
# single-sequence mode. Neither is wanted for a panel: the first re-searches the
# same sequence once per pair, the second throws the alignment away.
SINGLE_SEQUENCE = "empty"


def _protein_entry(protein: Protein, ids: list[str], msa: MsaView | None) -> dict:
    entry: dict[str, object] = {
        "id": ids if len(ids) > 1 else ids[0],
        "sequence": protein.sequence,
    }
    if msa is not None:
        path = msa.unpaired(protein)
        if path is not None:
            entry["msa"] = str(path)
    return {"protein": entry}


def _ligand_entry(ligand, ids: list[str]) -> dict:
    entry: dict[str, object] = {"id": ids if len(ids) > 1 else ids[0]}
    if ligand.smiles:
        entry["smiles"] = ligand.smiles
    else:
        entry["ccd"] = ligand.ccd
    return {"ligand": entry}


def build(job: Job, msa: MsaView | None = None) -> tuple[dict, list[list[str]]]:
    """Assemble the YAML document and the chain identifiers used."""
    chains = assign_chains(job)
    sequences = []
    for component, ids in zip(job.components, chains, strict=True):
        if isinstance(component.entity, Protein):
            sequences.append(_protein_entry(component.entity, ids, msa))
        else:
            sequences.append(_ligand_entry(component.entity, ids))

    document: dict[str, object] = {"version": 1, "sequences": sequences}

    if job.affinity_binder:
        binder_ids = [
            ids[0]
            for component, ids in zip(job.components, chains, strict=True)
            if component.entity.name == job.affinity_binder
        ]
        document["properties"] = [{"affinity": {"binder": binder_ids[0]}}]
    return document, chains


@register(NAME)
def write(job: Job, outdir: Path, msa: MsaView | None = None) -> WrittenJob:
    document, chains = build(job, msa)
    # Affinity jobs go to their own directory; see the module docstring.
    directory = Path(outdir) / ("affinity" if job.affinity_binder else "")
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{job.name}.yaml"
    path.write_text(yaml.safe_dump(document, sort_keys=False, default_flow_style=False))
    return WrittenJob(
        job=job,
        engine=NAME,
        path=path,
        chains=chains,
        msa_used=any("msa" in entry.get("protein", {}) for entry in document["sequences"]),
    )
