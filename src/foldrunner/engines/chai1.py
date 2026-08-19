"""Chai-1 FASTA writer.

Chai-1 reads one record per chain with the entity type in the header, so copies
are written out as repeated records rather than as a count.

Alignments are supplied through a directory of ``<sha256>.aligned.pqt`` files,
one per unique chain sequence, where the hash is taken over the uppercased
sequence — the same address the cache already uses.

The parquet is **not** written here. Parquet is only compatible in one
direction: a file written by a newer pyarrow fails to load on an older one with
``Repetition level histogram size mismatch``, and the environment that reads it
is the engine's, whose versions this package does not control. So the a3m files
are written instead, and the generated driver script converts them with Chai-1's
own function while running inside Chai-1's own environment.
"""

from __future__ import annotations

from pathlib import Path

from foldrunner.engines.base import (
    EngineDependencyError,
    MsaView,
    WrittenJob,
    assign_chains,
    register,
)
from foldrunner.ir import Job, Protein
from foldrunner.msa.a3m import Alignment, species_id

NAME = "chai1"

# One invocation per complex.
RUN_SCOPE = "job"

COMMAND = "python {script} {input} {out} {msa_dir} {msa_a3m_dir}"
# Driven through a generated Python script rather than `chai-lab fold`:
# the CLI never writes the full PAE to disk, and it cannot be recovered
# afterwards. See write_runner_script below.

# Columns Chai-1 expects in an aligned.pqt. "pairing_key" is what it matches on
# when merging alignments across chains; rows left blank there are never paired.
PQT_COLUMNS = ("sequence", "source_database", "pairing_key", "comment")

def _pandas():
    """Import pandas, naming the extra that provides it.

    Only this writer needs it, so it is optional; a bare ImportError here would
    leave the reader guessing which dependency to install.
    """
    try:
        import pandas
    except ImportError as error:  # pragma: no cover - environment dependent
        raise EngineDependencyError(
            "Chai-1 alignments are parquet files, which needs pandas and pyarrow. "
            'Install them with: pip install "foldrunner[chai]"'
        ) from error
    return pandas


QUERY_SOURCE = "query"
PAIRED_SOURCE = "uniprot"
UNPAIRED_SOURCE = "uniref90"


def _record(entity, label: str = "") -> str:
    """One FASTA record. ``label`` disambiguates copies of one entity.

    Chai-1 rejects a file where two records share a name, so the copies of a
    homomer cannot simply repeat it. The chain identifier is used because it
    is the same one recorded in the manifest, which keeps the results
    traceable back to the entity.
    """
    name = f"{entity.name}_{label}" if label else entity.name
    if isinstance(entity, Protein):
        return f">protein|name={name}\n{entity.sequence}\n"
    if entity.smiles:
        return f">ligand|name={name}\n{entity.smiles}\n"
    return f">ligand|name={name}\n{entity.ccd}\n"


def build(job: Job) -> tuple[str, list[list[str]]]:
    chains = assign_chains(job)
    parts = []
    for component, labels in zip(job.components, chains, strict=True):
        # A single copy keeps the plain entity name; only copies need a suffix.
        if component.count == 1:
            parts.append(_record(component.entity))
        else:
            parts.extend(_record(component.entity, label) for label in labels)
    return "".join(parts), chains


def aligned_frame(paired: Alignment | None, unpaired: Alignment | None):
    """Build the dataframe behind an ``aligned.pqt``.

    Rows from the pairing source keep their taxon in ``pairing_key`` so Chai-1
    can match them across chains; rows from unpaired sources leave it empty and
    contribute only within their own chain.
    """
    pd = _pandas()

    rows: list[dict[str, str]] = []
    query = (paired or unpaired)
    if query is None:
        raise ValueError("need at least one alignment to build an aligned.pqt")
    rows.append(
        {
            "sequence": query.rows[0],
            "source_database": QUERY_SOURCE,
            "pairing_key": "",
            "comment": "query",
        }
    )
    if paired is not None:
        for header, row in zip(paired.headers[1:], paired.rows[1:], strict=True):
            taxon = species_id(header)
            rows.append(
                {
                    "sequence": row,
                    "source_database": PAIRED_SOURCE,
                    "pairing_key": taxon or "",
                    "comment": header.split()[0] if header.split() else "",
                }
            )
    if unpaired is not None:
        for header, row in zip(unpaired.headers[1:], unpaired.rows[1:], strict=True):
            rows.append(
                {
                    "sequence": row,
                    "source_database": UNPAIRED_SOURCE,
                    "pairing_key": "",
                    "comment": header.split()[0] if header.split() else "",
                }
            )
    return pd.DataFrame(rows, columns=list(PQT_COLUMNS))


def write_source_a3m(protein: Protein, msa: MsaView, directory: Path) -> Path | None:
    """Write one sequence's alignments as a3m, one directory per sequence.

    The file names name their source, which is what the conversion maps onto
    Chai-1's own database labels.
    """
    written = False
    target = Path(directory) / protein.key
    sources = (
        (msa.paired, f"{PAIRED_SOURCE}.a3m"),
        (msa.unpaired, f"{UNPAIRED_SOURCE}.a3m"),
    )
    for kind, name in sources:
        source = kind(protein)
        if source is None:
            continue
        target.mkdir(parents=True, exist_ok=True)
        (target / name).write_text(source.read_text())
        written = True
    return target if written else None


def write_aligned_pqt(protein: Protein, msa: MsaView, directory: Path) -> Path | None:
    """Materialise one sequence's alignments as ``<sha256>.aligned.pqt``.

    Kept for callers whose reader and writer share an environment; the runner
    path uses :func:`write_source_a3m` instead, for the reason in the module
    docstring.
    """
    paired_path = msa.paired(protein)
    unpaired_path = msa.unpaired(protein)
    if paired_path is None and unpaired_path is None:
        return None

    from foldrunner.msa.a3m import parse_a3m

    paired = parse_a3m(paired_path.read_text()) if paired_path else None
    unpaired = parse_a3m(unpaired_path.read_text()) if unpaired_path else None
    frame = aligned_frame(paired, unpaired)
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f"{protein.key}.aligned.pqt"
    frame.to_parquet(target, index=False)
    return target


@register(NAME)
def write(job: Job, outdir: Path, msa: MsaView | None = None) -> WrittenJob:
    text, chains = build(job)
    directory = Path(outdir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{job.name}.fasta"
    path.write_text(text)

    msa_used = False
    extra: dict[str, str] = {"script": str(write_runner_script(directory))}
    if msa is not None:
        a3m_dir = directory / "msas_a3m"
        written = [write_source_a3m(p, msa, a3m_dir) for p in job.proteins]
        if any(w is not None for w in written):
            msa_used = True
            extra["msa_directory"] = str(directory / "msas")
            extra["msa_a3m_directory"] = str(a3m_dir)
    return WrittenJob(
        job=job, engine=NAME, path=path, chains=chains, msa_used=msa_used, extra=extra
    )


# Chai-1 is driven through this script rather than `chai-lab fold`. The CLI
# writes only the ranking summary; the full PAE exists solely on the object
# returned by run_inference, so a CLI run loses it permanently.
RUNNER_SCRIPT = '''\
"""Run Chai-1, converting alignments in place and keeping the full PAE.

Both steps have to happen here rather than upstream. The parquet Chai-1 reads is
only compatible in one direction, so it is built with the pyarrow that lives in
this environment. And the full PAE exists only on the object run_inference
returns — the command line writes the ranking summary and nothing else, so a CLI
run loses it for good.
"""

import sys
from pathlib import Path

import numpy as np
import torch
from chai_lab.chai1 import run_inference

fasta, out, msa_dir = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
a3m_dir = Path(sys.argv[4]) if len(sys.argv) > 4 else None
out.mkdir(parents=True, exist_ok=True)

if a3m_dir is not None and a3m_dir.is_dir():
    from chai_lab.data.parsing.msas.aligned_pqt import (
        expected_basename,
        merge_multi_a3m_to_aligned_dataframe,
    )
    from chai_lab.data.parsing.msas.data_source import MSADataSource

    SOURCES = {
        "uniprot.a3m": MSADataSource.UNIPROT,
        "uniref90.a3m": MSADataSource.UNIREF90,
    }
    target = Path(msa_dir)
    target.mkdir(parents=True, exist_ok=True)
    for entry in sorted(a3m_dir.iterdir()):
        if not entry.is_dir():
            continue
        files = {entry / n: s for n, s in SOURCES.items() if (entry / n).is_file()}
        if not files:
            continue
        frame = merge_multi_a3m_to_aligned_dataframe(files)
        query = frame.iloc[0]["sequence"].replace("-", "").upper()
        frame.to_parquet(target / expected_basename(query), index=False)
    print(f"converted {len(list(target.glob('*.aligned.pqt')))} alignments")

candidates = run_inference(
    fasta_file=fasta,
    output_dir=out,
    msa_directory=Path(msa_dir) if msa_dir not in ("", "-") else None,
    device=torch.device("cuda:0" if torch.cuda.is_available() else "cpu"),
    use_esm_embeddings=True,
)

pae = getattr(candidates, "pae", None)
if pae is not None:
    array = pae.detach().cpu().numpy() if hasattr(pae, "detach") else np.asarray(pae)
    np.savez_compressed(out / "pae.npz", pae=array)
    print("wrote", out / "pae.npz", array.shape)
else:
    print("warning: this chai-lab build returned no pae attribute", file=sys.stderr)
'''


def write_runner_script(directory: Path) -> Path:
    """Place the driver script next to the inputs."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "run_chai.py"
    path.write_text(RUNNER_SCRIPT)
    return path
