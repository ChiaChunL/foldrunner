"""Resuming a run, and the partial alignment reuse AlphaFold2 allows."""

from foldrunner.cli import main
from foldrunner.engines import get_engine
from foldrunner.engines.af2_multimer import BFD_A3M, STOCKHOLM_SLOTS, msa_directory
from foldrunner.engines.base import MsaView
from foldrunner.ir import Component, Job, Protein
from foldrunner.manifest import completed_jobs, read_manifest
from foldrunner.msa.cache import MsaCache, Recipe

A = Protein("A", "AQVINTFDGV")
B = Protein("B", "KKAVINGEQI")


def build_cache(tmp_path):
    cache = MsaCache(tmp_path / "cache")
    recipe = Recipe(backend="test", databases=("uniref30",), tool_version="18")
    for protein in (A, B):
        cache.store(
            recipe,
            protein.key,
            protein.sequence,
            unpaired=f">query\n{protein.sequence}\n>ERR1_2\n{protein.sequence}\n",
            paired=f">query\n{protein.sequence}\n>UniRef100_X_9606/\n{protein.sequence}\n",
        )
    return MsaView(cache=cache, recipe=recipe)


def write_library(tmp_path):
    fasta = tmp_path / "lib.fasta"
    fasta.write_text(f">A\n{A.sequence}\n>B\n{B.sequence}\n")
    return fasta


# --- resume -----------------------------------------------------------------

def test_resume_skips_what_the_manifest_records(tmp_path):
    fasta, out = write_library(tmp_path), tmp_path / "out"
    assert main(["write", str(fasta), "-o", str(out), "--engines", "boltz2"]) == 0
    first = {p.name for p in (out / "boltz2").glob("*.yaml")}
    assert completed_jobs(out / "manifest.tsv", "boltz2") == {"A__A", "A__B", "B__B"}

    for path in (out / "boltz2").glob("*.yaml"):
        path.unlink()
    assert main(["write", str(fasta), "-o", str(out), "--engines", "boltz2", "--resume"]) == 0
    # Nothing was rewritten, because the manifest already recorded all three.
    assert {p.name for p in (out / "boltz2").glob("*.yaml")} == set()
    assert first  # the first run really did produce files


def test_resume_writes_only_what_is_missing(tmp_path):
    fasta, out = write_library(tmp_path), tmp_path / "out"
    main(["write", str(fasta), "-o", str(out), "--engines", "boltz2"])

    # Drop one job from the manifest, as an interrupted run would leave it.
    rows = read_manifest(out / "manifest.tsv")
    kept = [r for r in rows if r["job"] != "A__B"]
    header = "\t".join(rows[0].keys())
    body = "\n".join("\t".join(r[k] for k in rows[0]) for r in kept)
    (out / "manifest.tsv").write_text(f"{header}\n{body}\n")
    for path in (out / "boltz2").glob("*.yaml"):
        path.unlink()

    main(["write", str(fasta), "-o", str(out), "--engines", "boltz2", "--resume"])
    assert {p.stem for p in (out / "boltz2").glob("*.yaml")} == {"A__B"}


def test_resume_keeps_the_earlier_rows_in_the_manifest(tmp_path):
    fasta, out = write_library(tmp_path), tmp_path / "out"
    main(["write", str(fasta), "-o", str(out), "--engines", "boltz2"])
    before = len(read_manifest(out / "manifest.tsv"))
    main(["write", str(fasta), "-o", str(out), "--engines", "boltz2", "--resume"])
    assert len(read_manifest(out / "manifest.tsv")) == before


def test_a_fresh_run_without_resume_is_unaffected(tmp_path):
    fasta, out = write_library(tmp_path), tmp_path / "out"
    main(["write", str(fasta), "-o", str(out), "--engines", "boltz2"])
    main(["write", str(fasta), "-o", str(out), "--engines", "boltz2"])
    assert len(read_manifest(out / "manifest.tsv")) == 6


# --- AlphaFold2 partial reuse ----------------------------------------------

def test_af2_alignment_directories_are_indexed_from_a_over_unique_sequences(tmp_path):
    """Verified against a real run: msas/ is indexed from A over distinct
    sequences, while the chains in the predicted structure start at B."""
    msa = build_cache(tmp_path)
    job = Job("A__B", (Component(A), Component(B)))
    written = get_engine("af2_multimer")(job, tmp_path / "af2", msa)
    root = msa_directory(tmp_path / "af2", job)
    assert written.msa_used
    assert sorted(p.name for p in root.iterdir()) == ["A", "B"]
    assert (root / "A" / BFD_A3M).is_file()


def test_af2_homodimer_gets_one_alignment_directory_not_two(tmp_path):
    """AlphaFold2 searches once per sequence, so copies share a directory."""
    msa = build_cache(tmp_path)
    job = Job("A__A", (Component(A, 2),))
    get_engine("af2_multimer")(job, tmp_path / "af2", msa)
    root = msa_directory(tmp_path / "af2", job)
    assert [p.name for p in root.iterdir()] == ["A"]


def test_af2_structure_chains_still_start_at_b(tmp_path):
    job = Job("A__A", (Component(A, 2),))
    written = get_engine("af2_multimer")(job, tmp_path / "af2")
    assert written.chains == [["B", "C"]]


def test_af2_records_which_searches_it_cannot_supply(tmp_path):
    """The Stockholm slots are still searched by AlphaFold2; say so rather than
    imply the alignment was fully reused."""
    msa = build_cache(tmp_path)
    job = Job("A__B", (Component(A), Component(B)))
    written = get_engine("af2_multimer")(job, tmp_path / "af2", msa)
    assert written.extra["msa_partial"] == ",".join(STOCKHOLM_SLOTS)
    assert "uniprot_hits.sto" in written.extra["msa_partial"]


def test_af2_without_a_cache_writes_no_msa_directory(tmp_path):
    job = Job("A__B", (Component(A), Component(B)))
    written = get_engine("af2_multimer")(job, tmp_path / "af2", None)
    assert not written.msa_used
    assert not msa_directory(tmp_path / "af2", job).exists()


def test_af2_command_stages_alignments_into_the_output_tree():
    """--use-precomputed-msas reads from <outdir>/<name>/msas/, not from beside
    the input, so the alignments have to be copied there before the run."""
    from foldrunner.engines.af2_multimer import COMMAND

    assert "{msa_dir}" in COMMAND
    assert "{out}/{job}" in COMMAND
    assert COMMAND.index("cp -rn") < COMMAND.index("{executable}")
