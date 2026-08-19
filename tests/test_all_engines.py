"""Format checks for every engine, including the disagreements between them."""

import json

import pytest
import yaml

from foldrunner.engines import ENGINES, get_engine
from foldrunner.engines.af_server import LIGAND_ALLOWLIST
from foldrunner.engines.af_server import write_batch as af_server_batch
from foldrunner.engines.base import MsaView
from foldrunner.engines.colabfold import build_a3m, query_sequence, write_queries_csv
from foldrunner.engines.protenix import write_batch as protenix_batch
from foldrunner.engines.seedfold import write_batch as seedfold_batch
from foldrunner.ir import Component, IRError, Job, Ligand, Protein
from foldrunner.msa.cache import MsaCache, Recipe

A = Protein("A", "AQVINTFDGV")
B = Protein("B", "KKAVINGEQI")

PAIRED = ">query\n{seq}\n>UniRef100_X_9606/\t1\t2\t3\n{seq}\n"
UNPAIRED = ">query\n{seq}\n>ERR1_2\t1\t2\t3\n{seq}\n"


@pytest.fixture
def msa(tmp_path):
    cache = MsaCache(tmp_path / "cache")
    recipe = Recipe(backend="test", databases=("uniref30",), tool_version="18")
    for protein in (A, B):
        cache.store(
            recipe,
            protein.key,
            protein.sequence,
            unpaired=UNPAIRED.format(seq=protein.sequence),
            paired=PAIRED.format(seq=protein.sequence),
        )
    return MsaView(cache=cache, recipe=recipe)


@pytest.fixture
def hetero():
    return Job("A__B", (Component(A), Component(B)), seeds=(2066,))


@pytest.fixture
def homo():
    return Job("A__A", (Component(A, 2),), seeds=(2066,))


def test_every_engine_writes_something(tmp_path, hetero, msa):
    for name in ENGINES:
        written = get_engine(name)(hetero, tmp_path / name, msa)
        assert written.path.is_file(), name
        assert written.path.stat().st_size > 0, name


def test_all_engines_agree_on_the_job_name(tmp_path, hetero, msa):
    for name in ENGINES:
        assert get_engine(name)(hetero, tmp_path / name, msa).job.name == "A__B"


# --- AlphaFold 3 ------------------------------------------------------------

def test_af3_uses_its_own_dialect_and_id_lists(tmp_path, homo):
    document = json.loads(get_engine("af3")(homo, tmp_path).path.read_text())
    assert document["dialect"] == "alphafold3"
    assert document["sequences"][0]["protein"]["id"] == ["A", "B"]
    assert document["modelSeeds"] == [2066]


def test_af3_references_alignments_by_path_not_by_value(tmp_path, hetero, msa):
    document = json.loads(get_engine("af3")(hetero, tmp_path, msa).path.read_text())
    entry = document["sequences"][0]["protein"]
    assert entry["unpairedMsaPath"].endswith("unpaired.a3m")
    assert entry["pairedMsaPath"].endswith("paired.a3m")
    assert "unpairedMsa" not in entry


# --- Protenix ---------------------------------------------------------------

def test_protenix_input_is_always_a_list(tmp_path, hetero):
    document = json.loads(get_engine("protenix")(hetero, tmp_path).path.read_text())
    assert isinstance(document, list) and len(document) == 1


def test_protenix_uses_counts_and_current_msa_fields(tmp_path, homo, msa):
    entry = json.loads(get_engine("protenix")(homo, tmp_path, msa).path.read_text())[0]
    chain = entry["sequences"][0]["proteinChain"]
    assert chain["count"] == 2
    assert "unpairedMsaPath" in chain
    assert "msa" not in chain  # the precomputed_msa_dir form is deprecated


def test_protenix_batches_a_whole_panel_into_one_file(tmp_path, hetero, homo, msa):
    written = protenix_batch([hetero, homo], tmp_path, msa)
    assert len({w.path for w in written}) == 1
    assert len(json.loads(written[0].path.read_text())) == 2


# --- AlphaFold Server -------------------------------------------------------

def test_af_server_uses_the_web_dialect_and_string_seeds(tmp_path, homo, msa):
    batch = json.loads(get_engine("af_server")(homo, tmp_path, msa).path.read_text())
    assert isinstance(batch, list)
    entry = batch[0]
    assert entry["dialect"] == "alphafoldserver"
    assert entry["modelSeeds"] == ["2066"]
    assert entry["sequences"][0]["proteinChain"]["count"] == 2


def test_af_server_embeds_alignments_because_there_is_no_filesystem(tmp_path, hetero, msa):
    entry = json.loads(get_engine("af_server")(hetero, tmp_path, msa).path.read_text())[0]
    chain = entry["sequences"][0]["proteinChain"]
    assert chain["unpairedMsa"].startswith(">query")
    assert "unpairedMsaPath" not in chain


def test_af_server_splits_one_job_per_seed(tmp_path, msa):
    job = Job("A__B", (Component(A), Component(B)), seeds=(1, 2, 3))
    written = af_server_batch([job], tmp_path, msa)
    entries = json.loads(written[0].path.read_text())
    assert len(entries) == 3
    assert [e["modelSeeds"] for e in entries] == [["1"], ["2"], ["3"]]
    assert {e["name"] for e in entries} == {"A__B_s1", "A__B_s2", "A__B_s3"}


def test_af_server_batches_are_sized_to_the_daily_quota(tmp_path, msa):
    jobs = [Job(f"J{i}", (Component(A),), seeds=(1,)) for i in range(25)]
    written = af_server_batch(jobs, tmp_path, msa, quota=20)
    assert len({w.path for w in written}) == 2
    assert len(json.loads(sorted({w.path for w in written})[0].read_text())) == 20


def test_af_server_refuses_a_ligand_outside_its_allowlist(tmp_path):
    job = Job("x", (Component(A), Component(Ligand("drug", smiles="CCO"))))
    with pytest.raises(IRError, match="never custom SMILES"):
        get_engine("af_server")(job, tmp_path)


def test_af_server_accepts_an_allowlisted_ligand(tmp_path):
    assert "ATP" in LIGAND_ALLOWLIST
    job = Job("x", (Component(A), Component(Ligand("atp", ccd="ATP"))))
    entry = json.loads(get_engine("af_server")(job, tmp_path).path.read_text())[0]
    assert entry["sequences"][1]["ligand"]["ligand"] == "CCD_ATP"


# --- ColabFold --------------------------------------------------------------

def test_colabfold_joins_chains_with_a_colon(hetero, homo):
    assert query_sequence(hetero) == f"{A.sequence}:{B.sequence}"
    assert query_sequence(homo) == f"{A.sequence}:{A.sequence}"


def test_colabfold_a3m_header_gives_lengths_and_copies(homo, msa):
    text = build_a3m(homo, msa)
    assert text.splitlines()[0] == f"#{len(A.sequence)}\t2"


def test_colabfold_a3m_header_lists_each_unique_chain(hetero, msa):
    text = build_a3m(hetero, msa)
    assert text.splitlines()[0] == f"#{len(A.sequence)},{len(B.sequence)}\t1,1"


def test_colabfold_query_row_is_the_concatenated_chains(hetero, msa):
    lines = build_a3m(hetero, msa).splitlines()
    assert lines[1] == ">101"
    assert lines[2] == A.sequence + B.sequence


def test_colabfold_csv_carries_the_alignment_path(tmp_path, hetero, msa):
    written = [get_engine("colabfold")(hetero, tmp_path, msa)]
    path = write_queries_csv(written, tmp_path)
    header, row = path.read_text().splitlines()[:2]
    assert header == "id,sequence,a3mpath"
    assert row.endswith(".a3m")


# --- AlphaFold2-Multimer ----------------------------------------------------

def test_af2_writes_one_record_per_chain_starting_at_b(tmp_path, homo):
    written = get_engine("af2_multimer")(homo, tmp_path)
    text = written.path.read_text()
    assert text.count(">") == 2
    assert ":" not in text  # the colon form belongs to ColabFold
    assert written.chains == [["B", "C"]]


# --- SeedFold ---------------------------------------------------------------

def test_seedfold_single_form_is_a_bare_entity_list(tmp_path, homo):
    document = json.loads(get_engine("seedfold")(homo, tmp_path).path.read_text())
    assert document == [{"entity": "Protein", "copies": 2, "sequence": A.sequence}]


def test_seedfold_batch_form_wraps_entities_in_a_job(tmp_path, hetero, homo):
    written = seedfold_batch([hetero, homo], tmp_path)
    entries = json.loads(written[0].path.read_text())
    assert entries[0]["job_name"] == "A__B"
    assert entries[0]["model"] == "SeedFold_v1.0.0"
    assert len(entries[0]["entities"]) == 2


def test_seedfold_cannot_take_an_alignment(tmp_path, hetero, msa):
    """The one engine that always searches for itself; comparisons must say so."""
    assert get_engine("seedfold")(hetero, tmp_path, msa).msa_used is False


def test_seedfold_puts_smiles_in_the_same_field_as_sequences(tmp_path):
    job = Job("x", (Component(A), Component(Ligand("drug", smiles="CCO"))))
    document = json.loads(get_engine("seedfold")(job, tmp_path).path.read_text())
    assert document[1] == {"entity": "Ligand-Smiles", "copies": 1, "sequence": "CCO"}


# --- cross-engine -----------------------------------------------------------

def test_the_same_homodimer_takes_four_different_shapes(tmp_path, homo, msa):
    """Copies are a count in some engines and repeated records in others."""
    af3 = json.loads(get_engine("af3")(homo, tmp_path / "af3").path.read_text())
    boltz = yaml.safe_load(get_engine("boltz2")(homo, tmp_path / "b").path.read_text())
    protenix = json.loads(get_engine("protenix")(homo, tmp_path / "p").path.read_text())
    chai = get_engine("chai1")(homo, tmp_path / "c").path.read_text()

    assert af3["sequences"][0]["protein"]["id"] == ["A", "B"]
    assert boltz["sequences"][0]["protein"]["id"] == ["A", "B"]
    assert protenix[0]["sequences"][0]["proteinChain"]["count"] == 2
    assert chai.count(">protein") == 2


def test_one_engine_failing_does_not_cost_the_panel_the_others(tmp_path, monkeypatch):
    """A missing optional dependency or a refused ligand must not abort the run."""
    from foldrunner.cli import main
    from foldrunner.engines import base

    fasta = tmp_path / "lib.fasta"
    fasta.write_text(f">A\n{A.sequence}\n>B\n{B.sequence}\n")

    def explode(job, outdir, msa=None):
        raise base.EngineDependencyError("pandas is not installed")

    monkeypatch.setitem(base.ENGINES, "chai1", explode)
    out = tmp_path / "out"
    assert main(["write", str(fasta), "-o", str(out), "--engines", "boltz2,chai1,af3"]) == 0
    assert len(list((out / "boltz2").glob("*.yaml"))) == 3
    assert len(list((out / "af3").glob("*.json"))) == 3
    assert not (out / "chai1").exists()


def test_a_missing_optional_dependency_names_the_extra(tmp_path, monkeypatch):
    import builtins

    from foldrunner.engines import chai1

    real = builtins.__import__

    def blocked(name, *rest):
        if name == "pandas":
            raise ImportError("no pandas")
        return real(name, *rest)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(chai1.EngineDependencyError, match=r'foldrunner\[chai\]'):
        chai1.aligned_frame(None, None)


def test_colabfold_keeps_alignments_out_of_the_query_directory(tmp_path, hetero, msa):
    """colabfold_batch reads every .fasta and .a3m in a directory it is given, so a
    flat layout would make it enumerate each complex twice."""
    written = get_engine("colabfold")(hetero, tmp_path, msa)
    assert written.msa_used
    assert list(tmp_path.glob("*.a3m")) == []
    assert (tmp_path / "msas" / "A__B.a3m").is_file()
    assert sorted(p.suffix for p in tmp_path.iterdir() if p.is_file()) == [".fasta"]


def test_colabfold_csv_still_points_at_the_moved_alignment(tmp_path, hetero, msa):
    written = [get_engine("colabfold")(hetero, tmp_path, msa)]
    path = write_queries_csv(written, tmp_path)
    row = path.read_text().splitlines()[1]
    assert "msas/A__B.a3m" in row


def test_protenix_command_asks_for_the_data_interface_scores_need():
    """Without this flag Protenix writes no token_pair_pae, and every PAE-based
    interface metric is lost for the whole panel."""
    from foldrunner.engines.protenix import COMMAND

    assert "--need_atom_confidence True" in COMMAND


def test_chai_command_avoids_the_cli_that_discards_pae():
    from foldrunner.engines.chai1 import COMMAND

    assert "chai-lab fold" not in COMMAND
    assert "{script}" in COMMAND


# --- failures that only a real run exposes -----------------------------------

def test_af3_supplies_templates_alongside_alignments(tmp_path, hetero, msa):
    """Supplying alignments skips the data pipeline, which then cannot fill
    templates either; AlphaFold 3 refuses the job after loading the model."""
    document = json.loads(get_engine("af3")(hetero, tmp_path, msa).path.read_text())
    for entry in document["sequences"]:
        assert entry["protein"]["templates"] == []


def test_af3_without_alignments_leaves_templates_to_the_pipeline(tmp_path, hetero):
    document = json.loads(get_engine("af3")(hetero, tmp_path).path.read_text())
    assert "templates" not in document["sequences"][0]["protein"]


def test_chai_gives_every_copy_a_distinct_name(tmp_path, homo):
    """Chai-1 refuses a file where two records share a name."""
    text = get_engine("chai1")(homo, tmp_path).path.read_text()
    names = [line.split("name=")[1] for line in text.splitlines() if line.startswith(">")]
    assert len(names) == 2
    assert len(set(names)) == 2
    assert all(n.startswith("A") for n in names)


def test_chai_keeps_the_plain_name_for_a_single_copy(tmp_path, hetero):
    text = get_engine("chai1")(hetero, tmp_path).path.read_text()
    names = [line.split("name=")[1] for line in text.splitlines() if line.startswith(">")]
    assert names == ["A", "B"]


def test_colabfold_command_honours_a_configured_executable():
    """localcolabfold installs outside PATH; the config has to be able to say so."""
    from foldrunner.engines.colabfold import COMMAND, DEFAULT_EXECUTABLE

    assert "{executable}" in COMMAND
    assert DEFAULT_EXECUTABLE == "colabfold_batch"
