import pytest
import yaml

from foldrunner.engines import get_engine
from foldrunner.engines.base import MsaView, assign_chains, chain_labels
from foldrunner.engines.boltz2 import build as boltz_build
from foldrunner.enumerate import enumerate_jobs, unique_sequences
from foldrunner.ir import Component, Job, Ligand, Protein, sequence_key
from foldrunner.msa.cache import MsaCache, Recipe

BARNASE = Protein("BARNASE", "AQVINTFDGVADYLQ")
BARSTAR = Protein("BARSTAR", "KKAVINGEQIRSISD")

PAIRED_TEXT = ">query\n{seq}\n>UniRef100_A_9606/\t1\t2\t3\n{seq}\n"
UNPAIRED_TEXT = ">query\n{seq}\n>ERR1_2\t1\t2\t3\n{seq}\n"


@pytest.fixture
def msa(tmp_path):
    cache = MsaCache(tmp_path / "msa")
    recipe = Recipe(backend="test", databases=("uniref30",), tool_version="18")
    for protein in (BARNASE, BARSTAR):
        cache.store(
            recipe,
            protein.key,
            protein.sequence,
            unpaired=UNPAIRED_TEXT.format(seq=protein.sequence),
            paired=PAIRED_TEXT.format(seq=protein.sequence),
        )
    return MsaView(cache=cache, recipe=recipe)


def test_chain_labels_roll_over_past_z():
    labels = chain_labels()
    first = [next(labels) for _ in range(26)]
    assert first[0] == "A" and first[-1] == "Z"
    assert next(labels) == "AA"


def test_af2_style_writers_can_start_at_b():
    job = Job("x", (Component(BARNASE), Component(BARSTAR)))
    assert assign_chains(job, start="B") == [["B"], ["C"]]


def test_boltz_expresses_a_homodimer_as_one_entity_with_two_ids():
    job = Job("A__A", (Component(BARNASE, 2),))
    document, _ = boltz_build(job)
    assert document["sequences"][0]["protein"]["id"] == ["A", "B"]
    assert len(document["sequences"]) == 1


def test_chai_expands_a_homodimer_into_repeated_records(tmp_path):
    job = Job("A__A", (Component(BARNASE, 2),))
    written = get_engine("chai1")(job, tmp_path)
    assert written.path.read_text().count(">protein|name=BARNASE") == 2


def test_chai_headers_carry_the_entity_type(tmp_path):
    job = Job(
        "x",
        (Component(BARNASE), Component(Ligand("lig", smiles="CCO"))),
        affinity_binder="lig",
    )
    text = get_engine("chai1")(job, tmp_path).path.read_text()
    assert ">protein|name=BARNASE" in text
    assert ">ligand|name=lig" in text


def test_affinity_jobs_are_kept_out_of_the_plain_boltz_directory(tmp_path):
    plain = Job("plain", (Component(BARNASE), Component(BARSTAR)))
    binding = Job(
        "binding",
        (Component(BARNASE), Component(Ligand("lig", smiles="CCO"))),
        affinity_binder="lig",
    )
    write = get_engine("boltz2")
    plain_out = write(plain, tmp_path)
    binding_out = write(binding, tmp_path)
    assert plain_out.path.parent != binding_out.path.parent
    document = yaml.safe_load(binding_out.path.read_text())
    assert document["properties"] == [{"affinity": {"binder": "B"}}]


def test_boltz_points_at_the_cached_alignment(tmp_path, msa):
    job = Job("A__B", (Component(BARNASE), Component(BARSTAR)))
    written = get_engine("boltz2")(job, tmp_path, msa)
    document = yaml.safe_load(written.path.read_text())
    assert written.msa_used
    assert document["sequences"][0]["protein"]["msa"].endswith("unpaired.a3m")


def test_chai_writes_one_alignment_directory_per_unique_sequence(tmp_path, msa):
    job = Job("A__B", (Component(BARNASE), Component(BARSTAR)))
    written = get_engine("chai1")(job, tmp_path, msa)
    dirs = sorted(d.name for d in (tmp_path / "msas_a3m").iterdir() if d.is_dir())
    assert dirs == sorted(p.key for p in (BARNASE, BARSTAR))
    assert written.extra["msa_directory"].endswith("msas")
    assert written.extra["msa_a3m_directory"].endswith("msas_a3m")


def test_chai_alignments_are_left_as_a3m_for_the_engine_to_convert(tmp_path, msa):
    """Parquet written by a newer pyarrow fails to load on an older one, and the
    reader's environment is the engine's, not this package's."""
    job = Job("A__B", (Component(BARNASE), Component(BARSTAR)))
    get_engine("chai1")(job, tmp_path, msa)
    assert list(tmp_path.rglob("*.aligned.pqt")) == []
    names = {f.name for f in (tmp_path / "msas_a3m" / BARNASE.key).iterdir()}
    assert names == {"uniprot.a3m", "uniref90.a3m"}


def test_chai_frame_keeps_taxa_only_on_pairable_rows(tmp_path, msa):
    """The frame builder is still exercised directly; only the runner path
    defers the write to the engine's own environment."""
    from foldrunner.engines.chai1 import aligned_frame
    from foldrunner.msa.a3m import parse_a3m

    paired = parse_a3m(PAIRED_TEXT.format(seq=BARNASE.sequence))
    unpaired = parse_a3m(UNPAIRED_TEXT.format(seq=BARNASE.sequence))
    frame = aligned_frame(paired, unpaired)
    assert list(frame.columns) == ["sequence", "source_database", "pairing_key", "comment"]
    assert set(frame.loc[frame.source_database == "uniprot", "pairing_key"]) == {"9606"}
    assert set(frame.loc[frame.source_database == "uniref90", "pairing_key"]) == {""}


def test_one_alignment_set_feeds_both_engines(tmp_path, msa):
    """The panel's whole reason for existing: search once, write many."""
    jobs = enumerate_jobs([BARNASE, BARSTAR], "all")
    assert len(jobs) == 3
    assert len(unique_sequences(jobs)) == 2

    for job in jobs:
        boltz = get_engine("boltz2")(job, tmp_path / "boltz2", msa)
        chai = get_engine("chai1")(job, tmp_path / "chai1", msa)
        assert boltz.msa_used and chai.msa_used

    assert len(list((tmp_path / "boltz2").glob("*.yaml"))) == 3
    assert len(list((tmp_path / "chai1").glob("*.fasta"))) == 3
    # Two sequences searched, three complexes written.
    a3m = tmp_path / "chai1" / "msas_a3m"
    assert len([d for d in a3m.iterdir() if d.is_dir()]) == 2


def test_manifest_rows_record_what_each_engine_named_the_chains(tmp_path, msa):
    job = Job("A__B", (Component(BARNASE), Component(BARSTAR)))
    rows = get_engine("boltz2")(job, tmp_path, msa).chain_rows()
    assert rows == [
        ("A__B", "boltz2", 0, "BARNASE", "A"),
        ("A__B", "boltz2", 1, "BARSTAR", "B"),
    ]


def test_cache_key_matches_chai_naming_convention():
    assert sequence_key(BARNASE.sequence) == BARNASE.key
