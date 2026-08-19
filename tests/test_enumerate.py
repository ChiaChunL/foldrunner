import pytest

from foldrunner.enumerate import (
    dedupe,
    enumerate_jobs,
    expected_job_count,
    iter_chunks,
    parse_fasta,
    unique_sequences,
)
from foldrunner.ir import IRError, Protein

LIB = [Protein("A", "MKVL"), Protein("B", "PQRS"), Protein("C", "WYFG")]


def test_all_mode_includes_self_pairs():
    jobs = enumerate_jobs(LIB, "all")
    assert [j.name for j in jobs] == ["A__A", "A__B", "A__C", "B__B", "B__C", "C__C"]
    assert len(jobs) == expected_job_count(3, "all")


def test_hetero_mode_drops_self_pairs():
    assert [j.name for j in enumerate_jobs(LIB, "hetero")] == ["A__B", "A__C", "B__C"]


def test_monomer_mode_gives_one_job_per_sequence():
    jobs = enumerate_jobs(LIB, "monomer")
    assert all(j.is_monomer for j in jobs)


def test_self_pair_becomes_a_homodimer():
    job = next(j for j in enumerate_jobs(LIB, "all") if j.name == "A__A")
    assert len(job.components) == 1
    assert job.components[0].count == 2


def test_msa_cost_grows_with_the_library_not_the_pairs():
    jobs = enumerate_jobs(LIB, "all")
    assert len(jobs) == 6
    assert len(unique_sequences(jobs)) == 3


def test_pair_names_are_order_independent():
    forward = enumerate_jobs([LIB[0], LIB[1]], "hetero")[0].name
    backward = enumerate_jobs([LIB[1], LIB[0]], "hetero")[0].name
    assert forward == backward == "A__B"


def test_bipartite_crosses_two_libraries_without_duplicates():
    jobs = enumerate_jobs([LIB[0]], "bipartite", preys=[LIB[1], LIB[0]])
    assert sorted(j.name for j in jobs) == ["A__A", "A__B"]


def test_bipartite_requires_preys():
    with pytest.raises(IRError, match="needs a second library"):
        enumerate_jobs(LIB, "bipartite")


def test_dedupe_collapses_identical_sequences():
    kept, duplicates = dedupe([Protein("A", "MKVL"), Protein("A2", "mkvl")])
    assert len(kept) == 1
    assert duplicates[0].dropped == "A2"


def test_parse_fasta_takes_first_header_field():
    proteins = parse_fasta(">sp|P1|NAME desc\nMKVL\n>B\nPQ\nRS\n")
    assert [p.name for p in proteins] == ["sp", "B"]
    assert proteins[1].sequence == "PQRS"


def test_chunking_splits_for_quota_limited_engines():
    jobs = enumerate_jobs(LIB, "all")
    assert [len(c) for c in iter_chunks(jobs, 4)] == [4, 2]
