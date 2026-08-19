import pytest

from foldrunner.ir import Component, IRError, Job, Ligand, Protein, sequence_key


def test_sequence_key_is_case_insensitive_and_matches_chai_naming():
    assert Protein("a", "mkvl").key == Protein("b", "MKVL").key
    assert len(sequence_key("MKVL")) == 64


def test_whitespace_is_stripped_from_sequences():
    assert Protein("p", "MK VL\nAA").sequence == "MKVLAA"


def test_non_amino_acid_characters_are_rejected():
    with pytest.raises(IRError, match="non-amino-acid"):
        Protein("p", "MKVL123")


def test_names_are_reduced_to_safe_characters():
    assert Protein("sp|P04637|P53 HUMAN", "MKVL").name == "sp_P04637_P53_HUMAN"


def test_self_pair_is_one_entity_with_two_copies():
    p = Protein("A", "MKVL")
    job = Job("A__A", (Component(p, 2),))
    assert job.is_homomer
    assert job.n_chains == 2
    assert len(job.unique_protein_keys) == 1
    assert job.n_residues == 8


def test_hetero_pair_has_two_components():
    job = Job("A__B", (Component(Protein("A", "MKVL")), Component(Protein("B", "PQRS"))))
    assert not job.is_homomer
    assert job.n_chains == 2
    assert len(job.unique_protein_keys) == 2


def test_ligand_needs_smiles_or_ccd():
    with pytest.raises(IRError, match="smiles or ccd"):
        Ligand("x")


def test_ccd_prefix_is_normalised():
    assert Ligand("atp", ccd="CCD_atp").ccd == "ATP"


def test_affinity_binder_must_name_a_ligand_in_the_job():
    p = Protein("A", "MKVL")
    with pytest.raises(IRError, match="affinity binder"):
        Job("A", (Component(p),), affinity_binder="nope")


def test_ligands_do_not_count_toward_residues():
    job = Job(
        "A",
        (Component(Protein("A", "MKVL")), Component(Ligand("lig", smiles="CCO"))),
        affinity_binder="lig",
    )
    assert job.n_residues == 4
    assert job.n_chains == 2
