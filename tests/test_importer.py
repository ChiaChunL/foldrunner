import os

import pytest

from foldrunner.ir import sequence_key
from foldrunner.msa.cache import MsaCache, Recipe
from foldrunner.msa.importer import discover, import_entry, import_tree

SEQ = "MKVLAAGIVG"
PAIRED = f">query\n{SEQ}\n>UniRef100_A_9606/\t1\t2\t3\n{SEQ}\n>UniRef100_B_10090/\t1\t2\t3\n{SEQ}\n"
UNPAIRED = f">query\n{SEQ}\n>ERR1_2\t1\t2\t3\n{SEQ}\n"
UNIREF = f">query_0\n{SEQ}\n>UniRef100_C\t1\t2\t3\n{SEQ}\n"
BFD = f">query_0\n{SEQ}\n>ERR9_9\t1\t2\t3\n{SEQ}\n"


@pytest.fixture
def recipe():
    return Recipe(backend="imported", databases=("uniref100",), tool_version="x")


def make_protenix(root, index, seq=SEQ):
    directory = root / str(index)
    directory.mkdir(parents=True)
    (directory / "pairing.a3m").write_text(PAIRED.replace(SEQ, seq))
    (directory / "non_pairing.a3m").write_text(UNPAIRED.replace(SEQ, seq))
    return directory


def test_key_is_recovered_from_the_alignment_itself(tmp_path, recipe):
    """The first record of an a3m is the query, so no index-to-FASTA map is needed."""
    make_protenix(tmp_path, 7)
    cache = MsaCache(tmp_path / "cache")
    key, written = import_entry(next(discover(tmp_path)), cache, recipe)
    assert key == sequence_key(SEQ)
    assert written


def test_position_named_directories_collapse_by_content(tmp_path, recipe):
    for index in (0, 1, 2):
        make_protenix(tmp_path / "chunk", index)
    report = import_tree(tmp_path / "chunk", MsaCache(tmp_path / "cache"), recipe)
    assert report.scanned == 3
    assert len(report.keys) == 1
    assert len(report.collisions) == 2


def test_distinct_sequences_stay_distinct(tmp_path, recipe):
    make_protenix(tmp_path / "chunk", 0, "MKVLAAGIVG")
    make_protenix(tmp_path / "chunk", 1, "PQRSTVWYAC")
    report = import_tree(tmp_path / "chunk", MsaCache(tmp_path / "cache"), recipe)
    assert report.imported == 2
    assert len(report.keys) == 2


def test_import_is_idempotent(tmp_path, recipe):
    make_protenix(tmp_path / "chunk", 0)
    cache = MsaCache(tmp_path / "cache")
    first = import_tree(tmp_path / "chunk", cache, recipe)
    second = import_tree(tmp_path / "chunk", cache, recipe)
    assert first.imported == 1 and second.imported == 0
    assert second.skipped == 1


def test_depth_and_species_are_recorded(tmp_path, recipe):
    make_protenix(tmp_path, 0)
    cache = MsaCache(tmp_path / "cache")
    import_tree(tmp_path, cache, recipe)
    entry = cache.read_meta(recipe, sequence_key(SEQ))
    assert entry.depth_paired == 3
    assert entry.depth_unpaired == 2
    assert entry.n_species == 2


def test_colabfold_layout_merges_its_per_database_files(tmp_path, recipe):
    directory = tmp_path / "0"
    directory.mkdir()
    (directory / "uniref.a3m").write_text(UNIREF)
    (directory / "bfd.mgnify30.metaeuk30.smag30.a3m").write_text(BFD)
    cache = MsaCache(tmp_path / "cache")
    import_tree(tmp_path, cache, recipe)
    # One query row kept, plus one hit from each database.
    assert cache.read_meta(recipe, sequence_key(SEQ)).depth_unpaired == 3


def test_mismatched_pair_is_reported_not_imported(tmp_path, recipe):
    directory = tmp_path / "0"
    directory.mkdir()
    (directory / "pairing.a3m").write_text(PAIRED)
    (directory / "non_pairing.a3m").write_text(UNPAIRED.replace(SEQ, "PQRSTVWYAC"))
    report = import_tree(tmp_path, MsaCache(tmp_path / "cache"), recipe)
    assert report.imported == 0
    assert "different query sequences" in report.failed[0][1]


def test_a_broken_entry_does_not_stop_the_rest(tmp_path, recipe):
    make_protenix(tmp_path / "chunk", 0, "MKVLAAGIVG")
    bad = tmp_path / "chunk" / "1"
    bad.mkdir()
    (bad / "pairing.a3m").write_text("")
    report = import_tree(tmp_path / "chunk", MsaCache(tmp_path / "cache"), recipe)
    assert report.imported == 1
    assert len(report.failed) == 1


def test_dry_run_writes_nothing(tmp_path, recipe):
    make_protenix(tmp_path, 0)
    cache = MsaCache(tmp_path / "cache")
    report = import_tree(tmp_path, cache, recipe, dry_run=True)
    assert report.imported == 1
    assert cache.keys(recipe) == []


def test_link_mode_shares_bytes_with_the_source(tmp_path, recipe):
    directory = make_protenix(tmp_path, 0)
    cache = MsaCache(tmp_path / "cache")
    import_tree(tmp_path, cache, recipe, link=True)
    target = cache.entry_dir(recipe, sequence_key(SEQ)) / "paired.a3m"
    assert os.stat(directory / "pairing.a3m").st_ino == os.stat(target).st_ino


def test_limit_stops_early(tmp_path, recipe):
    for index in range(5):
        make_protenix(tmp_path / "chunk", index, "MKVLAAGIV" + "ACDEF"[index])
    report = import_tree(tmp_path / "chunk", MsaCache(tmp_path / "cache"), recipe, limit=2)
    assert report.scanned == 2


def test_query_is_read_without_loading_the_whole_file(tmp_path):
    from foldrunner.msa.a3m import read_query

    path = tmp_path / "big.a3m"
    path.write_text(f">query\n{SEQ}\n" + "".join(f">hit{i}\n{SEQ}\n" for i in range(10_000)))
    assert read_query(path) == SEQ


def test_scan_counts_records_and_taxa_by_streaming(tmp_path):
    from foldrunner.msa.a3m import scan

    path = tmp_path / "p.a3m"
    path.write_text(PAIRED)
    assert scan(path) == (3, 2)


def test_an_empty_query_record_is_rejected(tmp_path):
    from foldrunner.msa.a3m import read_query

    path = tmp_path / "bad.a3m"
    path.write_text(">query\n")
    with pytest.raises(ValueError, match="empty query record"):
        read_query(path)


def make_colabfold(root, index, seq=SEQ):
    directory = root / str(index)
    directory.mkdir(parents=True)
    (directory / "uniref.a3m").write_text(UNIREF.replace(SEQ, seq))
    (directory / "bfd.mgnify30.metaeuk30.smag30.a3m").write_text(BFD.replace(SEQ, seq))
    return directory


def test_dry_run_counts_distinct_sequences_not_directories(tmp_path, recipe):
    """Nothing is written, so novelty has to come from the keys already seen."""
    for index in (0, 1, 2):
        make_protenix(tmp_path / "chunk", index)
    report = import_tree(tmp_path / "chunk", MsaCache(tmp_path / "c"), recipe, dry_run=True)
    assert report.scanned == 3
    assert report.imported == 1
    assert len(report.keys) == 1


def test_the_same_sequence_from_two_tools_is_flagged_not_silently_merged(tmp_path, recipe):
    make_protenix(tmp_path / "a", 0)
    make_colabfold(tmp_path / "b", 0)
    report = import_tree(tmp_path, MsaCache(tmp_path / "c"), recipe, dry_run=True)
    assert len(report.conflicts) == 1
    assert "different source layout" in report.summary()
    assert "must not share a recipe" in report.summary()


def test_repeats_from_one_tool_are_plain_duplicates(tmp_path, recipe):
    for index in (0, 1):
        make_protenix(tmp_path / "chunk", index)
    report = import_tree(tmp_path / "chunk", MsaCache(tmp_path / "c"), recipe, dry_run=True)
    assert report.collisions and not report.conflicts


def test_layout_restriction_keeps_one_source(tmp_path, recipe):
    make_protenix(tmp_path / "a", 0)
    make_colabfold(tmp_path / "b", 0)
    report = import_tree(
        tmp_path, MsaCache(tmp_path / "c"), recipe, layout="protenix", dry_run=True
    )
    assert report.scanned == 1
    assert not report.conflicts
