import pytest

from foldrunner.ir import sequence_key
from foldrunner.msa.cache import PAIRED, UNPAIRED, CacheError, MsaCache, Recipe, looks_like_gzip

PAIRED_TEXT = ">query\nMKVL\n>UniRef100_A_9606/\nMKVI\n>UniRef100_B_10090/\nMKVV\n"
UNPAIRED_TEXT = ">query\nMKVL\n>ERR1_2\nMKCC\n"


@pytest.fixture
def cache(tmp_path):
    return MsaCache(tmp_path)


@pytest.fixture
def recipe():
    return Recipe(backend="colabfold-api", databases=("uniref30",), tool_version="18")


def test_recipe_id_changes_with_the_databases(recipe):
    other = Recipe(backend="colabfold-api", databases=("uniref30", "envdb"), tool_version="18")
    assert recipe.recipe_id != other.recipe_id


def test_recipe_id_is_stable_across_calls(recipe):
    assert recipe.recipe_id == Recipe(
        backend="colabfold-api", databases=("uniref30",), tool_version="18"
    ).recipe_id


def test_store_records_depth_and_species_count(cache, recipe):
    key = sequence_key("MKVL")
    entry = cache.store(recipe, key, "MKVL", unpaired=UNPAIRED_TEXT, paired=PAIRED_TEXT)
    assert entry.depth_unpaired == 2
    assert entry.depth_paired == 3
    assert entry.n_species == 2


def test_entries_are_addressed_by_content_not_by_position(cache, recipe):
    key = sequence_key("mkvl")
    cache.store(recipe, key, "MKVL", unpaired=UNPAIRED_TEXT)
    assert cache.has(recipe, sequence_key("MKVL"), UNPAIRED)


def test_missing_reports_what_still_needs_work(cache, recipe):
    key = sequence_key("MKVL")
    cache.store(recipe, key, "MKVL", unpaired=UNPAIRED_TEXT)
    assert cache.missing(recipe, [key, "beef"], UNPAIRED) == ["beef"]


def test_different_recipes_do_not_share_entries(cache, recipe):
    key = sequence_key("MKVL")
    cache.store(recipe, key, "MKVL", unpaired=UNPAIRED_TEXT)
    other = Recipe(backend="local-mmseqs", databases=("uniref30",), tool_version="18")
    assert not cache.has(other, key, UNPAIRED)


def test_reading_an_absent_entry_names_the_cache(cache, recipe):
    with pytest.raises(CacheError, match="no unpaired.a3m cached"):
        cache.read(recipe, "beef", UNPAIRED)


def test_no_temporary_files_survive_a_write(cache, recipe):
    key = sequence_key("MKVL")
    cache.store(recipe, key, "MKVL", unpaired=UNPAIRED_TEXT, paired=PAIRED_TEXT)
    leftovers = list(cache.entry_dir(recipe, key).glob(".*tmp"))
    assert leftovers == []


def test_throttled_html_response_is_not_mistaken_for_an_archive():
    assert looks_like_gzip(b"\x1f\x8b\x08\x00rest")
    assert not looks_like_gzip(b"<html><body>Too many requests")


def test_round_trip_through_the_cache_preserves_taxa(cache, recipe):
    key = sequence_key("MKVL")
    cache.store(recipe, key, "MKVL", paired=PAIRED_TEXT)
    assert cache.read(recipe, key, PAIRED).species == [None, "9606", "10090"]
