"""The programmatic interface the CLI is a thin layer over."""

import pytest

import foldrunner as fr
from foldrunner.msa.backends.base import MsaResult
from foldrunner.msa.cache import MsaCache, Recipe

FASTA = ">A\nAQVINTFDGV\n>B\nKKAVINGEQI\n>C\nMEEPQSDPSV\n"


@pytest.fixture
def lib(tmp_path):
    path = tmp_path / "lib.fasta"
    path.write_text(FASTA)
    return path


class FakeBackend:
    """Answers with a canned alignment, counting how often it is asked."""

    def __init__(self):
        self.calls = 0
        self.sequences = []

    @property
    def recipe(self):
        return Recipe(backend="fake", databases=("uniref30",), tool_version="1")

    def fetch(self, sequences, *, paired=False):
        self.calls += 1
        self.sequences.extend(sequences)
        return [
            MsaResult(
                sequence=s,
                unpaired=None if paired else f">query\n{s}\n>ERR1_2\n{s}\n",
                paired=f">query\n{s}\n>UniRef100_X_9606/\n{s}\n" if paired else None,
            )
            for s in sequences
        ]


def test_library_collapses_duplicate_sequences(tmp_path):
    path = tmp_path / "dup.fasta"
    path.write_text(">A\nAQVINTFDGV\n>A2\naqvintfdgv\n")
    assert len(fr.library(path)) == 1


def test_panel_enumerates_all_pairs_including_self(lib):
    jobs = fr.panel(lib, "all")
    assert [j.name for j in jobs] == ["A__A", "A__B", "A__C", "B__B", "B__C", "C__C"]


def test_cost_reports_the_saving_that_motivates_the_package(lib):
    jobs = fr.panel(lib, "all")
    summary = fr.cost(jobs, engines=["boltz2", "af3"])
    assert summary.complexes == 6
    assert summary.sequences == 3
    assert summary.searches_saved == 3
    assert summary.files == 12


def test_search_asks_for_each_sequence_once_not_each_pair(lib, tmp_path):
    jobs = fr.panel(lib, "all")
    backend = FakeBackend()
    searched = fr.search(jobs, tmp_path / "msa", backend, paired=False)
    assert searched == 3            # three sequences, not six complexes
    assert backend.calls == 1       # and one batched request
    assert sorted(set(backend.sequences)) == sorted(backend.sequences)


def test_search_skips_what_is_already_cached(lib, tmp_path):
    jobs = fr.panel(lib, "all")
    backend = FakeBackend()
    fr.search(jobs, tmp_path / "msa", backend, paired=False)
    assert fr.search(jobs, tmp_path / "msa", backend, paired=False) == 0


def test_write_produces_inputs_a_manifest_and_a_snapshot(lib, tmp_path):
    jobs = fr.panel(lib, "all", seeds=(2066,))
    written = fr.write(jobs, tmp_path / "panel", engines=["boltz2", "chai1"])
    assert len(written) == 12
    assert (tmp_path / "panel" / "manifest.tsv").is_file()
    assert (tmp_path / "panel" / "panel.json").is_file()
    assert len(list((tmp_path / "panel" / "boltz2").glob("*.yaml"))) == 6


def test_write_with_a_cache_wires_the_alignments_in(lib, tmp_path):
    jobs = fr.panel(lib, "all")
    backend = FakeBackend()
    fr.search(jobs, tmp_path / "msa", backend)
    written = fr.write(
        jobs,
        tmp_path / "panel",
        engines=["boltz2"],
        cache=tmp_path / "msa",
        recipe=backend.recipe,
    )
    assert all(item.msa_used for item in written)


def test_a_cache_without_its_recipe_is_refused(lib, tmp_path):
    with pytest.raises(ValueError, match="recipe"):
        fr.write(fr.panel(lib), tmp_path / "panel", cache=tmp_path / "msa")


def test_run_plan_builds_invocations_and_scripts(lib, tmp_path):
    jobs = fr.panel(lib, "all", seeds=(2066,))
    written = fr.write(jobs, tmp_path / "panel", engines=["boltz2", "af_server"])
    invocations = fr.run_plan(
        written, tmp_path / "results", scripts=tmp_path / "scripts", config=None
    )
    assert any(i.manual for i in invocations)          # af_server is a web form
    assert any(not i.manual for i in invocations)      # boltz2 is runnable
    assert (tmp_path / "scripts" / "run_all.sh").is_file()


def test_import_alignments_reports_what_it_found(tmp_path):
    source = tmp_path / "old" / "0"
    source.mkdir(parents=True)
    (source / "pairing.a3m").write_text(">query\nAQVINTFDGV\n>UniRef100_X_9606/\nAQVINTFDGV\n")
    recipe = Recipe(backend="imported", databases=("uniref100",), tool_version="x")
    report = fr.import_alignments(tmp_path / "old", tmp_path / "msa", recipe)
    assert report.imported == 1
    assert MsaCache(tmp_path / "msa").keys(recipe)


def test_the_package_exposes_a_version():
    assert fr.__version__


# --- the cache has to be able to identify itself ----------------------------

def test_search_records_the_recipe_beside_the_alignments(lib, tmp_path):
    """A recipe id hashes every field, so retyping a few flags cannot reproduce
    one; the cache has to say what it holds."""
    from foldrunner.msa.cache import MsaCache

    backend = FakeBackend()
    fr.search(fr.panel(lib, "all"), tmp_path / "msa", backend, paired=False)
    stored = MsaCache(tmp_path / "msa").stored_recipes()
    assert [r.recipe_id for r in stored] == [backend.recipe.recipe_id]


def test_storing_the_second_half_keeps_the_first(lib, tmp_path):
    """The two halves arrive in separate passes; the later one must not report
    the earlier as empty."""
    from foldrunner.msa.cache import MsaCache

    jobs = fr.panel(lib, "all")
    backend = FakeBackend()
    fr.search(jobs, tmp_path / "msa", backend, paired=True)
    cache = MsaCache(tmp_path / "msa")
    for key in cache.keys(backend.recipe):
        entry = cache.read_meta(backend.recipe, key)
        assert entry.depth_unpaired > 0, "unpaired depth was clobbered"
        assert entry.depth_paired > 0
