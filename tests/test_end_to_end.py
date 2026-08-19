"""The whole chain, on the shipped example library, without a GPU.

Every failure this project has actually hit came from a stage the unit tests do
not reach: a writer raising and taking the panel down with it, a generated
script that does not parse, a panel whose recorded paths only work on the
machine that wrote it. None of those need a model to catch — they need the
pipeline run end to end, which is what this does.
"""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import pytest

from foldrunner.cli import main
from foldrunner.engines import ENGINES
from foldrunner.manifest import PANEL

LIBRARY = Path(__file__).resolve().parents[1] / "examples" / "data" / "library.fasta"

# Three sequences, so six complexes: three homodimers on the diagonal and three
# hetero pairs.
COMPLEXES = 6


@pytest.fixture(scope="module")
def panel(tmp_path_factory):
    out = tmp_path_factory.mktemp("panel")
    assert main(["write", str(LIBRARY), "-o", str(out), "--seeds", "2066,318"]) == 0
    return out


def test_the_example_library_ships_with_the_package():
    assert LIBRARY.is_file(), "examples/data/library.fasta is referenced by the README"


def test_plan_reports_the_alignment_saving(capsys):
    assert main(["plan", str(LIBRARY)]) == 0
    out = capsys.readouterr().out
    assert f"complexes      {COMPLEXES}" in out
    assert "MSA searches   3" in out


def test_every_registered_engine_produces_input(panel):
    for name in ENGINES:
        directory = panel / name
        assert directory.is_dir(), f"{name} wrote nothing"
        assert any(directory.iterdir()), f"{name} produced an empty directory"


def test_the_manifest_covers_every_engine_and_complex(panel):
    with (panel / "manifest.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["engine"] for r in rows} == set(ENGINES)
    for name in ENGINES:
        jobs = {r["job"] for r in rows if r["engine"] == name}
        assert len(jobs) == COMPLEXES, f"{name} covers {len(jobs)} complexes"


def test_the_panel_snapshot_records_no_absolute_paths(panel):
    snapshot = json.loads((panel / PANEL).read_text())
    offenders = [
        item["path"] for item in snapshot["written"] if Path(item["path"]).is_absolute()
    ]
    assert offenders == [], "a panel with absolute paths cannot be run elsewhere"


def test_generated_scripts_are_valid_shell(panel, tmp_path):
    """A script that does not parse fails only once the queue reaches it."""
    assert main(["run", str(panel), "-o", str(tmp_path / "results"), "--runner", "script"]) == 0
    scripts = sorted((panel / "scripts").glob("*.sh"))
    assert scripts, "no scripts were written"
    for script in scripts:
        check = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True)
        assert check.returncode == 0, f"{script.name} is not valid shell: {check.stderr}"


def test_web_only_engines_are_listed_for_manual_upload(panel, tmp_path):
    main(["run", str(panel), "-o", str(tmp_path / "results"), "--runner", "script"])
    notes = (panel / "scripts" / "MANUAL.md").read_text()
    assert "af_server" in notes
    assert "seedfold" in notes


def test_a_homodimer_is_one_entity_with_two_copies_everywhere(panel):
    """Six of the eight engines say this with a count; two repeat the record."""
    import yaml

    name = "BARNASE__BARNASE"
    boltz = yaml.safe_load((panel / "boltz2" / f"{name}.yaml").read_text())
    assert boltz["sequences"][0]["protein"]["id"] == ["A", "B"]

    af3 = json.loads((panel / "af3" / f"{name}.json").read_text())
    assert af3["sequences"][0]["protein"]["id"] == ["A", "B"]

    protenix = json.loads((panel / "protenix" / "panel.json").read_text())
    entry = next(e for e in protenix if e["name"] == name)
    assert entry["sequences"][0]["proteinChain"]["count"] == 2

    chai = (panel / "chai1" / f"{name}.fasta").read_text()
    assert chai.count(">protein") == 2


def test_a_second_write_is_idempotent(panel, tmp_path):
    """Re-running write must not duplicate manifest rows."""
    before = (panel / "manifest.tsv").read_text().count("\n")
    assert main(["write", str(LIBRARY), "-o", str(panel), "--seeds", "2066,318"]) == 0
    assert (panel / "manifest.tsv").read_text().count("\n") == before


def test_write_uses_the_recipe_the_cache_recorded(tmp_path):
    """The documented flow is search then write; it must not need the caller to
    retype parameters that only the backend knows."""
    from foldrunner.msa.backends.base import MsaResult
    from foldrunner.msa.cache import MsaCache, Recipe

    class Backend:
        @property
        def recipe(self):
            return Recipe(
                backend="protenix-api",
                databases=("uniref30", "colabfold_envdb"),
                tool_version="mmseqs2-service",
                params=(("pairing_strategy", "greedy"),),
            )

        def fetch(self, sequences, *, paired=False):
            return [
                MsaResult(
                    sequence=s,
                    unpaired=None if paired else f">query\n{s}\n>ERR1_2\n{s}\n",
                    paired=f">query\n{s}\n>UniRef100_X_9606/\n{s}\n" if paired else None,
                )
                for s in sequences
            ]

    import foldrunner as fr

    jobs = fr.panel(LIBRARY, "all")
    backend = Backend()
    cache = tmp_path / "msa"
    fr.search(jobs, cache, backend)
    MsaCache(cache).write_recipe(backend.recipe)

    # No --backend/--databases/--tool-version: the cache knows.
    out = tmp_path / "panel"
    assert main(["write", str(LIBRARY), "-o", str(out), "--cache", str(cache),
                 "--engines", "boltz2"]) == 0
    text = "\n".join(p.read_text() for p in (out / "boltz2").glob("*.yaml"))
    assert "msa:" in text, "the cached alignments were not attached"


def test_a_cache_that_matches_nothing_is_an_error_not_a_silent_downgrade(tmp_path):
    """Writing inputs that quietly carry no alignment produces a panel that
    either wastes a search per complex or refuses to start."""
    empty = tmp_path / "msa"
    (empty / "protenix-api-deadbeef0000").mkdir(parents=True)
    (empty / "protenix-api-deadbeef0000" / "recipe.json").write_text(
        '{"backend": "protenix-api", "databases": [], "tool_version": "", "params": []}'
    )
    code = main(["write", str(LIBRARY), "-o", str(tmp_path / "panel"),
                 "--cache", str(empty), "--engines", "boltz2"])
    assert code == 2
