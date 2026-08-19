"""Turning a written panel into commands."""

import pytest

from foldrunner.cli import main
from foldrunner.config import Config, EngineConfig, load
from foldrunner.manifest import PANEL, read_panel
from foldrunner.run import engine_scope, plan, write_scripts

CONFIG = """
[defaults]
conda_sh = "/opt/conda/etc/profile.d/conda.sh"

[engines.boltz2]
conda_env = "boltz2"

[engines.chai1]
conda_env = "chai1"
env = { HF_ENDPOINT = "https://mirror.example" }

[engines.protenix]
conda_env = "protenix"
env = { PROTENIX_ROOT_DIR = "/opt/protenix_data" }

[engines.af3]
container = "docker"
image = "alphafold3"
gpus = "all"

[engines.af2_multimer]
executable = "/opt/af2/run.sh"
"""


@pytest.fixture
def panel(tmp_path):
    fasta = tmp_path / "lib.fasta"
    fasta.write_text(">A\nAQVINTFDGV\n>B\nKKAVINGEQI\n")
    out = tmp_path / "panel"
    main(["write", str(fasta), "-o", str(out), "--seeds", "2066,318"])
    return out


@pytest.fixture
def conf(tmp_path):
    path = tmp_path / "foldrunner.toml"
    path.write_text(CONFIG)
    return load(path)


def replay(panel, engines=()):
    from foldrunner.cli import _replay

    return _replay(panel, list(engines))


# --- what the snapshot has to carry -----------------------------------------

def test_panel_snapshot_records_the_seeds(panel):
    assert read_panel(panel / PANEL)["seeds"] == [2066, 318]


def test_seeds_survive_into_the_command(panel, conf):
    """Reconstructing from the manifest alone would lose them and quietly run seed 0."""
    invocations = plan(replay(panel), conf, panel / "results", ["protenix"])
    assert "-s 2066,318" in invocations[0].script


def test_engine_extras_survive_into_the_command(panel, conf):
    """Chai-1 is driven by a generated script; losing its path would make the
    command try to execute the FASTA as Python."""
    invocations = plan(replay(panel), conf, panel / "results", ["chai1"])
    assert "run_chai.py" in invocations[0].script


# --- how much one invocation covers -----------------------------------------

def test_directory_engines_run_once_for_the_whole_panel(panel, conf):
    assert engine_scope("boltz2") == "directory"
    invocations = plan(replay(panel), conf, panel / "results", ["boltz2"])
    # Two seeds, one invocation each, not one per complex per seed.
    assert len(invocations) == 2
    assert {i.seed for i in invocations} == {2066, 318}


def test_per_complex_engines_run_once_each(panel, conf):
    invocations = plan(replay(panel), conf, panel / "results", ["af3"])
    assert len(invocations) == 3
    assert {i.label for i in invocations} == {"A__A", "A__B", "B__B"}


def test_batched_engines_run_once(panel, conf):
    invocations = plan(replay(panel), conf, panel / "results", ["protenix"])
    assert len(invocations) == 1
    assert invocations[0].label == "3 complexes"


def test_single_seed_engines_get_one_output_directory_per_seed(panel, conf):
    invocations = plan(replay(panel), conf, panel / "results", ["boltz2"])
    outs = {str(i.out) for i in invocations}
    assert len(outs) == 2
    assert any(o.endswith("_seed318") for o in outs)


# --- how the environment is entered -----------------------------------------

def test_conda_is_activated_after_being_sourced(panel, conf):
    script = plan(replay(panel), conf, panel / "results", ["boltz2"])[0].script
    lines = script.splitlines()
    assert lines[0].startswith("source /opt/conda")
    assert "conda activate boltz2" in lines[1]


def test_environment_variables_are_exported(panel, conf):
    script = plan(replay(panel), conf, panel / "results", ["chai1"])[0].script
    assert "export HF_ENDPOINT=https://mirror.example" in script


def test_container_engines_bind_paths_at_their_own_location(panel, conf):
    """Alignment paths inside an AlphaFold 3 job are absolute, so a container
    that mounts the input elsewhere cannot resolve them."""
    script = plan(replay(panel), conf, panel / "results", ["af3"])[0].script
    assert "docker run --rm --gpus all" in script
    mount = [p for p in script.split() if ":" in p and p.count("/") > 1][0]
    left, right = mount.rsplit(":", 1)
    assert left == right


def test_container_engines_pass_env_with_e_not_export(tmp_path, panel):
    conf = Config(engines={"af3": EngineConfig(name="af3", container="docker", image="af3",
                                               env={"K": "V"})})
    script = plan(replay(panel), conf, tmp_path / "r", ["af3"])[0].script
    assert "-e K=V" in script
    assert "export K=" not in script


def test_a_container_without_an_image_is_refused(panel, tmp_path):
    conf = Config(engines={"af3": EngineConfig(name="af3", container="docker")})
    with pytest.raises(ValueError, match="no image"):
        plan(replay(panel), conf, tmp_path / "r", ["af3"])


# --- web-form engines -------------------------------------------------------

def test_web_engines_are_reported_not_invented(panel, conf):
    invocations = plan(replay(panel), conf, panel / "results", ["af_server", "seedfold"])
    assert len(invocations) == 2  # one upload each, not one per complex
    assert all(i.manual for i in invocations)
    assert all("web form" in i.reason for i in invocations)
    assert all(not i.script for i in invocations)


# --- script output ----------------------------------------------------------

def test_scripts_are_written_and_executable(panel, conf, tmp_path):
    invocations = plan(replay(panel), conf, panel / "results")
    paths = write_scripts(invocations, tmp_path / "scripts")
    scripts = [p for p in paths if p.suffix == ".sh"]
    assert scripts
    assert all(p.stat().st_mode & 0o111 for p in scripts)
    assert (tmp_path / "scripts" / "run_all.sh").is_file()


def test_manual_steps_get_their_own_note(panel, conf, tmp_path):
    invocations = plan(replay(panel), conf, panel / "results")
    write_scripts(invocations, tmp_path / "scripts")
    notes = (tmp_path / "scripts" / "MANUAL.md").read_text()
    assert "af_server" in notes and "seedfold" in notes


def test_run_command_writes_scripts(panel, tmp_path):
    path = tmp_path / "foldrunner.toml"
    path.write_text(CONFIG)
    assert main(["run", str(panel), "--config", str(path)]) == 0
    assert (panel / "scripts" / "run_all.sh").is_file()


def test_run_refuses_a_directory_that_was_never_written(tmp_path):
    assert main(["run", str(tmp_path)]) == 1


# --- what a failure leaves behind -------------------------------------------

def test_a_failure_keeps_the_whole_output_not_just_the_last_line(tmp_path):
    """Re-running a GPU job to read a message it already printed is expensive."""
    from foldrunner.run import Invocation, execute

    script = "\n".join(f"echo line{i}" for i in range(20)) + "\nexit 3"
    invocation = Invocation(engine="boltz2", jobs=("A__B",), out=tmp_path / "out", script=script)
    result = execute(invocation, log_dir=tmp_path / "logs")

    assert not result.ok
    assert result.returncode == 3
    assert result.log is not None and result.log.is_file()
    body = result.log.read_text()
    assert "line0" in body and "line19" in body     # the whole run, not the tail
    assert "--- script ---" in body                  # and what produced it
    assert "exit 3" in body
    assert result.tail(3).splitlines() == ["line17", "line18", "line19"]


def test_a_run_that_produced_a_structure_is_logged_and_counted(tmp_path):
    from foldrunner.run import Invocation, execute

    out = tmp_path / "out"
    invocation = Invocation(
        engine="boltz2",
        jobs=("A__B",),
        out=out,
        script=f"mkdir -p {out} && echo data > {out}/A__B_model_0.cif && echo fine",
    )
    result = execute(invocation, log_dir=tmp_path / "logs")
    assert result.ok
    assert result.produced == 1
    assert "fine" in result.log.read_text()


def test_exiting_zero_without_writing_anything_is_not_success(tmp_path):
    """Several engines catch their own per-complex failures and still exit 0, so
    a panel can report success having written no structure at all."""
    from foldrunner.run import Invocation, execute

    invocation = Invocation(
        engine="boltz2", jobs=("A__B",), out=tmp_path / "out", script="echo nothing to do"
    )
    result = execute(invocation, log_dir=tmp_path / "logs")
    assert result.returncode == 0
    assert result.produced == 0
    assert result.silent_failure
    assert not result.ok


def test_manual_invocations_are_skipped_not_run(tmp_path):
    from foldrunner.run import Invocation, execute

    invocation = Invocation(
        engine="af_server", jobs=("A__B",), out=tmp_path, manual=True, reason="upload it"
    )
    result = execute(invocation, log_dir=tmp_path / "logs")
    assert result.skipped and result.ok
    assert result.log is None


def test_a_configured_executable_replaces_the_default(tmp_path):
    """localcolabfold installs outside PATH, so the config has to be able to
    point at it; a command template that hardcodes the name ignores that."""
    from foldrunner.engines import get_engine
    from foldrunner.ir import Component, Job, Protein
    from foldrunner.run import build_invocation, engine_executable

    assert engine_executable("colabfold") == "colabfold_batch"
    job = Job(
        "A__B",
        (Component(Protein("A", "AQVINTFDGV")), Component(Protein("B", "KKAVINGEQI"))),
    )
    written = get_engine("colabfold")(job, tmp_path / "cf")
    entry = EngineConfig(name="colabfold", executable="/opt/localcolabfold/colabfold_batch")
    invocation = build_invocation(written, entry, tmp_path / "out")
    assert "/opt/localcolabfold/colabfold_batch" in invocation.script
    assert "\ncolabfold_batch" not in invocation.script


def test_a_container_runs_as_the_invoking_user_by_default(tmp_path):
    """A container running as root leaves results their owner cannot remove."""
    from foldrunner.config import EngineConfig
    from foldrunner.engines import get_engine
    from foldrunner.ir import Component, Job, Protein
    from foldrunner.run import build_invocation

    job = Job("A__B", (Component(Protein("A", "AQVINTFDGV")),))
    written = get_engine("af3")(job, tmp_path / "af3")
    entry = EngineConfig(name="af3", container="docker", image="alphafold3")
    assert '--user "$(id -u):$(id -g)"' in build_invocation(written, entry, tmp_path).script


def test_an_image_that_needs_root_can_opt_out(tmp_path):
    from foldrunner.config import EngineConfig
    from foldrunner.engines import get_engine
    from foldrunner.ir import Component, Job, Protein
    from foldrunner.run import build_invocation

    job = Job("A__B", (Component(Protein("A", "AQVINTFDGV")),))
    written = get_engine("af3")(job, tmp_path / "af3")
    entry = EngineConfig(name="af3", container="docker", image="alphafold3", user="")
    assert "--user" not in build_invocation(written, entry, tmp_path).script


def test_an_explicit_user_is_passed_through(tmp_path):
    from foldrunner.config import EngineConfig
    from foldrunner.engines import get_engine
    from foldrunner.ir import Component, Job, Protein
    from foldrunner.run import build_invocation

    job = Job("A__B", (Component(Protein("A", "AQVINTFDGV")),))
    written = get_engine("af3")(job, tmp_path / "af3")
    entry = EngineConfig(name="af3", container="docker", image="alphafold3", user="1000:1000")
    assert "--user 1000:1000" in build_invocation(written, entry, tmp_path).script


def test_a_panel_moved_to_another_machine_runs_against_its_new_location(tmp_path):
    """Configuring a host means the panel is run somewhere else; absolute paths
    recorded at write time would point at the machine it was written on."""
    import json
    import shutil
    from pathlib import Path

    from foldrunner.cli import main
    from foldrunner.manifest import PANEL

    fasta = tmp_path / "lib.fasta"
    fasta.write_text(">A\nAQVINTFDGV\n>B\nKKAVINGEQI\n")
    here = tmp_path / "here"
    main(["write", str(fasta), "-o", str(here), "--engines", "boltz2,chai1"])

    snapshot = json.loads((here / PANEL).read_text())
    assert not any(
        Path(item["path"]).is_absolute() for item in snapshot["written"]
    ), "panel snapshot must not record absolute paths"

    # Move the panel, as copying it to another machine would.
    there = tmp_path / "there"
    shutil.move(str(here), str(there))
    assert main(["run", str(there), "-o", str(tmp_path / "out"), "--runner", "script"]) == 0

    text = "\n".join(p.read_text() for p in (there / "scripts").glob("*.sh"))
    assert str(there) in text
    assert str(here) not in text
