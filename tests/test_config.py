"""Site configuration: where the engines are installed."""

import pytest

from foldrunner.config import ENV_VAR, Config, ConfigError, EngineConfig, find_config, load

SAMPLE = """
[defaults]
conda_sh = "/opt/conda/etc/profile.d/conda.sh"
output_root = "/scratch/results"
env = { SHARED = "1" }

[engines.boltz2]
conda_env = "boltz2"

[engines.chai1]
conda_env = "chai1"
env = { HF_ENDPOINT = "https://mirror.example", SHARED = "2" }

[engines.af3]
container = "docker"
image = "alphafold3"
gpus = "all"
mounts = ["{input_dir}:{input_dir}"]

[engines.af2_multimer]
executable = "/opt/af2/run.sh"

[engines.colabfold]
enabled = false
"""


@pytest.fixture
def config_file(tmp_path):
    path = tmp_path / "foldrunner.toml"
    path.write_text(SAMPLE)
    return path


def test_absent_config_is_not_an_error(tmp_path, monkeypatch):
    monkeypatch.delenv(ENV_VAR, raising=False)
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("foldrunner.config.SEARCH_PATHS", (tmp_path / "nope.toml",))
    conf = load()
    assert isinstance(conf, Config)
    assert conf.engines == {}


def test_explicit_missing_path_is_an_error(tmp_path):
    with pytest.raises(ConfigError, match="does not exist"):
        find_config(tmp_path / "absent.toml")


def test_environment_variable_is_honoured(config_file, monkeypatch):
    monkeypatch.setenv(ENV_VAR, str(config_file))
    assert find_config() == config_file


def test_engine_specific_env_overrides_the_shared_one(config_file):
    entry = load(config_file).for_engine("chai1")
    assert entry.env["SHARED"] == "2"
    assert entry.env["HF_ENDPOINT"] == "https://mirror.example"


def test_engines_inherit_the_shared_env(config_file):
    assert load(config_file).for_engine("boltz2").env["SHARED"] == "1"


def test_conda_is_sourced_before_anything_else(config_file):
    """A non-interactive shell has no conda on PATH, so this has to come first."""
    entry = load(config_file).for_engine("boltz2")
    assert entry.preamble[0] == "source /opt/conda/etc/profile.d/conda.sh"


def test_an_unconfigured_engine_still_resolves(config_file):
    entry = load(config_file).for_engine("protenix")
    assert entry.enabled
    assert entry.conda_env is None
    assert entry.preamble  # the shared preamble still applies


def test_disabled_engines_are_reported(config_file):
    conf = load(config_file)
    assert not conf.for_engine("colabfold").enabled
    assert "colabfold" not in conf.configured()


def test_unknown_engine_key_is_rejected(tmp_path):
    path = tmp_path / "foldrunner.toml"
    path.write_text('[engines.boltz2]\nconda_environment = "typo"\n')
    with pytest.raises(ConfigError, match="unknown key"):
        load(path)


def test_unknown_top_level_section_is_rejected(tmp_path):
    path = tmp_path / "foldrunner.toml"
    path.write_text('[enginez.boltz2]\nconda_env = "x"\n')
    with pytest.raises(ConfigError, match="unknown top-level"):
        load(path)


def test_describe_says_where_the_engine_lives(config_file):
    conf = load(config_file)
    assert "docker alphafold3" in conf.for_engine("af3").describe()
    assert "conda boltz2" in conf.for_engine("boltz2").describe()
    assert EngineConfig(name="x").describe() == "on PATH"


def test_a_host_key_is_rejected_rather_than_silently_ignored(tmp_path):
    """Remote execution was removed; a config carrying it should say so instead
    of running everything locally without a word."""
    config = tmp_path / "foldrunner.toml"
    config.write_text('[engines.boltz2]\nhost = "gpu-box"\n')
    with pytest.raises(ConfigError, match="unknown key"):
        load(config)


def test_the_shipped_example_config_parses():
    conf = load("examples/foldrunner.toml")
    assert conf.conda_sh
    assert "boltz2" in conf.engines
    # The web-form engines have nothing to configure.
    assert "af_server" not in conf.engines
    assert "seedfold" not in conf.engines
