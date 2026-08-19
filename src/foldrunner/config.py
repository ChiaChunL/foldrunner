"""Where the engines are installed, and how to enter their environments.

Every engine on a real machine sits behind something: a conda environment that a
non-interactive shell does not have on its PATH, a container that needs the input
directory mounted at the same path it has on the host, an environment variable
without which a model download goes to a host that is unreachable, or another
machine entirely. None of that belongs in the writers, and none of it is
portable, so it lives in a config file that travels with the site rather than
with the package.

The engine modules own the *shape* of their command; this file owns *where*.
"""

from __future__ import annotations

import os
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

CONFIG_NAME = "foldrunner.toml"
ENV_VAR = "FOLDRUNNER_CONFIG"

# Searched in order; the first that exists wins.
SEARCH_PATHS = (
    Path.cwd() / CONFIG_NAME,
    Path.home() / ".config" / "foldrunner" / CONFIG_NAME,
    Path.home() / f".{CONFIG_NAME}",
)


class ConfigError(RuntimeError):
    """The configuration cannot be used as written."""


def _load_toml(path: Path) -> dict[str, Any]:
    try:
        import tomllib
    except ModuleNotFoundError:  # Python 3.10
        try:
            import tomli as tomllib  # type: ignore[no-redef]
        except ModuleNotFoundError as error:
            raise ConfigError(
                "reading a config file on Python 3.10 needs tomli: "
                'pip install "foldrunner[config]"'
            ) from error
    with path.open("rb") as handle:
        return tomllib.load(handle)


@dataclass
class EngineConfig:
    """How to invoke one engine on this machine."""

    name: str
    command: str | None = None
    conda_env: str | None = None
    executable: str | None = None
    env: dict[str, str] = field(default_factory=dict)
    preamble: list[str] = field(default_factory=list)
    container: str | None = None
    image: str | None = None
    mounts: list[str] = field(default_factory=list)
    gpus: str | None = None
    # Who the container runs as. Unset means the invoking user, because a
    # container that runs as root leaves results its owner cannot delete or
    # move. Set it to "" to let the image decide, or to an explicit uid:gid.
    user: str | None = None
    extra_args: str = ""
    enabled: bool = True

    def describe(self) -> str:
        where = []
        if self.conda_env:
            where.append(f"conda {self.conda_env}")
        if self.container:
            where.append(f"{self.container} {self.image or '?'}")
        if self.executable:
            where.append(self.executable)
        return ", ".join(where) if where else "on PATH"


@dataclass
class Config:
    """Site configuration: defaults plus one entry per engine."""

    preamble: list[str] = field(default_factory=list)
    env: dict[str, str] = field(default_factory=dict)
    conda_sh: str | None = None
    output_root: str | None = None
    engines: dict[str, EngineConfig] = field(default_factory=dict)
    path: Path | None = None

    def for_engine(self, name: str) -> EngineConfig:
        """Config for one engine, with the defaults folded in.

        Defaults come first so an engine can override a variable the site sets
        globally, and the shared preamble runs before the engine's own.
        """
        specific = self.engines.get(name, EngineConfig(name=name))
        preamble = list(self.preamble)
        if self.conda_sh:
            # A non-interactive shell is not a login shell, so conda is not on
            # PATH and `conda activate` is not even defined until this is run.
            line = f"source {shlex.quote(self.conda_sh)}"
            if line not in preamble:
                preamble.insert(0, line)
        preamble += [p for p in specific.preamble if p not in preamble]
        return EngineConfig(
            name=name,
            command=specific.command,
            conda_env=specific.conda_env,
            executable=specific.executable,
            env={**self.env, **specific.env},
            preamble=preamble,
            container=specific.container,
            image=specific.image,
            mounts=list(specific.mounts),
            gpus=specific.gpus,
            user=specific.user,
            extra_args=specific.extra_args,
            enabled=specific.enabled,
        )

    def configured(self) -> list[str]:
        return sorted(name for name, entry in self.engines.items() if entry.enabled)


def find_config(explicit: str | Path | None = None) -> Path | None:
    """Locate a config file: explicit path, environment variable, then defaults."""
    if explicit:
        path = Path(explicit)
        if not path.is_file():
            raise ConfigError(f"config file {path} does not exist")
        return path
    from_env = os.environ.get(ENV_VAR)
    if from_env:
        path = Path(from_env)
        if not path.is_file():
            raise ConfigError(f"{ENV_VAR} points at {path}, which does not exist")
        return path
    for candidate in SEARCH_PATHS:
        if candidate.is_file():
            return candidate
    return None


def load(explicit: str | Path | None = None) -> Config:
    """Read the config, or return an empty one when there is none.

    An absent config is not an error: writing inputs needs none, and running
    without one still works for engines that are already on PATH.
    """
    path = find_config(explicit)
    if path is None:
        return Config()

    raw = _load_toml(path)
    unknown = set(raw) - {"defaults", "engines"}
    if unknown:
        raise ConfigError(
            f"{path}: unknown top-level section(s) {sorted(unknown)}; "
            "expected [defaults] and [engines.<name>]"
        )

    defaults = raw.get("defaults", {})
    config = Config(
        preamble=list(defaults.get("preamble", [])),
        env={str(k): str(v) for k, v in defaults.get("env", {}).items()},
        conda_sh=defaults.get("conda_sh"),
        output_root=defaults.get("output_root"),
        path=path,
    )

    known = {f.name for f in EngineConfig.__dataclass_fields__.values()} - {"name"}
    for name, entry in raw.get("engines", {}).items():
        if not isinstance(entry, dict):
            raise ConfigError(f"{path}: [engines.{name}] must be a table")
        bad = set(entry) - known
        if bad:
            raise ConfigError(
                f"{path}: [engines.{name}] has unknown key(s) {sorted(bad)}; "
                f"known keys are {sorted(known)}"
            )
        config.engines[name] = EngineConfig(
            name=name,
            command=entry.get("command"),
            conda_env=entry.get("conda_env"),
            executable=entry.get("executable"),
            env={str(k): str(v) for k, v in entry.get("env", {}).items()},
            preamble=list(entry.get("preamble", [])),
            container=entry.get("container"),
            image=entry.get("image"),
            mounts=list(entry.get("mounts", [])),
            gpus=entry.get("gpus"),
            user=entry.get("user"),
            extra_args=entry.get("extra_args", ""),
            enabled=bool(entry.get("enabled", True)),
        )
    return config
