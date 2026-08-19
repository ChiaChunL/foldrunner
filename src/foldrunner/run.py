"""Turning written inputs into commands, and running them.

The engine modules know the shape of their invocation, the config knows where
the engine is installed, and this module puts the two together. Keeping them
apart is what lets the same panel run on a workstation and on a shared GPU box
behind a conda environment without any of the writers changing. To run it on
another machine, copy the panel there and run it there: the snapshot records
relative paths so that it survives the move.

Two engines have no command at all. AlphaFold Server and SeedFold are web forms;
their inputs are uploaded by hand. They are reported as such rather than given
an invented command line.
"""

from __future__ import annotations

import importlib
import shlex
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from foldrunner.config import Config, EngineConfig
from foldrunner.engines.base import WrittenJob

RUNNERS = ("script", "local")


@dataclass
class Invocation:
    """One command to run, or one thing a person has to do."""

    engine: str
    jobs: tuple[str, ...]
    out: Path
    script: str = ""
    manual: bool = False
    reason: str = ""

    seed: int | None = None

    @property
    def label(self) -> str:
        base = self.jobs[0] if len(self.jobs) == 1 else f"{len(self.jobs)} complexes"
        return f"{base}_seed{self.seed}" if self.seed is not None else base


@dataclass
class RunResult:
    invocation: Invocation
    returncode: int = 0
    stdout: str = ""
    stderr: str = ""
    skipped: bool = False
    log: Path | None = None
    produced: int = 0

    @property
    def ok(self) -> bool:
        """Whether the run actually did anything.

        A zero exit is not enough. Several engines catch their own per-complex
        failures and still exit cleanly, so a panel can report success having
        written no structure at all — which nobody notices until the scoring
        step finds an empty directory.
        """
        if self.skipped:
            return True
        return self.returncode == 0 and self.produced > 0

    @property
    def silent_failure(self) -> bool:
        """Exited cleanly, produced nothing."""
        return not self.skipped and self.returncode == 0 and self.produced == 0

    def tail(self, lines: int = 15) -> str:
        """The end of the combined output, which is where the cause usually is."""
        combined = (self.stdout + self.stderr).strip().splitlines()
        return "\n".join(combined[-lines:])


def engine_command(name: str) -> str | None:
    """The invocation template an engine module declares."""
    module = importlib.import_module(f"foldrunner.engines.{name}")
    return getattr(module, "COMMAND", None)


def engine_scope(name: str) -> str:
    """How much of a panel one invocation of this engine covers.

    ``job`` is one complex, ``file`` is whatever single file the writer produced,
    and ``directory`` means the engine reads the whole directory itself — for
    those, one invocation per complex would run every complex once per complex.
    """
    module = importlib.import_module(f"foldrunner.engines.{name}")
    return getattr(module, "RUN_SCOPE", "job")


def engine_executable(name: str) -> str:
    """The program an engine module expects to find on PATH by default."""
    module = importlib.import_module(f"foldrunner.engines.{name}")
    return getattr(module, "DEFAULT_EXECUTABLE", "")


def _context(written: WrittenJob, out: Path, config: EngineConfig) -> dict[str, str]:
    job = written.job
    return {
        "engine": written.engine,
        "job": job.name,
        "input": str(written.path),
        "input_dir": str(written.path.parent),
        "out": str(out),
        "seed": str(job.seeds[0]) if job.seeds else "0",
        "seeds": ",".join(str(s) for s in job.seeds) or "0",
        "msa_dir": written.extra.get("msa_directory", "-"),
        "msa_a3m_dir": written.extra.get("msa_a3m_directory", "-"),
        "script": written.extra.get("script", ""),
        "executable": config.executable or engine_executable(written.engine),
        "extra_args": config.extra_args,
        "image": config.image or "",
    }


def render_command(template: str, context: dict[str, str]) -> str:
    try:
        return " ".join(template.format(**context).split())
    except KeyError as error:
        raise KeyError(
            f"command template refers to {error}, which is not one of "
            f"{sorted(context)}"
        ) from None


def _containerise(command: str, config: EngineConfig, context: dict[str, str]) -> str:
    """Wrap a command for a container runtime.

    Mounts default to binding each referenced directory at its own path. A
    container that sees the input at a different path than the file names
    written into it cannot resolve them, and alignment paths inside an
    AlphaFold 3 job are absolute.
    """
    if not config.container:
        return command
    if not config.image:
        raise ValueError(f"[engines.{config.name}] sets container but no image")

    mounts = config.mounts or [
        f"{context['input_dir']}:{context['input_dir']}",
        f"{context['out']}:{context['out']}",
    ]
    parts = [config.container, "run", "--rm"]
    if config.gpus:
        parts += ["--gpus", config.gpus]
    # Without this the container writes as root and the results belong to a user
    # who cannot move or delete them — which only shows up after a long run.
    if config.user is None:
        parts += ["--user", '"$(id -u):$(id -g)"']
    elif config.user:
        parts += ["--user", config.user]
    for mount in mounts:
        parts += ["-v", mount.format(**context)]
    for key, value in config.env.items():
        parts += ["-e", f"{key}={value}"]
    parts.append(config.image)
    return " ".join(parts) + " " + command


def build_invocation(
    written: WrittenJob,
    config: EngineConfig,
    out_root: Path,
    jobs: tuple[str, ...] | None = None,
    seed: int | None = None,
) -> Invocation:
    """Assemble the shell snippet that runs one engine on one input."""
    name = written.engine
    template = config.command or engine_command(name)
    stem = written.job.name if jobs is None else "panel"
    if seed is not None:
        stem = f"{stem}_seed{seed}"
    out = Path(out_root) / name / stem
    labels = jobs or (written.job.name,)

    if template is None:
        return Invocation(
            engine=name,
            jobs=labels,
            out=out,
            manual=True,
            reason=(
                f"{name} is submitted through a web form; upload {written.path} "
                "and download the results into the output directory"
            ),
        )

    context = _context(written, out, config)
    if seed is not None:
        context["seed"] = str(seed)
    command = render_command(template, context)

    lines: list[str] = []
    lines += config.preamble
    if config.conda_env and not config.container:
        lines.append(f"conda activate {shlex.quote(config.conda_env)}")
    if not config.container:
        # Inside a container these are passed with -e instead.
        lines += [f"export {key}={shlex.quote(value)}" for key, value in config.env.items()]
    lines.append(f"mkdir -p {shlex.quote(str(out))}")
    lines.append(_containerise(command, config, context))

    return Invocation(engine=name, jobs=labels, out=out, script="\n".join(lines), seed=seed)


def plan(
    written: Iterable[WrittenJob],
    config: Config,
    out_root: Path,
    engines: Iterable[str] | None = None,
) -> list[Invocation]:
    """Build one invocation per unit of work.

    Engines that take a whole panel in one file produce one invocation for it,
    not one per complex.
    """
    wanted = set(engines) if engines else None
    grouped: dict[tuple[str, str], list[WrittenJob]] = {}
    for item in written:
        if wanted and item.engine not in wanted:
            continue
        scope = engine_scope(item.engine)
        if scope == "directory":
            key = str(item.path.parent)
        elif scope in ("file", "manual"):
            # Both are one file covering many complexes: a batched panel, or a
            # single upload.
            key = str(item.path)
        else:
            key = f"{item.path}::{item.job.name}"
        grouped.setdefault((item.engine, key), []).append(item)

    invocations = []
    for (name, _), items in grouped.items():
        entry = config.for_engine(name)
        if not entry.enabled:
            continue
        batched = len(items) > 1
        labels = tuple(dict.fromkeys(i.job.name for i in items))
        template = entry.command or engine_command(name)
        seeds = items[0].job.seeds

        # An engine whose template takes a single {seed} runs once per seed, to
        # its own output directory. Passing only the first would silently drop
        # the rest, and the output carries no seed in its file names to notice
        # it afterwards.
        if template and "{seed}" in template and len(seeds) > 1:
            for seed in seeds:
                invocations.append(
                    build_invocation(
                        items[0],
                        entry,
                        out_root,
                        jobs=labels if batched else None,
                        seed=seed,
                    )
                )
            continue

        invocations.append(
            build_invocation(items[0], entry, out_root, jobs=labels if batched else None)
        )
    return invocations


def write_scripts(invocations: Iterable[Invocation], directory: Path) -> list[Path]:
    """Write one shell script per invocation, plus a runner that calls them all."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    runnable = [i for i in invocations if not i.manual]

    for index, invocation in enumerate(runnable):
        path = directory / f"{index:03d}_{invocation.engine}_{invocation.label}.sh".replace(
            " ", "_"
        )
        path.write_text(f"#!/usr/bin/env bash\nset -euo pipefail\n\n{invocation.script}\n")
        path.chmod(0o755)
        written.append(path)

    if written:
        index_path = directory / "run_all.sh"
        body = "\n".join(f'"$HERE/{p.name}"' for p in written)
        index_path.write_text(
            "#!/usr/bin/env bash\n"
            "set -euo pipefail\n"
            'HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"\n\n'
            f"{body}\n"
        )
        index_path.chmod(0o755)
        written.append(index_path)

    manual = [i for i in invocations if i.manual]
    if manual:
        notes = "\n".join(f"- {i.engine}: {i.reason}" for i in manual)
        (directory / "MANUAL.md").write_text(
            "# Submitted by hand\n\n"
            "These engines are web forms and have no command line.\n\n"
            f"{notes}\n"
        )
        written.append(directory / "MANUAL.md")
    return written


def _decode(stream) -> str:
    if stream is None:
        return ""
    return stream.decode() if isinstance(stream, bytes) else stream


# Extensions an engine writes a structure to.
STRUCTURE_SUFFIXES = (".cif", ".pdb")


def _count_outputs(directory: Path) -> int:
    root = Path(directory)
    if not root.is_dir():
        return 0
    return sum(1 for p in root.rglob("*") if p.suffix.lower() in STRUCTURE_SUFFIXES)


def execute(
    invocation: Invocation,
    timeout: float | None = None,
    log_dir: Path | None = None,
):
    """Run one invocation on this machine.

    There is deliberately no remote mode. Driving an engine over ssh needs the
    panel visible at the same path on both the machine reading it and the
    machine running it, and sites draw that line differently. Installing the
    package where the engine already lives and running it there is one rsync
    and no new failure modes — the panel records relative paths, so it moves.

    The whole output is kept, not just the tail. An engine that fails usually
    explains itself several lines before the exception, and re-running a GPU job
    to read the message it already printed is an expensive way to learn nothing
    new.
    """
    if invocation.manual:
        return RunResult(invocation=invocation, skipped=True)

    before = _count_outputs(invocation.out)

    command = ["bash", "-c", invocation.script]

    try:
        process = subprocess.run(
            command, capture_output=True, text=True, timeout=timeout, check=False
        )
        returncode, stdout, stderr = process.returncode, process.stdout, process.stderr
    except subprocess.TimeoutExpired as expired:
        returncode = 124
        stdout = _decode(expired.stdout)
        stderr = _decode(expired.stderr) + f"\nfoldrunner: timed out after {timeout}s"

    log = None
    if log_dir is not None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)
        log = log_dir / f"{invocation.engine}_{invocation.label}.log".replace(" ", "_")
        log.write_text(
            f"# {' '.join(command[:2])}\n"
            f"# exit {returncode}\n\n"
            f"--- script ---\n{invocation.script}\n\n"
            f"--- stdout ---\n{stdout}\n"
            f"--- stderr ---\n{stderr}\n"
        )

    return RunResult(
        invocation=invocation,
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
        log=log,
        produced=_count_outputs(invocation.out) - before,
    )
