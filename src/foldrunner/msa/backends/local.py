"""Local MMseqs2 search through colabfold_search.

**Not available yet.** The code below is exercised against a stub, never against
real databases, so `foldrunner search --local-db` refuses to use it. Until it has
been run for real, the supported route for a local search is to run
``colabfold_search`` directly and adopt the output with ``foldrunner import``,
which is a path that has been exercised at proteome scale.

What follows is the groundwork for lifting that restriction.

Worth the setup once a panel outgrows the public services: those are shared
academic resources sized for a few thousand alignments a day across all users,
and a proteome-scale library exceeds that on its own.

The search is run once over the whole batch. mmseqs scans the databases per run,
not per query, so a loop over sequences repeats the expensive part for every one
of them; handing it the full set is the difference between one pass and N.

Database and tool versions are recorded in the recipe. Alignments produced here
are not interchangeable with the services' — depth differs with the database
release — and depth moves the confidence scores, so the two are kept apart in
the cache rather than pooled.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from foldrunner.msa.backends.base import MsaResult, SearchError
from foldrunner.msa.cache import Recipe

# colabfold_search writes one file per query, named after the FASTA header.
UNPAIRED_SUFFIX = ".a3m"
PAIRED_SUFFIXES = (".paired.a3m", ".env.paired.a3m")

# Setting up the databases takes roughly this much disk; the number is worth
# surfacing before someone starts a download that will not fit.
DATABASE_FOOTPRINT_GB = 940


@dataclass
class LocalBackend:
    """Runs colabfold_search against locally installed databases."""

    database_dir: Path
    executable: str = "colabfold_search"
    mmseqs: str | None = None
    threads: int = 0
    use_env: bool = True
    pairing_strategy: str = "greedy"
    db_load_mode: int = 0
    tool_version: str = ""
    database_version: str = ""
    extra_args: tuple[str, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        self.database_dir = Path(self.database_dir)
        if not self.database_dir.is_dir():
            raise SearchError(
                f"database directory {self.database_dir} does not exist. "
                f"Set it up with colabfold's setup_databases.sh; it needs about "
                f"{DATABASE_FOOTPRINT_GB} GB."
            )
        if shutil.which(self.executable) is None:
            raise SearchError(f"{self.executable} is not on PATH")

    @property
    def recipe(self) -> Recipe:
        databases = ("uniref30", "colabfold_envdb") if self.use_env else ("uniref30",)
        return Recipe(
            backend="local-mmseqs",
            databases=databases,
            # Matching the service's mmseqs release is what makes local and
            # remote alignments comparable; a mismatch changes the depth.
            tool_version=self.tool_version or "unspecified",
            params=(
                ("pairing_strategy", self.pairing_strategy),
                ("database_version", self.database_version or "unspecified"),
            ),
        )

    def _command(self, query: Path, out: Path, paired: bool) -> list[str]:
        command = [self.executable, str(query), str(self.database_dir), str(out)]
        if self.mmseqs:
            command += ["--mmseqs", self.mmseqs]
        if self.threads:
            command += ["--threads", str(self.threads)]
        command += ["--db-load-mode", str(self.db_load_mode)]
        if paired:
            command += ["--pairing_strategy", "0" if self.pairing_strategy == "greedy" else "1"]
            if self.use_env:
                command += ["--use-env-pairing", "1"]
        else:
            command += ["--use-env", "1" if self.use_env else "0"]
        return command + list(self.extra_args)

    def fetch(self, sequences: list[str], *, paired: bool = False) -> list[MsaResult]:
        unique: list[str] = []
        for sequence in sequences:
            if sequence not in unique:
                unique.append(sequence)

        with tempfile.TemporaryDirectory() as workspace:
            root = Path(workspace)
            query = root / "query.fasta"
            # Headers are positional here only because the file is written and
            # consumed inside this call; nothing positional escapes it.
            query.write_text(
                "".join(f">q{index}\n{sequence}\n" for index, sequence in enumerate(unique))
            )
            out = root / "out"
            out.mkdir()

            command = self._command(query, out, paired)
            process = subprocess.run(command, capture_output=True, text=True, check=False)
            if process.returncode != 0:
                raise SearchError(
                    f"{self.executable} exited {process.returncode}: "
                    f"{process.stderr.strip()[-500:]}"
                )

            results: dict[str, MsaResult] = {}
            for index, sequence in enumerate(unique):
                text = self._read_output(out, f"q{index}", paired)
                result = MsaResult(sequence=sequence)
                if paired:
                    result.paired = text
                else:
                    result.unpaired = text
                results[sequence] = result
        return [results[sequence] for sequence in sequences]

    def _read_output(self, out: Path, stem: str, paired: bool) -> str:
        suffixes = PAIRED_SUFFIXES if paired else (UNPAIRED_SUFFIX,)
        for suffix in suffixes:
            candidate = out / f"{stem}{suffix}"
            if candidate.is_file():
                return candidate.read_text()
        available = sorted(p.name for p in out.iterdir())[:5]
        raise SearchError(
            f"no output for {stem} in {out}; found {available}. "
            "Check that the database directory matches the search mode."
        )
