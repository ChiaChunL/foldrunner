"""The local backend is tested through a stub search program."""

import os
import stat

import pytest

from foldrunner.msa.backends import SearchError
from foldrunner.msa.backends.local import LocalBackend

SEQ_A = "MKVLAAGIVG"
SEQ_B = "PQRSTVWYAC"

STUB = """#!/bin/sh
# arguments: query database out ...
out="$3"
for header in $(grep '^>' "$1" | tr -d '>'); do
  {
    echo ">query"
    echo "MKVLAAGIVG"
    echo ">UniRef100_X_9606/"
    echo "MKVLAAGIVG"
  } > "$out/${header}__SUFFIX__"
done
"""


@pytest.fixture
def stub_search(tmp_path, monkeypatch):
    def make(suffix=".a3m", exit_code=0):
        script = tmp_path / "colabfold_search_stub"
        body = STUB.replace("__SUFFIX__", suffix)
        if exit_code:
            body += f"exit {exit_code}\n"
        script.write_text(body)
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
        return script

    return make


@pytest.fixture
def database(tmp_path):
    directory = tmp_path / "db"
    directory.mkdir()
    return directory


def test_a_missing_database_directory_says_how_big_it_needs_to_be(tmp_path):
    with pytest.raises(SearchError, match="940 GB"):
        LocalBackend(database_dir=tmp_path / "absent")


def test_a_missing_executable_is_reported(database):
    with pytest.raises(SearchError, match="not on PATH"):
        LocalBackend(database_dir=database, executable="definitely-not-installed")


def test_one_search_covers_the_whole_batch(database, stub_search):
    script = stub_search()
    backend = LocalBackend(database_dir=database, executable=script.name)
    results = backend.fetch([SEQ_A, SEQ_B])
    assert len(results) == 2
    assert all(r.unpaired for r in results)


def test_repeated_sequences_are_searched_once(database, stub_search):
    script = stub_search()
    backend = LocalBackend(database_dir=database, executable=script.name)
    results = backend.fetch([SEQ_A, SEQ_B, SEQ_A])
    assert [r.sequence for r in results] == [SEQ_A, SEQ_B, SEQ_A]


def test_pairing_reads_the_paired_output(database, stub_search):
    script = stub_search(suffix=".paired.a3m")
    backend = LocalBackend(database_dir=database, executable=script.name)
    results = backend.fetch([SEQ_A], paired=True)
    assert results[0].paired is not None
    assert results[0].unpaired is None


def test_a_failing_search_surfaces_its_error(database, stub_search):
    script = stub_search(exit_code=3)
    backend = LocalBackend(database_dir=database, executable=script.name)
    with pytest.raises(SearchError, match="exited 3"):
        backend.fetch([SEQ_A])


def test_missing_output_names_what_was_produced(database, stub_search):
    script = stub_search(suffix=".unexpected")
    backend = LocalBackend(database_dir=database, executable=script.name)
    with pytest.raises(SearchError, match="no output for q0"):
        backend.fetch([SEQ_A])


def test_pairing_strategy_is_passed_through(database, stub_search):
    script = stub_search()
    backend = LocalBackend(database_dir=database, executable=script.name)
    command = backend._command("q.fasta", "out", paired=True)
    assert "--pairing_strategy" in command
    assert command[command.index("--pairing_strategy") + 1] == "0"
    strict = LocalBackend(
        database_dir=database, executable=script.name, pairing_strategy="complete"
    )
    command = strict._command("q.fasta", "out", paired=True)
    assert command[command.index("--pairing_strategy") + 1] == "1"


def test_local_and_remote_alignments_are_kept_apart(database, stub_search):
    from foldrunner.msa.backends import WebBackend

    script = stub_search()
    local = LocalBackend(database_dir=database, executable=script.name).recipe
    remote = WebBackend(contact="a@b.c").recipe
    assert local.backend == "local-mmseqs"
    assert local.recipe_id != remote.recipe_id


def test_database_version_is_part_of_the_recipe(database, stub_search):
    script = stub_search()
    one = LocalBackend(database_dir=database, executable=script.name, database_version="2024")
    two = LocalBackend(database_dir=database, executable=script.name, database_version="2026")
    assert one.recipe.recipe_id != two.recipe.recipe_id


def test_the_cli_refuses_local_search_rather_than_running_it_unverified(tmp_path, capsys):
    """Groundwork, not a feature: an unverified search sitting between the panel
    and every score computed from it is the worst place for a silent difference."""
    from foldrunner.cli import main

    fasta = tmp_path / "lib.fasta"
    fasta.write_text(">A\nMKVLAAGIVG\n")
    assert main(["search", str(fasta), "--cache", str(tmp_path / "msa"),
                 "--local-db", str(tmp_path)]) == 2
    message = capsys.readouterr().err
    assert "not available yet" in message
    assert "foldrunner import" in message
