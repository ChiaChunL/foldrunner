import re

import pytest

from foldrunner.cli import main
from foldrunner.engines import get_engine
from foldrunner.ir import Component, Job, Protein
from foldrunner.manifest import completed_jobs, read_manifest, rows_for, write_manifest

A = Protein("A", "AQVINTFDGV")
B = Protein("B", "KKAVINGEQI")


def test_manifest_records_the_chain_each_engine_assigned(tmp_path):
    job = Job("A__B", (Component(A), Component(B)))
    rows = []
    for engine in ("boltz2", "chai1"):
        rows.extend(rows_for(get_engine(engine)(job, tmp_path / engine)))
    path = write_manifest(tmp_path / "manifest.tsv", rows)
    table = read_manifest(path)
    assert {r["engine"] for r in table} == {"boltz2", "chai1"}
    assert [r["chain"] for r in table if r["engine"] == "boltz2"] == ["A", "B"]


def test_homodimer_contributes_one_row_per_copy(tmp_path):
    job = Job("A__A", (Component(A, 2),))
    rows = rows_for(get_engine("boltz2")(job, tmp_path))
    assert len(rows) == 2
    assert [r.chain for r in rows] == ["A", "B"]
    assert {r.entity for r in rows} == {"A"}


def test_paired_depth_travels_with_the_job(tmp_path):
    job = Job("A__B", (Component(A), Component(B)))
    rows = rows_for(get_engine("boltz2")(job, tmp_path), n_paired=1681)
    assert all(r.n_paired == 1681 for r in rows)


def test_completed_jobs_supports_resuming(tmp_path):
    jobs = [Job("A__B", (Component(A), Component(B))), Job("A__A", (Component(A, 2),))]
    rows = [r for j in jobs for r in rows_for(get_engine("boltz2")(j, tmp_path))]
    path = write_manifest(tmp_path / "manifest.tsv", rows)
    assert completed_jobs(path, "boltz2") == {"A__B", "A__A"}
    assert completed_jobs(path, "chai1") == set()


def test_completed_jobs_on_a_fresh_run_is_empty(tmp_path):
    assert completed_jobs(tmp_path / "absent.tsv", "boltz2") == set()


def test_cli_plan_reports_alignment_cost(tmp_path, capsys):
    fasta = tmp_path / "lib.fasta"
    fasta.write_text(">A\nAQVINTFDGV\n>B\nKKAVINGEQI\n>C\nMEEPQSDPSV\n")
    assert main(["plan", str(fasta)]) == 0
    out = capsys.readouterr().out
    assert "complexes      6" in out
    assert "MSA searches   3" in out


def test_cli_write_produces_inputs_and_a_manifest(tmp_path):
    fasta = tmp_path / "lib.fasta"
    fasta.write_text(">A\nAQVINTFDGV\n>B\nKKAVINGEQI\n")
    out = tmp_path / "out"
    assert main(["write", str(fasta), "-o", str(out), "--engines", "boltz2,chai1"]) == 0
    assert len(list((out / "boltz2").glob("*.yaml"))) == 3
    assert len(read_manifest(out / "manifest.tsv")) == 12


@pytest.mark.parametrize(
    "argv",
    [
        ["plan", "LIB"],
        ["search", "LIB", "--cache", "C", "--contact", "a@b.c"],
        ["search", "LIB", "--cache", "C", "--local-db", "D"],
        ["import", "S", "--cache", "C"],
        ["write", "LIB", "-o", "O"],
    ],
)
def test_every_subcommand_supplies_the_arguments_its_handler_reads(argv, tmp_path):
    """Guards against a parser option being dropped while the handler still reads it."""
    import inspect

    from foldrunner.cli import build_parser

    args = build_parser().parse_args(argv)
    source = inspect.getsource(args.func)
    for name in set(re.findall(r"args\.([a-z_]+)", source)):
        assert hasattr(args, name), f"{argv[0]} reads args.{name} but does not define it"
