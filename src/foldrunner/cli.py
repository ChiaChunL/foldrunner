"""Command line interface."""

from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

from foldrunner import __version__
from foldrunner import config as config_module
from foldrunner.engines import ENGINES, colabfold, get_engine
from foldrunner.engines.base import MsaView
from foldrunner.enumerate import (
    dedupe,
    enumerate_jobs,
    load_fasta,
    unique_sequences,
)
from foldrunner.manifest import (
    PANEL,
    Row,
    completed_jobs,
    read_manifest,
    read_panel,
    rebase,
    rows_for,
    write_manifest,
    write_panel,
)
from foldrunner.msa.backends import (
    COLABFOLD_HOST,
    PROTENIX_HOST,
    RateLimited,
    SearchError,
    WebBackend,
)
from foldrunner.msa.cache import MsaCache, Recipe
from foldrunner.msa.importer import import_tree
from foldrunner.run import RUNNERS, engine_command, execute, plan, write_scripts


def _coerce(entry: dict[str, str]) -> dict[str, object]:
    """Turn a manifest row back into the types :class:`Row` holds."""
    return {
        **entry,
        "component": int(entry["component"]),
        "count": int(entry["count"]),
        "n_chains": int(entry["n_chains"]),
        "n_residues": int(entry["n_residues"]),
        "msa_used": entry["msa_used"] == "1",
        "n_paired": int(entry["n_paired"]),
    }

# Above this many unique sequences the public MSA services stop being a
# reasonable choice: they are shared academic resources sized for a few thousand
# alignments a day across all users, and large panels get throttled.
API_SEQUENCE_WARNING = 50


class ConfigProblem(RuntimeError):
    """Something the caller has to decide before the run can continue."""


def _resolve_recipe(store: MsaCache, args: argparse.Namespace) -> Recipe:
    """Work out which alignments in the cache to use.

    A recipe id is a hash of every field, so retyping a few flags is not enough
    to reproduce one — `search` records the recipe beside the alignments and
    this reads it back. Flags remain available for a cache built elsewhere.
    """
    stored = store.stored_recipes()

    chosen = getattr(args, "recipe", None)
    if chosen:
        found = store.read_recipe(chosen)
        if found is None:
            listing = (
                "\n".join(f"    {r.recipe_id}  {r.describe()}" for r in stored)
                or "    (none)"
            )
            raise ConfigProblem(f"no recipe {chosen} in {store.root}. It holds:\n{listing}")
        return found

    asked = bool(args.backend != "colabfold-api" or args.databases or args.tool_version)

    if not asked:
        if len(stored) == 1:
            return stored[0]
        if not stored:
            raise ConfigProblem(
                f"{store.root} holds no recipe. Run `foldrunner search` against it, "
                "or name the alignments with --backend/--databases/--tool-version."
            )
        listing = "\n".join(f"    {r.recipe_id}  {r.describe()}" for r in stored)
        raise ConfigProblem(
            f"{store.root} holds {len(stored)} sets of alignments; pick one with "
            f"--recipe:\n{listing}"
        )

    wanted = _recipe(args)
    if store.keys(wanted):
        return wanted
    known = "\n".join(f"    {r.recipe_id}  {r.describe()}" for r in stored) or "    (none)"
    raise ConfigProblem(
        f"no alignments under {wanted.recipe_id} in {store.root}.\n"
        f"  The cache holds:\n{known}\n"
        "  Recipe ids cover the search parameters too, so flags alone may not "
        "reproduce one; omit them to use what the cache recorded."
    )


def _recipe(args: argparse.Namespace) -> Recipe:
    return Recipe(
        backend=args.backend,
        databases=tuple(args.databases.split(",")) if args.databases else (),
        tool_version=args.tool_version,
    )


def _load_library(path: Path) -> list:
    proteins = load_fasta(path)
    kept, duplicates = dedupe(proteins)
    for duplicate in duplicates:
        print(
            f"note: {duplicate.dropped} has the same sequence as {duplicate.kept}; "
            "keeping the first",
            file=sys.stderr,
        )
    return kept


def cmd_plan(args: argparse.Namespace) -> int:
    """Report what a run would cost before anything is submitted."""
    library = _load_library(args.fasta)
    preys = _load_library(args.preys) if args.preys else None
    jobs = enumerate_jobs(library, args.mode, preys=preys)
    sequences = unique_sequences(jobs)

    print(f"library        {len(library)} unique sequences")
    if preys:
        print(f"preys          {len(preys)} unique sequences")
    print(f"mode           {args.mode}")
    print(f"complexes      {len(jobs)}")
    print(f"MSA searches   {len(sequences)}   (one per unique sequence, not per complex)")
    if len(jobs) > len(sequences):
        print(f"               searching per complex instead would cost {len(jobs)}")
    print(f"engine files   {len(jobs) * max(len(args.engines.split(',')), 1)}")

    if len(sequences) > API_SEQUENCE_WARNING and args.backend.endswith("-api"):
        print(
            f"\nwarning: {len(sequences)} sequences exceeds what the public MSA services\n"
            f"         are meant to absorb (threshold {API_SEQUENCE_WARNING}). They are shared\n"
            "         academic resources and will throttle a panel this size.\n"
            "         Use a local search backend, or import alignments you already have.",
            file=sys.stderr,
        )
    if args.verbose:
        print()
        for job in jobs:
            print(f"  {job.name}\t{job.n_chains} chains\t{job.n_residues} aa")
    return 0


def cmd_import(args: argparse.Namespace) -> int:
    """Re-address alignments from an earlier run into the cache."""
    cache = MsaCache(args.cache)
    recipe = _recipe(args)
    print(f"recipe  {recipe.recipe_id}  ({recipe.describe()})")
    print(f"source  {args.source}")
    if args.dry_run:
        print("mode    dry run, nothing will be written")
    elif args.link:
        print("mode    hard links, no bytes copied")

    report = import_tree(
        args.source,
        cache,
        recipe,
        layout=args.layout,
        dry_run=args.dry_run,
        link=args.link,
        limit=args.limit,
    )
    print()
    print(report.summary())
    if args.show_failures:
        for directory, reason in report.failed:
            print(f"  failed: {directory} — {reason}", file=sys.stderr)
    return 0


BATCH_WRITERS = {
    "protenix": "foldrunner.engines.protenix",
    "af_server": "foldrunner.engines.af_server",
    "seedfold": "foldrunner.engines.seedfold",
}


def cmd_search(args: argparse.Namespace) -> int:
    """Compute the alignments a panel needs, one per unique sequence."""
    library = _load_library(args.fasta)
    preys = _load_library(args.preys) if args.preys else None
    jobs = enumerate_jobs(library, args.mode, preys=preys)
    sequences = unique_sequences(jobs)

    if args.local_db:
        # The local search is implemented but has only ever been exercised
        # against a stub. Offering it as though it worked would put an unverified
        # search between the panel and every score computed from it, which is the
        # one place a silent difference does the most damage.
        print(
            "error: --local-db is not available yet.\n"
            "  The local search is written but has never been run against real\n"
            "  databases, so it is not offered as though it were.\n"
            "  Run colabfold_search yourself and adopt the result instead:\n"
            "    colabfold_search library.fasta /path/to/db out/\n"
            "    foldrunner import out/ --cache msa/ --backend local-mmseqs \\\n"
            "        --databases uniref30,colabfold_envdb --tool-version <version>",
            file=sys.stderr,
        )
        return 2

    if not args.contact:
        print(
            "error: --contact is required: the alignment services ask callers to "
            "identify themselves.",
            file=sys.stderr,
        )
        return 1
    backend = WebBackend(
        host=args.host,
        fallback=args.fallback or None,
        contact=args.contact,
        use_env=not args.no_env,
        pairing_strategy=args.pairing_strategy,
    )
    cache = MsaCache(args.cache)
    recipe = backend.recipe
    cache.write_recipe(recipe)
    print(f"recipe   {recipe.recipe_id}  ({recipe.describe()})")

    kinds = ["unpaired"] + (["paired"] if not args.unpaired_only else [])
    for kind in kinds:
        filename = "paired.a3m" if kind == "paired" else "unpaired.a3m"
        todo = [
            protein
            for key, protein in sequences.items()
            if not cache.has(recipe, key, filename)
        ]
        print(f"\n{kind}: {len(sequences) - len(todo)} cached, {len(todo)} to search")
        if not todo:
            continue
        if len(sequences) > API_SEQUENCE_WARNING and not args.local_db:
            print(
                f"warning: {len(sequences)} sequences is more than the public services\n"
                "         are sized for; expect throttling and consider a local search.",
                file=sys.stderr,
            )
        for start in range(0, len(todo), args.batch_size):
            batch = todo[start : start + args.batch_size]
            try:
                results = backend.fetch(
                    [p.sequence for p in batch], paired=(kind == "paired")
                )
            except RateLimited as error:
                print(f"\nthrottled: {error}", file=sys.stderr)
                return 2
            except SearchError as error:
                print(f"\nsearch failed: {error}", file=sys.stderr)
                return 1
            for protein, result in zip(batch, results, strict=True):
                cache.store(
                    recipe,
                    protein.key,
                    protein.sequence,
                    unpaired=result.unpaired,
                    paired=result.paired,
                    source=args.host,
                )
            print(f"  {min(start + len(batch), len(todo))}/{len(todo)}")
    return 0


def cmd_write(args: argparse.Namespace) -> int:
    """Generate engine inputs for every complex in the panel."""
    library = _load_library(args.fasta)
    preys = _load_library(args.preys) if args.preys else None
    jobs = enumerate_jobs(
        library,
        args.mode,
        preys=preys,
        seeds=tuple(int(s) for s in args.seeds.split(",")) if args.seeds else (),
    )

    view = None
    if args.cache:
        store = MsaCache(args.cache)
        try:
            recipe = _resolve_recipe(store, args)
        except ConfigProblem as error:
            print(f"error: {error}", file=sys.stderr)
            return 2
        view = MsaView(cache=store, recipe=recipe)
        print(f"alignments {recipe.recipe_id}  ({recipe.describe()})")

    engines = [name.strip() for name in args.engines.split(",") if name.strip()]
    for name in engines:
        get_engine(name)

    outdir = Path(args.out)

    def paired_depth(job) -> int:
        if view is None:
            return 0
        found = [
            entry.depth_paired
            for entry in (view.cache.read_meta(view.recipe, p.key) for p in job.proteins)
            if entry is not None
        ]
        return min(found) if found else 0

    manifest_path = outdir / "manifest.tsv"
    rows = []
    if args.resume and manifest_path.is_file():
        # Keep what an earlier run recorded, so its rows survive into the
        # rewritten manifest rather than being replaced by this run's subset.
        rows = [Row(**_coerce(entry)) for entry in read_manifest(manifest_path)]

    without_msa = 0
    failures: list[tuple[str, str]] = []
    done: list[str] = []
    skipped_total = 0
    all_written: list = []
    for name in engines:
        todo = jobs
        if args.resume:
            already = completed_jobs(manifest_path, name)
            todo = [job for job in jobs if job.name not in already]
            skipped_total += len(jobs) - len(todo)
            if not todo:
                done.append(name)
                continue
        # An engine that cannot be written — a missing optional dependency, a
        # ligand its format refuses — must not cost the panel the other seven.
        try:
            if name in BATCH_WRITERS and not args.no_batch:
                # These engines take a whole panel in one file, which is also
                # how their daily upload quotas are meant to be filled.
                module = importlib.import_module(BATCH_WRITERS[name])
                produced = module.write_batch(todo, outdir / name, view)
            else:
                produced = [get_engine(name)(job, outdir / name, view) for job in todo]
            if name == "colabfold":
                colabfold.write_queries_csv(produced, outdir / name)
        except Exception as error:  # noqa: BLE001 - reported per engine, never fatal
            failures.append((name, str(error)))
            continue
        done.append(name)
        all_written.extend(produced)
        for written in produced:
            if view is not None and not written.msa_used:
                without_msa += 1
            rows.extend(rows_for(written, n_paired=paired_depth(written.job)))

    manifest = write_manifest(manifest_path, rows)
    write_panel(
        outdir / PANEL,
        all_written,
        seeds=tuple(int(s) for s in args.seeds.split(",")) if args.seeds else (),
        mode=args.mode,
        recipe=view.recipe.recipe_id if view else "",
        cache=str(args.cache) if args.cache else "",
    )
    print(f"complexes  {len(jobs)}")
    if args.resume:
        print(f"resumed    {skipped_total} already written, skipped")
    print(f"engines    {', '.join(done) if done else 'none'}")
    print(f"files      {len(jobs) * len(done)} under {outdir}")
    print(f"manifest   {manifest}")
    for name, reason in failures:
        print(f"\nskipped {name}: {reason}", file=sys.stderr)
    if view is not None and without_msa == len(jobs) * len(done) and done:
        print(
            "\nerror: the cache attached no alignment to any input.\n"
            "  Every engine would search for itself, once per complex, which is\n"
            "  what the cache exists to avoid — and some refuse to start.\n"
            f"  Check that {args.cache} holds this library's sequences.",
            file=sys.stderr,
        )
        return 2
    if without_msa:
        print(
            f"\nwarning: {without_msa} of {len(jobs) * len(engines)} inputs have no cached\n"
            "         alignment and will make the engine search for itself, once per\n"
            "         complex. Run the import or a search backend first.",
            file=sys.stderr,
        )
    return 1 if failures and not done else 0


def _replay(outdir: Path, engines: list[str]):
    """Rebuild what `write` produced, from the panel snapshot it left behind."""
    from foldrunner.engines.base import WrittenJob
    from foldrunner.ir import Component, Job, Ligand, Protein

    snapshot = read_panel(outdir / PANEL)
    seeds = tuple(snapshot.get("seeds", ()))
    out = []
    for item in snapshot["written"]:
        if engines and item["engine"] not in engines:
            continue
        components = tuple(
            Component(
                Protein(c["entity"], "A") if c["protein"] else Ligand(c["entity"], smiles="C"),
                c["count"],
            )
            for c in item["components"]
        )
        out.append(
            WrittenJob(
                job=Job(item["job"], components, seeds=seeds),
                engine=item["engine"],
                path=Path(rebase(item["path"], outdir)),
                chains=[list(x) for x in item["chains"]],
                msa_used=item["msa_used"],
                extra={k: rebase(v, outdir) for k, v in item.get("extra", {}).items()},
            )
        )
    return out


def cmd_env(args: argparse.Namespace) -> int:
    """Show which engines this machine is configured to run."""
    try:
        conf = config_module.load(args.config)
    except config_module.ConfigError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    if conf.path is None:
        print("no config file found; looked at:")
        for candidate in config_module.SEARCH_PATHS:
            print(f"  {candidate}")
        print(f"\nset {config_module.ENV_VAR} or pass --config to use one elsewhere.")
    else:
        print(f"config     {conf.path}")
    if conf.conda_sh:
        print(f"conda      {conf.conda_sh}")
    if conf.output_root:
        print(f"output     {conf.output_root}")

    print("\nengines:  [x] configured  [ ] on PATH  [-] disabled  [~] web upload")
    for name in sorted(ENGINES):
        entry = conf.for_engine(name)
        template = engine_command(name)
        if template is None:
            note = "web upload, no command"
        elif name in conf.engines:
            note = entry.describe()
        else:
            note = "not configured (will be run as found on PATH)"
        if template is None:
            flag = "~"          # nothing to run
        elif not entry.enabled:
            flag = "-"          # switched off in the config
        elif name in conf.engines:
            flag = "x"          # configured and ready
        else:
            flag = " "          # will be tried as found on PATH
        print(f"  [{flag}] {name:14s} {note}")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    """Run, or emit scripts for, a panel that has already been written."""
    try:
        conf = config_module.load(args.config)
    except config_module.ConfigError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    outdir = Path(args.panel)
    if not (outdir / PANEL).is_file():
        print(
            f"error: no {PANEL} in {outdir}; run `foldrunner write` first",
            file=sys.stderr,
        )
        return 1

    engines = [n.strip() for n in args.engines.split(",") if n.strip()] if args.engines else []
    written = _replay(outdir, engines)
    out_root = Path(args.out or conf.output_root or outdir / "results")
    invocations = plan(written, conf, out_root, engines or None)

    runnable = [i for i in invocations if not i.manual]
    manual = [i for i in invocations if i.manual]
    print(f"panel      {outdir}")
    print(f"results    {out_root}")
    print(f"runnable   {len(runnable)} invocations")
    if manual:
        print(f"by hand    {len(manual)} ({', '.join(sorted({i.engine for i in manual}))})")

    if args.runner == "script":
        paths = write_scripts(invocations, outdir / "scripts")
        print(f"\nwrote {len(paths)} files into {outdir / 'scripts'}")
        print(f"run them with: {outdir / 'scripts' / 'run_all.sh'}")
        return 0

    if args.dry_run:
        for invocation in runnable:
            print(f"\n--- {invocation.engine} / {invocation.label} ---")
            print(invocation.script)
        return 0

    log_dir = outdir / "logs"
    failed = 0
    for invocation in runnable:
        print(f"\n>>> {invocation.engine} / {invocation.label}")
        result = execute(invocation, timeout=args.timeout, log_dir=log_dir)
        if result.silent_failure:
            failed += 1
            print(
                f"    exited 0 but wrote no structure into {invocation.out}",
                file=sys.stderr,
            )
            if result.log:
                print(f"    full output: {result.log}", file=sys.stderr)
        elif result.ok:
            print(f"    done ({invocation.out}, {result.produced} structures)")
        else:
            failed += 1
            print(f"    failed rc={result.returncode}", file=sys.stderr)
            for line in result.tail().splitlines():
                print(f"    | {line}", file=sys.stderr)
            if result.log:
                print(f"    full output: {result.log}", file=sys.stderr)
            if not args.keep_going:
                return 1
    for invocation in manual:
        print(f"\nby hand: {invocation.reason}")
    return 1 if failed else 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="foldrunner",
        description=(
            "Enumerate all-vs-all protein pairs, reuse one alignment per sequence, "
            "and write native inputs for every supported prediction engine."
        ),
    )
    parser.add_argument("--version", action="version", version=f"foldrunner {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    def add_recipe_args(p: argparse.ArgumentParser) -> None:
        p.add_argument("--backend", default="colabfold-api", help="alignment source name")
        p.add_argument("--databases", default="", help="comma-separated database names")
        p.add_argument("--tool-version", default="", help="version of the search tool")

    def add_panel_args(p: argparse.ArgumentParser) -> None:
        p.add_argument("fasta", type=Path, help="protein library in FASTA format")
        p.add_argument(
            "--mode", default="all", help="all, hetero, monomer or bipartite (default: all)"
        )
        p.add_argument("--preys", type=Path, help="second library, for bipartite mode")

    plan = sub.add_parser("plan", help="report job and alignment counts without writing")
    add_panel_args(plan)
    add_recipe_args(plan)
    plan.add_argument("--engines", default="boltz2,chai1")
    plan.add_argument("-v", "--verbose", action="store_true", help="list every complex")
    plan.set_defaults(func=cmd_plan)

    imp = sub.add_parser("import", help="adopt alignments from an earlier run")
    imp.add_argument("source", type=Path, help="directory to scan")
    imp.add_argument("--cache", type=Path, required=True, help="alignment cache directory")
    imp.add_argument("--layout", default="auto", help="auto, protenix or colabfold")
    imp.add_argument("--link", action="store_true", help="hard link instead of copying")
    imp.add_argument("--dry-run", action="store_true", help="scan without writing")
    imp.add_argument("--limit", type=int, help="stop after this many entries")
    imp.add_argument("--show-failures", action="store_true")
    add_recipe_args(imp)
    imp.set_defaults(func=cmd_import)

    search = sub.add_parser("search", help="compute the alignments a panel needs")
    add_panel_args(search)
    search.add_argument("--cache", type=Path, required=True, help="alignment cache directory")
    search.add_argument("--contact", default="", help="email the services can reach you at")
    search.add_argument(
        "--local-db", type=Path, help="search locally against this database directory"
    )
    search.add_argument("--threads", type=int, default=0, help="threads for a local search")
    search.add_argument("--database-version", default="", help="local database release")
    search.add_argument("--tool-version", default="", help="version of the search tool")
    search.add_argument("--host", default=PROTENIX_HOST, help="primary alignment service")
    search.add_argument("--fallback", default=COLABFOLD_HOST, help="service to try if throttled")
    search.add_argument("--batch-size", type=int, default=20, help="sequences per request")
    search.add_argument("--pairing-strategy", default="greedy", help="greedy or complete")
    search.add_argument("--unpaired-only", action="store_true", help="skip the pairing search")
    search.add_argument("--no-env", action="store_true", help="search UniRef only")
    search.set_defaults(func=cmd_search)

    write = sub.add_parser("write", help="generate engine inputs")
    add_panel_args(write)
    add_recipe_args(write)
    write.add_argument("-o", "--out", type=Path, required=True, help="output directory")
    write.add_argument("--cache", type=Path, help="alignment cache directory")
    write.add_argument("--recipe", help="which alignments in the cache to use, by id")
    write.add_argument(
        "--engines",
        default=",".join(sorted(ENGINES)),
        help=f"comma-separated engines (available: {', '.join(sorted(ENGINES))})",
    )
    write.add_argument("--seeds", default="", help="comma-separated integer seeds")
    write.add_argument(
        "--no-batch", action="store_true", help="one file per complex even for batch engines"
    )
    write.add_argument(
        "--resume", action="store_true", help="skip complexes the manifest already records"
    )
    write.set_defaults(func=cmd_write)

    env = sub.add_parser("env", help="show which engines this machine can run")
    env.add_argument("--config", type=Path, help="config file to read")
    env.set_defaults(func=cmd_env)

    run = sub.add_parser("run", help="run a written panel, or emit scripts for it")
    run.add_argument("panel", type=Path, help="directory produced by `foldrunner write`")
    run.add_argument("--config", type=Path, help="config file to read")
    run.add_argument("--engines", default="", help="restrict to these engines")
    run.add_argument("-o", "--out", type=Path, help="where results go")
    run.add_argument(
        "--runner",
        default="script",
        choices=RUNNERS,
        help="script: write shell scripts (default); local: run them here",
    )
    run.add_argument("--dry-run", action="store_true", help="print commands instead of running")
    run.add_argument("--keep-going", action="store_true", help="continue after a failure")
    run.add_argument("--timeout", type=float, help="seconds per invocation")
    run.set_defaults(func=cmd_run)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
