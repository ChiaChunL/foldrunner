# foldrunner: all-vs-all prediction panels with one alignment per sequence

<img src="https://raw.githubusercontent.com/ChiaChunL/foldrunner/main/docs/assets/foldrunner_banner.png" alt="foldrunner" width="100%">

| Testing | [![CI](https://github.com/ChiaChunL/foldrunner/actions/workflows/ci.yml/badge.svg)](https://github.com/ChiaChunL/foldrunner/actions/workflows/ci.yml) |
|---|---|
| Package | [![PyPI Latest Release](https://img.shields.io/pypi/v/foldrunner.svg)](https://pypi.org/project/foldrunner/) [![Python versions](https://img.shields.io/pypi/pyversions/foldrunner.svg)](https://pypi.org/project/foldrunner/) [![PyPI total downloads](https://img.shields.io/pepy/dt/foldrunner.svg?label=total%20downloads)](https://pepy.tech/projects/foldrunner) |
| Meta | [![License - BSD 3-Clause](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](LICENSE) [![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff) |

Download statistics: [PyPI totals (Pepy)](https://pepy.tech/projects/foldrunner).

## 🧬 What is it?

`foldrunner` takes a protein library and produces the input files for a complete
all-vs-all prediction panel across **AlphaFold 3, AlphaFold Server, Boltz-2,
Protenix, Chai-1, ColabFold, AlphaFold2-Multimer and SeedFold** — computing each
multiple sequence alignment once and handing the same one to every engine.

A panel of *N* sequences contains *N(N+1)/2* complexes but only *N* distinct
sequences. The engines cannot see that: left to itself AlphaFold 3 runs its data
pipeline once per job, so a 100-sequence panel searches the same sequence up to
100 times. That difference is the point of the package.

```
library.fasta ──► enumerate ──► N(N+1)/2 complexes
                      │
                      └──► N alignments  (searched once, addressed by sequence hash)
                                │
      ┌───────────┬─────────────┼──────────────┬────────────┐
  AlphaFold3   Boltz-2      Protenix        Chai-1     ColabFold  …
      └───────────┴─────────────┼──────────────┴────────────┘
                          manifest.tsv ──► foldmetrics
```

Scores from the resulting predictions are read by
[foldmetrics](https://github.com/ChiaChunL/foldmetrics), which computes ipSAE,
pDockQ2, LIS and DockQ over the same panel. It reads each engine's output
tree untouched, so no file has to be renamed or moved first:

```bash
pip install "foldrunner[metrics]"     # pulls in foldmetrics >= 0.1.7
fmx score results/ --by-target summary.tsv
```

Checked against the untouched output of all six engines this package can run:
foldmetrics 0.1.7 resolves each one to the same six targets, and counts the
models each engine actually produced rather than the files it wrote.

| engine | targets | models | what could go wrong |
|---|---|---|---|
| AlphaFold 3 | 6 | 60 | the top-level copy of the best sample counted twice |
| Protenix | 6 | 60 | `seed_N/predictions/` swallowing the job name |
| Chai-1 | 6 | 30 | |
| AlphaFold2-Multimer | 6 | 30 | `.pdb` and `.cif` of one model counted as two |
| Boltz-2 | 6 | 12 | two seed directories read as two targets |
| ColabFold | 6 | 6 | a flat directory collapsing the whole panel into one |

Older releases do collapse a panel into a single row on some of these layouts,
so the floor is not cosmetic.

## 📦 Installation

```bash
pip install foldrunner
```

With Chai-1 alignment support (writes `.aligned.pqt` parquet files):

```bash
pip install "foldrunner[chai]"
```

Development install:

```bash
git clone https://github.com/ChiaChunL/foldrunner.git
cd foldrunner
pip install -e ".[dev]"
```

## ⚡ Quickstart

CLI (`foldrunner`, short alias `frx`):

```bash
# what would this panel cost, before anything is submitted?
frx plan library.fasta --mode all

# compute the alignments: once per unique sequence, not once per pair
frx search library.fasta --cache msa/ --contact you@example.org

# ... or adopt alignments an earlier run already produced
frx import old_msa_tree/ --cache msa/ --backend protenix-api --link

# write the inputs
frx write library.fasta -o panel/ --cache msa/ --engines af3,boltz2,protenix,chai1
```

`plan` reports the arithmetic that makes the rest worthwhile:

```
library        87 unique sequences
mode           all
complexes      3828
MSA searches   87   (one per unique sequence, not per complex)
               searching per complex instead would cost 3828
engine files   15312
```

Python — the CLI is a thin layer over this:

```python
import foldrunner as fr

jobs = fr.panel("library.fasta", "all", seeds=(2066,))

fr.cost(jobs, engines=["af3", "boltz2"])
# PanelCost(complexes=3828, sequences=87, engines=2)
#   .searches_saved -> 3741      .files -> 7656

fr.search(jobs, cache="./msa", backend=backend)     # one search per sequence
written = fr.write(jobs, "./panel", engines=["af3", "boltz2"],
                   cache="./msa", recipe=backend.recipe)
fr.run_plan(written, "./results", scripts="./panel/scripts")
```

## 🛠️ Command-line reference

| Command | Purpose |
|---|---|
| `plan` | report complex and alignment counts without writing anything |
| `search` | compute the alignments a panel needs, web service or local databases |
| `import` | re-address alignments from an earlier run into the cache |
| `write` | generate the input files, plus `manifest.tsv` |

| Option | Commands | Meaning |
|---|---|---|
| `--mode {all,hetero,monomer,bipartite}` | plan, search, write | which complexes to enumerate (default `all`) |
| `--preys FILE` | plan, search, write | second library, for `bipartite` |
| `--cache DIR` | search, import, write | alignment cache directory |
| `--engines LIST` | plan, write | comma-separated engines; default is all eight |
| `--seeds LIST` | write | comma-separated integer seeds |
| `-o, --out DIR` | write | output directory |
| `--no-batch` | write | one file per complex even for the batch engines |
| `--contact EMAIL` | search | the web services require callers to identify themselves |
| `--host` / `--fallback` | search | primary and standby alignment services |
| `--batch-size N` | search | sequences per request (default 20) |
| `--pairing-strategy {greedy,complete}` | search | how strictly taxa must be shared |
| `--unpaired-only` | search | skip the pairing search |
| `--link` | import | hard link instead of copying, for large trees |
| `--dry-run` / `--limit N` | import | scan without writing / stop early |
| `--backend` / `--databases` / `--tool-version` | import, write | which recipe the alignments belong to |

`frx <command> --help` prints the complete option list.

## 🔗 Supported engines

| Engine | Input | Copies expressed as | Keeps chain names | External alignment |
|---|---|---|---|---|
| AlphaFold 3 | JSON per complex | `id` list | yes | `unpairedMsaPath` / `pairedMsaPath` |
| AlphaFold Server | JSON batches | `count` | no (forces A/B) | `unpairedMsa`, embedded |
| Boltz-2 | YAML per complex | `id` list | yes | `msa:` path |
| Protenix | one JSON panel | `count` | no (forces A/B) | `unpairedMsaPath` / `pairedMsaPath` |
| Chai-1 | FASTA per complex | repeated records | — | `<sha256>.aligned.pqt` directory |
| ColabFold | FASTA or CSV | a3m header cardinality | — | concatenated a3m, or CSV `a3mpath` |
| AlphaFold2-Multimer | FASTA per complex | repeated records | no (starts at B) | `--use-precomputed-msas`, partial at best |
| SeedFold | JSON, single or batch | `copies` | — | not supported |

Every format was checked against upstream source rather than documentation prose;
[docs/engine-formats.md](docs/engine-formats.md) records each one with the
inconsistencies that are real and must not be "tidied up" — AlphaFold Server's
separate dialect and nineteen-code ligand allow-list, ColabFold joining chains
with `:` where AlphaFold2-Multimer uses separate records, and Boltz-2 refusing a
directory that mixes entries with and without an `affinity` property.

Every format was additionally checked by feeding the generated files to each
engine's *own* parser — `folding_input.Input.from_json`, `boltz.main.check_inputs`,
`SampleDictToFeatures`, `read_inputs`, `get_queries`, `parse_fasta` — rather than
only to this package's tests.

**Two engines cannot fully reuse an alignment, which any comparison has to say
out loud.** SeedFold accepts neither a precomputed alignment nor a seed, so it
always searches for itself and its sampling is not reproducible.
AlphaFold2-Multimer can at most reuse `bfd_uniref_hits.a3m`, the one file in its
precomputed set that is already a3m; the three Stockholm files it also looks for,
including the `uniprot_hits.sto` that drives pairing, are left for it to search,
because converting into them risks changing the alignment rather than reusing it.
The manifest records this as `msa_partial`. Whether even that much is reused
depends on the entry point: `run_alphafold.py` decides per file, but a site
wrapper may refuse to start unless the whole set is present, in which case
AlphaFold2 searches everything itself.

## 🧮 Enumeration

| Mode | Complexes | Use |
|---|---|---|
| `all` | *N(N+1)/2* | every unordered pair, self-pairs included |
| `hetero` | *N(N-1)/2* | pairs only |
| `monomer` | *N* | one job per sequence, as controls |
| `bipartite` | *baits × preys* | screening candidates against known partners |

A **self-pair is a homodimer**: one entity with `count=2`, not two components
that happen to share a sequence. Six of the eight engines express that directly;
Chai-1 and AlphaFold2-Multimer expand it into repeated records.

Sequences are deduplicated by content on load, so the same protein under two
labels is searched once and produces one job rather than two.

## 🧠 Alignments

Cache entries are keyed by `sha256(sequence.upper())` — the same address Chai-1
uses for `<hash>.aligned.pqt`, so a cache directory doubles as its
`msa_directory` with no conversion step. Addressing by content rather than by
position is what lets a cache survive a reordered library, a renamed directory,
or reuse by the next panel entirely.

```
msa/
└── protenix-api-6b0f3f937902/          ← recipe: backend + databases + version
    └── 3dd8c3f8…fcfaaa/                ← sha256 of the sequence
        ├── unpaired.a3m
        ├── paired.a3m
        └── meta.json
```

The **recipe** layer is not decoration. Alignments from different backends or
database releases have different depths, depth moves ipTM, and a panel that
silently mixes two sources produces differences that cannot be attributed
afterwards.

### Paired and unpaired

The two are fetched and stored separately because they come from different
places. Pairing needs a taxonomy-labelled source — hits arrive as
`UniRef100_<accession>_<taxid>/` — while the metagenomic databases that give
unpaired alignments their depth carry no taxon at all.

Pairing follows AlphaFold2-Multimer: group by taxon, drop taxa present in only
one chain, rank within a taxon by similarity to that chain's own query.

**This matters for any panel spanning more than one species.** Two human proteins
share orthologues in mouse, zebrafish and yeast, so their paired alignment runs
to hundreds of rows. A human protein and a viral one share no organism at all, so
theirs is nearly empty and the model has only its priors to work with. Their
interface scores are therefore not on the same footing, and sorting a mixed panel
by ipTM ranks pairing depth as much as it ranks interaction. `manifest.tsv`
carries `n_paired` per job so the two can be read apart —
[docs/msa-pairing.md](docs/msa-pairing.md) works through it.

### Where alignments come from

| Backend | When |
|---|---|
| Protenix / ColabFold services | up to a few dozen unique sequences; either can stand in for the other |
| a local MMseqs2 search you run yourself | anything larger — the public services are shared academic resources sized for a few thousand alignments a day across all users. Run `colabfold_search` and adopt the result with `foldrunner import`; `--local-db` is groundwork and refuses to run |
| `frx import` | alignments an earlier run already produced |

`plan` and `search` warn before a panel crosses into territory the public
services are not meant to absorb.

## 🚀 Running

`foldrunner write` produces inputs; `foldrunner run` turns them into commands.
The engine modules know the shape of each invocation, and a config file says
where the engines are installed — so moving a panel to another machine changes
the config, not the package.

```toml
# ./foldrunner.toml, ~/.config/foldrunner/foldrunner.toml, or $FOLDRUNNER_CONFIG

[defaults]
# A non-interactive shell is not a login shell: conda is not on PATH, and
# `conda activate` is not even defined until this is sourced.
conda_sh = "/opt/conda/etc/profile.d/conda.sh"

[engines.boltz2]
conda_env = "boltz2"
extra_args = "--cache /opt/boltz-2"

[engines.protenix]
conda_env = "protenix"
env = { PROTENIX_ROOT_DIR = "/opt/Protenix_data" }

[engines.chai1]
conda_env = "chai1"
# Chai-1 fetches its ESM-2 embeddings from Hugging Face even though its own
# weights come from elsewhere; this is where that download is redirected.
env = { HF_ENDPOINT = "https://hf-mirror.com" }

[engines.af3]
container = "docker"
image = "alphafold3"
gpus = "all"
# Alignment paths inside an AlphaFold 3 job are absolute, so the input has to be
# visible at the same path inside the container as outside it.
mounts = ["{input_dir}:{input_dir}", "/opt/af3_models:/models"]
extra_args = "--model_dir=/models --norun_data_pipeline"

[engines.af2_multimer]
executable = "/opt/af2/run_af2.sh"
# The multimer preset is not the default anywhere; the monomer pipeline rejects
# the multi-record FASTA this engine writes, after the job has started.
extra_args = "--model-preset multimer"
```

```bash
foldrunner env                 # what can this machine run?
foldrunner run panel -o results               # writes panel/scripts/run_all.sh
foldrunner run panel -o results --runner local
```

A worked configuration is in [examples/foldrunner.toml](examples/foldrunner.toml).

How much one invocation covers differs by engine and is not a detail: Boltz-2 and
ColabFold read a whole directory, so they run once for the panel rather than once
per complex; Protenix takes the panel and every seed in one file; AlphaFold 3,
Chai-1 and AlphaFold2-Multimer run per complex.

Chai-1 is driven through a generated Python script rather than `chai-lab fold`,
because the command line writes only the ranking summary and the full PAE exists
solely on the object the Python API returns — a CLI run loses it permanently.

AlphaFold Server and SeedFold have no command at all. They are web forms, and are
listed in `scripts/MANUAL.md` with the file to upload rather than given an
invented command line.

Failures keep their whole output under `panel/logs/`, script included. Re-running
a GPU job to read a message it already printed is an expensive way to learn
nothing new.

## 🧪 Examples

[examples/](examples/) ships a three-sequence library — barnase, barstar and
ubiquitin — so every command above runs as-is after a clone. It holds one known
binder, two known non-binders and three homodimers, which is enough to tell a
working run from a broken one before committing a real library to it.

## 📋 Manifest

`write` produces `manifest.tsv` alongside the inputs:

```
job         engine    component  entity   chain  count  n_chains  n_residues  msa_used  n_paired
A__A        boltz2    0          A        A      2      2         1334        1         1681
A__A        boltz2    0          A        B      2      2         1334        1         1681
A__B        protenix  0          A        A      1      2         1090        1         1444
```

Three engines discard the chain names they are given — Protenix and AlphaFold
Server relabel everything to A/B/C, AlphaFold2-Multimer starts at B — so results
coming back from different engines cannot be lined up against the same entities
without this table. It also carries the paired depth per job, and lets a run
resume where it stopped.

## 📁 Outputs and paths

- `write -o DIR` → one subdirectory per engine, plus `manifest.tsv`.
- Batch engines write differently by design: Protenix gets one `panel.json`,
  AlphaFold Server and SeedFold get `batch_000.json`, `batch_001.json`, … sized
  to a day's upload quota. `--no-batch` forces one file per complex.
- ColabFold additionally gets `queries.csv` with an `a3mpath` column, which is
  how it takes precomputed alignments.
- Boltz-2 jobs carrying an affinity property go to `boltz2/affinity/`. Boltz-2
  runs the affinity stage over every entry in a directory, so an entry without
  the property fails there; keeping them apart is the only way both run.
- Job names are `<A>__<B>`, order-independent and sanitized to `[A-Za-z0-9._-]`.
  The same name is used by every engine so results can be matched back.

## 💡 Conventions worth knowing

- **Entities, not chains.** The intermediate representation models an entity and
  a copy count. Chain identifiers are assigned by each writer and recorded, never
  carried in the job, because for three engines any name chosen up front is
  fiction.
- **Never let an engine build its own alignment.** It is not only slower; the
  alignment then differs between engines, and the comparison stops being one.
- **Alignments are referenced by path, not embedded**, everywhere except
  AlphaFold Server, which has no filesystem to point at. A protein's alignment
  runs to tens of megabytes and a panel repeats it across every complex the
  sequence appears in.
- **Throttling looks like a corrupt download.** The services answer a throttled
  request with an HTML page under a normal status code; unpacking it raises a
  gzip error far from its cause. foldrunner checks the archive magic first and
  reports it as throttling, which is worth waiting out.
- **Ligands add PAE tokens.** One per heavy atom, so a residue-count-based
  alignment of a PAE matrix will be wrong for any complex containing one.

## ✅ What has actually been run

Formats were checked against each engine's own parser, but a parser accepting an
input does not mean the engine can run it — every defect this package has had
passed format validation and failed after the job started. So this table
separates what has produced a structure from what has only been exercised
against a stand-in.

| | status |
|---|---|
| AlphaFold 3, Protenix, Boltz-2, Chai-1, ColabFold, AlphaFold2-Multimer | run end to end on GPU; 6-complex panel, 2 seeds, structures scored by foldmetrics |
| AlphaFold Server, SeedFold | inputs match the published schema and official examples; **never submitted** — they are web forms with no API |
| Alignments from the web services | real calls against both endpoints |
| `foldrunner import` | 22,804 real entries, hard-linked, no duplicated bytes |
| `--local-db` (built-in local search) | **not available** — written but never run against real databases, so it refuses rather than pretending. Run `colabfold_search` yourself and use `foldrunner import` |

To run a panel on another machine, install foldrunner there and copy the panel
across — the snapshot records relative paths, so it works from wherever it
lands. There is no remote-execution mode, because driving an engine over ssh
needs the panel visible at the same path on both machines and that assumption
does not hold often enough to be worth the failure modes.

## 📚 References

- Abramson J et al. [*Accurate structure prediction of biomolecular interactions
  with AlphaFold 3.*](https://doi.org/10.1038/s41586-024-07487-w) Nature 630,
  493–500 (2024).
- Evans R et al. [*Protein complex prediction with
  AlphaFold-Multimer.*](https://doi.org/10.1101/2021.10.04.463034) bioRxiv (2021).
  — the species-pairing scheme this package follows
- Mirdita M et al. [*ColabFold: making protein folding accessible to
  all.*](https://doi.org/10.1038/s41592-022-01488-1) Nat Methods 19, 679–682
  (2022). — the MMseqs2 alignment service and its a3m conventions
- Steinegger M, Söding J. [*MMseqs2 enables sensitive protein sequence searching
  for the analysis of massive data sets.*](https://doi.org/10.1038/nbt.3988)
  Nat Biotechnol 35, 1026–1028 (2017).
- Wohlwend J et al. [*Boltz-1: democratizing biomolecular interaction
  modeling.*](https://doi.org/10.1101/2024.11.19.624167) bioRxiv (2024).
- Chai Discovery. [*Chai-1: decoding the molecular interactions of
  life.*](https://doi.org/10.1101/2024.10.10.615955) bioRxiv (2024).
- Chen X et al. [*Protenix: advancing structure prediction through a
  comprehensive AlphaFold3 reproduction.*](https://doi.org/10.1101/2025.01.08.631967)
  bioRxiv (2025).
- ByteDance Seed. [*SeedFold: scaling biomolecular structure
  prediction.*](https://arxiv.org/abs/2512.24354) arXiv (2025).

## 📄 License

[BSD 3-Clause](LICENSE)
