# Examples

`data/library.fasta` holds three well-characterised sequences: barnase and
barstar, which form one of the best-studied protein complexes, and ubiquitin,
which binds neither. An all-vs-all panel over them therefore contains a known
positive, two known negatives, and three homodimers — enough to tell whether a
run is working before spending anything on a real library.

Every command below runs as-is from a clone.

```bash
# 6 complexes from 3 sequences, and the alignment count that follows from it
foldrunner plan data/library.fasta

# write inputs for two engines, no alignments yet
foldrunner write data/library.fasta -o panel --engines boltz2,chai1

# the same with alignments: one search per sequence, not per pair
foldrunner search data/library.fasta --cache msa --contact you@example.org
foldrunner write data/library.fasta -o panel --cache msa \
  --backend protenix-api --databases uniref30,colabfold_envdb \
  --tool-version mmseqs2-service --engines boltz2,chai1

# what would run here, and the scripts that would run it
foldrunner env
foldrunner run panel -o results
```

`foldrunner.toml` is a worked configuration for a machine with conda
environments, a container, and a second host. Copy it next to your panel, or to
`~/.config/foldrunner/foldrunner.toml`, and edit the paths.

Expect barnase–barstar to score far above the pairs involving ubiquitin. If it
does not, the alignments are the first thing to check: `manifest.tsv` reports
`n_paired` per complex, and a near-zero value there explains a low interface
score without the model being wrong about anything.
