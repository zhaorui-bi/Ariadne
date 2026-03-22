# Tutorials

## Tutorial 1. Reproduce the bundled example

The repository already contains an example protein input folder and a prepared reference folder:

- `input/`
- `tree/`

Run:

```bash
.venv/bin/python -m ariadne run \
  --protein-folder input \
  --reference-dir tree \
  --output-dir tmp_run_example
```

Expected stage directories:

- `tmp_run_example/01_discovery/`
- `tmp_run_example/02_filtering/`
- `tmp_run_example/03_classification/`
- `tmp_run_example/04_phylogeny/`

## Tutorial 2. Run stage by stage

### Discovery

```bash
.venv/bin/python -m ariadne discover \
  --protein-folder input \
  --hmm tmp_run_example/01_discovery/query.hmm \
  --output-dir tmp_manual/01_discovery
```

### Filtering

```bash
.venv/bin/python -m ariadne filter \
  --input-fasta tmp_manual/01_discovery/candidates.protein.faa \
  --output-dir tmp_manual/02_filtering
```

### Classification

```bash
.venv/bin/python -m ariadne classify \
  --candidates tmp_manual/02_filtering/candidates.filtered.faa \
  --reference-dir tree \
  --output-dir tmp_manual/03_classification
```

### Phylogeny

```bash
.venv/bin/python -m ariadne phylogeny \
  --candidates tmp_manual/02_filtering/candidates.filtered.faa \
  --reference-dir tree \
  --output-dir tmp_manual/04_phylogeny
```

## Tutorial 3. Interpret the classification stage

The most informative files are:

- `classification.tsv`
- `nearest_neighbors.tsv`
- `embedding.svg`
- `embedding_3d_sections.svg`

Suggested reading order:

1. open `classification.tsv` to inspect the predicted reference source for each candidate
2. read `nearest_neighbors.tsv` to see which references support that assignment
3. inspect `embedding.svg` for a global 2D view
4. inspect `embedding_3d_sections.svg` for publication-style projections

## Tutorial 4. Interpret the phylogeny stage

Focus on:

- `phylogeny_alignment.fasta`
- `phylogeny_sequence_map.tsv`
- `iqtree.treefile`
- `iqtree.iqtree`

The sequence map is especially useful because it connects the tree-safe identifiers back to the original headers.

## Tutorial 5. Use your own prebuilt models

If you already have a custom discovery HMM:

```bash
ariadne run \
  --protein-folder my_inputs/ \
  --reference-dir tree/ \
  --query-hmm my_query.hmm \
  --output-dir my_results/
```

If you already have a custom TPS HMM library:

```bash
ariadne run \
  --protein-folder my_inputs/ \
  --reference-dir tree/ \
  --tps-hmm-dir my_tps_hmms/ \
  --output-dir my_results/
```

## Tutorial 6. Skip phylogeny during fast iteration

If you are tuning thresholds and only want discovery, filtering, and classification:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --skip-phylogeny \
  --output-dir tmp_fast_iteration/
```

This is useful when you need to inspect `classification.tsv` and `embedding.svg` repeatedly without rerunning MAFFT and IQ-TREE.
