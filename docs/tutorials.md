# Tutorials

## Tutorial philosophy

The tutorials below are written as practical analysis pathways rather than isolated command snippets. The goal is to mirror how Ariadne is typically used in a research workflow: run the pipeline, inspect the classification layer, then move into phylogenetic interpretation.

## Tutorial 1. Reproduce the bundled example

The repository already contains:

- `input/` as an example protein input folder
- `tree/` as the prepared reference directory

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

This is the best starting point if you want to see the full Ariadne workflow on a known example before adapting it to your own dataset.

## Tutorial 2. Reconstruct the workflow stage by stage

Running stage by stage is useful when you want to inspect intermediate files or alter one stage without rerunning the rest.

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

This stage-wise decomposition is particularly helpful when you are tuning filtering thresholds or comparing different HMM sources.

## Tutorial 3. Read the classification outputs as a screening layer

The classification stage is often the most informative point for early interpretation.

Start with:

- `classification.tsv`
- `nearest_neighbors.tsv`
- `embedding.svg`
- `embedding_3d_sections.svg`
- `ceess_predictions.tsv`
- `ceess_candidates.tsv`
- `ceess_candidates.fasta`

Suggested reading order:

1. open `classification.tsv` to inspect the predicted source assignment for each candidate
2. read `nearest_neighbors.tsv` to understand which references support that assignment
3. inspect `embedding.svg` for a global two-dimensional view
4. inspect `embedding_3d_sections.svg` for a more presentation-ready geometric summary
5. read `ceess_predictions.tsv` to see which coral-like candidates were scored toward `cembrene A / cembrene B`
6. open `ceess_candidates.tsv` and `ceess_candidates.fasta` for the final CeeSs shortlist

In practice, this is the stage where you decide which candidates deserve deeper phylogenetic attention.

## Tutorial 4. Run the standalone coral TPS ESM model

If you want to inspect the supervised coral TPS type model independently of de novo candidate discovery:

```bash
.venv/bin/python -m ariadne esm-type \
  --xlsx TPS/TPS.xlsx \
  --output-dir tmp_esm_results
```

Recommended files to inspect:

- `esm_embedding.svg`
- `esm_metrics.tsv`
- `esm_confusion_matrix.tsv`
- `esm_predictions.tsv`

This standalone mode is useful when you want to evaluate how well the labeled coral TPS proteins in `TPS/TPS.xlsx` separate in ESM space before using that model inside the main classification workflow.

## Tutorial 5. Read the phylogeny outputs as an evolutionary layer

The phylogeny stage provides the alignment-driven context that follows the feature-space screening step.

Focus on:

- `phylogeny_alignment.fasta`
- `phylogeny_sequence_map.tsv`
- `iqtree.treefile`
- `iqtree.iqtree`

The sequence map is especially important because it connects the tree-safe identifiers back to the original FASTA headers.

This stage is the most useful when you want to compare screened candidates against known reference clades, examine local placement, or prepare figures for a manuscript or presentation.

## Tutorial 6. Use your own prebuilt models

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

These modes are useful when you want to preserve Ariadne's downstream classification and phylogeny logic while substituting your own HMM resources upstream.

## Tutorial 7. Iterate quickly without rebuilding the final tree

When your immediate goal is threshold tuning or candidate triage, the fastest practical loop is:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --skip-phylogeny \
  --output-dir tmp_fast_iteration/
```

This allows repeated inspection of:

- `classification.tsv`
- `nearest_neighbors.tsv`
- `embedding.svg`

without rerunning MAFFT and IQ-TREE at every iteration.

## Suggested workflow for manuscript preparation

If you are preparing figures or supplement-style analyses, a good working order is:

1. generate a complete run with `ariadne run`
2. screen candidates using `classification.tsv` and `embedding.svg`
3. interpret the final placement with `iqtree.treefile`
4. archive `pipeline_summary.tsv`, `classification.tsv`, and the phylogeny outputs together for reproducibility
