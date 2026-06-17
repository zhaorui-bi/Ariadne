# Tutorials

These tutorials are written as practical analysis pathways rather than isolated snippets. They mirror how Ariadne is used in practice: run the pipeline, read the classification layer, then move into phylogenetic interpretation.

## 1 · Reproduce the bundled example

The repository ships everything needed for a complete run:

- `input/` — an example protein input folder;
- `tree/` — the prepared multi-clade reference directory.

```bash
.venv/bin/python -m ariadne run \
  --protein-folder input \
  --reference-dir tree \
  --output-dir tmp_run_example
```

Expected stage directories:

```text
tmp_run_example/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
└── 04_phylogeny/
```

This is the best starting point for seeing the full workflow on a known example before adapting it to your own data.

## 2 · Reconstruct the workflow stage by stage

Running stage by stage lets you inspect intermediate files or change one stage without rerunning the rest.

```bash
# Discovery
.venv/bin/python -m ariadne discover \
  --protein-folder input \
  --hmm tmp_run_example/01_discovery/query.hmm \
  --output-dir tmp_manual/01_discovery

# Filtering
.venv/bin/python -m ariadne filter \
  --input-fasta tmp_manual/01_discovery/candidates.protein.faa \
  --output-dir tmp_manual/02_filtering

# Classification
.venv/bin/python -m ariadne classify \
  --candidates tmp_manual/02_filtering/candidates.filtered.faa \
  --reference-dir tree \
  --output-dir tmp_manual/03_classification

# Phylogeny
.venv/bin/python -m ariadne phylogeny \
  --candidates tmp_manual/02_filtering/candidates.filtered.faa \
  --reference-dir tree \
  --output-dir tmp_manual/04_phylogeny
```

This decomposition is particularly useful when tuning filtering thresholds or comparing different HMM sources.

## 3 · Read the classification outputs as a screening layer

Classification is usually the most informative point for early interpretation. Start with:

- `classification.tsv` — predicted source assignment for each candidate;
- `nearest_neighbors.tsv` — the references that support each assignment;
- `embedding.svg` — a global two-dimensional view;
- `embedding_3d_sections.svg` — a presentation-ready geometric summary;
- `ceess_predictions.tsv`, `ceess_candidates.tsv`, `ceess_candidates.fasta` — the CeeSs shortlist (when ESM scoring is enabled).

A productive reading order:

1. open `classification.tsv` to inspect the predicted source for each candidate;
2. read `nearest_neighbors.tsv` to see which references support it;
3. inspect `embedding.svg`, then `embedding_3d_sections.svg`, for the geometry;
4. read `ceess_predictions.tsv` to see which coral-like candidates score toward cembrene A / B;
5. open `ceess_candidates.tsv` / `.fasta` for the final shortlist.

This is the stage at which you decide which candidates deserve deeper phylogenetic attention.

## 4 · Run the standalone coral TPS ESM model

To inspect the supervised coral TPS type model on its own — independently of *de novo* discovery — call the Python API directly. (There is no `esm-type` subcommand; the standalone analysis lives in `ariadne.model`.)

```python
import ariadne as ad

outputs = ad.analyze_tps_types_with_esm(
    "TPS/TPS.xlsx",
    output_dir="tmp_esm_results",
    classifier_kind="mlp",   # or "logreg"
)
print(outputs)               # dict of written file paths
```

Files worth inspecting:

- `esm_embedding.svg` — labeled projection of the ESM embedding space;
- `esm_metrics.tsv` — cross-validated accuracy and macro-F1;
- `esm_confusion_matrix.tsv` — multi-class confusion matrix;
- `esm_predictions.tsv` — per-record cross-validated predictions.

This is useful for evaluating how cleanly the labeled coral TPS proteins in `TPS/TPS.xlsx` separate in ESM space before relying on the model inside the main classification workflow.

## 5 · Read the phylogeny outputs as an evolutionary layer

The phylogeny stage provides the alignment-driven context that follows feature-space screening. Focus on:

- `phylogeny_alignment.fasta` — the MAFFT alignment;
- `phylogeny_sequence_map.tsv` — maps tree-safe identifiers back to original FASTA headers;
- `iqtree.treefile` — the final phylogeny (Newick);
- `iqtree.iqtree` — IQ-TREE's model and inference summary.

The sequence map matters because it connects tree-safe identifiers back to the original headers. This stage is most useful when comparing screened candidates against known reference clades or preparing figures for a manuscript.

## 6 · Use your own prebuilt HMMs

Substitute your own HMM resources upstream while keeping Ariadne's downstream classification and phylogeny logic:

```bash
# custom discovery HMM
ariadne run --protein-folder my_inputs/ --reference-dir tree/ \
  --query-hmm my_query.hmm --output-dir my_results/

# custom TPS HMM library
ariadne run --protein-folder my_inputs/ --reference-dir tree/ \
  --tps-hmm-dir my_tps_hmms/ --output-dir my_results/
```

## 7 · Iterate quickly without rebuilding the tree

When the goal is threshold tuning or candidate triage, skip the phylogeny stage for a fast loop:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --skip-phylogeny \
  --output-dir tmp_fast_iteration/
```

This allows repeated inspection of `classification.tsv`, `nearest_neighbors.tsv`, and `embedding.svg` without rerunning MAFFT and IQ-TREE on every iteration.

## Suggested workflow for manuscript preparation

1. generate a complete run with `ariadne run`;
2. screen candidates using `classification.tsv` and `embedding.svg`;
3. interpret the final placement with `iqtree.treefile`;
4. archive `pipeline_summary.tsv`, `classification.tsv`, and the phylogeny outputs together for reproducibility.
