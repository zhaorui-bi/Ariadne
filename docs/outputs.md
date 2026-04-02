# Outputs

## Output root

A standard `ariadne run` creates:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
├── 04_phylogeny/
└── pipeline_summary.tsv
```

## `01_discovery`

Key files:

- `all_predicted_proteins.faa`
- `candidates.protein.faa`
- `candidates.orf.fna`
- `candidates.hits.tsv`

What to inspect:

- use `candidates.hits.tsv` to check which records passed HMM discovery
- use `candidates.protein.faa` as the starting point for downstream filtering

## `02_filtering`

Key files:

- `candidates.filtered.faa`
- `filter_report.tsv`
- `dedupe_clusters.tsv`
- `reference_matches.tsv`
- `manual_review.tsv`

What to inspect:

- `filter_report.tsv` explains the status of every candidate (`kept`, `removed`, `deduplicated_against`)
- `dedupe_clusters.tsv` records representative-member relationships
- `reference_matches.tsv` logs candidates whose sequences match one or more reference sequences — these candidates are **still kept** in `candidates.filtered.faa`; the file is for traceability only
- `candidates.filtered.faa` is the canonical handoff into classification and phylogeny

## `03_classification`

Key files:

- `tps_features.tsv`
- `embedding.tsv`
- `embedding.svg`
- `embedding_3d_sections.svg`
- `classification.tsv`
- `nearest_neighbors.tsv`
- `candidate_cluster_context.tsv`
- `global_context_tree.nwk`

What to inspect:

- `classification.tsv`
  best single-file summary of candidate predictions
- `nearest_neighbors.tsv`
  evidence for the predicted source assignment
- `embedding.svg`
  fast visual inspection of candidate placement
- `embedding_3d_sections.svg`
  publication-style multi-view embedding figure

### Optional CeeSs outputs (`TPS/TPS.xlsx` + ESM dependencies required)

When `--ceess-xlsx` points to a valid workbook and the optional ESM stack is installed, the following files are added:

- `ceess_predictions.tsv` — per coral-like candidate: `esm_type_prediction`, aggregated `P(CeeSs)`, and one `esm_type_probability_*` column per TPS type
- `ceess_candidates.tsv` — shortlisted candidates with `P(CeeSs) ≥ --ceess-threshold` (default 0.9)
- `ceess_candidates.fasta` — FASTA export of the shortlisted candidates
- `ceess_projection.tsv` — 2D LDA/PCA projected coordinates of the ESM embedding space (training references + candidates)
- `ceess_embedding.svg` — combined labeled-reference plus candidate projection figure
- `ceess_model_metrics.tsv` — cross-validated ESM classifier metrics on `TPS.xlsx` (accuracy, macro-F1, per-type precision/recall)
- `ceess_model_confusion_matrix.tsv` — multi-class confusion matrix over labeled TPS types
- `ceess_group_confusion_matrix.tsv` — binary CeeSs vs non-CeeSs confusion matrix
- `type_score_hits/` — per-type TSV of candidates scoring above 0.95 for each TPS type
- `type_score_fastas/` — per-type FASTA of the same high-confidence candidates

`classification.tsv` also carries the CeeSs columns for every coral-like candidate: `is_coral_like`, `esm_type_prediction`, `esm_ceess_label`, `esm_ceess_probability`, `is_ceess_candidate`, and the per-type `esm_type_probability_*` columns.

## `04_phylogeny`

Key files:

- `phylogeny_input.fasta`
- `phylogeny_alignment.fasta`
- `phylogeny_sequence_map.tsv`
- `iqtree.treefile`
- `iqtree.iqtree`
- `iqtree.log`
- `phylogeny_preview.svg`

What to inspect:

- `phylogeny_sequence_map.tsv`
  maps tree-safe identifiers back to original sequence headers
- `iqtree.treefile`
  final phylogeny in Newick format
- `iqtree.iqtree`
  model and inference summary from IQ-TREE
- `phylogeny_preview.svg`
  compact SVG preview automatically rendered from the final IQ-TREE tree

## Reading order

If you are new to Ariadne, this reading order usually works best:

1. `pipeline_summary.tsv`
2. `03_classification/classification.tsv`
3. `03_classification/embedding.svg`
4. `04_phylogeny/phylogeny_preview.svg`
5. `04_phylogeny/iqtree.treefile`
6. `04_phylogeny/iqtree.iqtree`

## Representative preview

<figure class="paper-figure">
  <img src="assets/latest_embedding.svg" alt="Bundled classification embedding preview">
  <figcaption>
    Figure 2. Classification embedding from a representative run. 100 candidates discovered → 36 retained after filtering → 36 classified as coral-like → 5 CeeSs candidates shortlisted (P(CeeSs) ≥ 0.9, ESM2-650M MLP, CV accuracy 77.6%).
  </figcaption>
</figure>
