# Outputs

A standard README-aligned `ariadne run` produces discovery, filtering, and classification/visualization directories plus a top-level summary:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
└── pipeline_summary.tsv
```

Use `--skip-phylogeny` with `ariadne run` to keep the output layout aligned with the current four-stage workflow.

## `01_discovery`

| File | Contents |
|---|---|
| `all_predicted_proteins.faa` | all input/ORF-derived proteins searched |
| `candidates.protein.faa` | proteins that passed discovery and move into filtering |
| `candidates.orf.fna` | nucleotide ORFs in transcriptome mode |
| `candidates.hits.tsv` | per-record discovery hit table |

**Inspect first:** `candidates.hits.tsv` to confirm which records passed discovery.

## `02_filtering`

| File | Contents |
|---|---|
| `candidates.filtered.faa` | canonical handoff into classification and visualization |
| `filter_report.tsv` | status of every candidate (`kept`, `removed`, `deduplicated_against`) |
| `dedupe_clusters.tsv` | representative to member relationships from near-duplicate collapsing |
| `reference_matches.tsv` | candidates matching a reference sequence, kept in the filtered output and logged for traceability |
| `manual_review.tsv` | borderline cases flagged for manual inspection |

## `03_classification`

| File | Contents |
|---|---|
| `classification.tsv` | best single-file summary of candidate predictions |
| `nearest_neighbors.tsv` | reference evidence for each predicted assignment |
| `tps_features.tsv` | per-sequence feature matrix |
| `embedding.tsv` | PCA/LDA coordinates and labels |
| `embedding_variance.tsv` | variance / projection metadata |
| `embedding.svg` | 2-D visual inspection plot |
| `embedding_3d_sections.svg` | publication-style multi-view embedding figure |
| `candidate_cluster_context.tsv` | local cluster context per candidate |
| `assignment_summary.tsv` | count and confidence summary by predicted source |

**Inspect first:** `classification.tsv`, then `nearest_neighbors.tsv`, then `embedding.svg`.

### Optional CeeSs Outputs

When `--ceess-xlsx TPS.xlsx` points to a valid workbook and the ESM stack is installed, the following are added:

| File | Contents |
|---|---|
| `ceess_predictions.tsv` | per candidate: predicted TPS type, aggregated `P(CeeSs)`, and per-type probabilities |
| `ceess_candidates.tsv` | shortlist with `P(CeeSs) >= --ceess-threshold` (default `0.9`) |
| `ceess_candidates.fasta` | FASTA of the shortlist |
| `ceess_projection.tsv` | 2-D LDA/PCA coordinates of the ESM embedding |
| `ceess_embedding.svg` | labeled-reference plus candidate projection figure |
| `ceess_model_metrics.tsv` | cross-validated metrics on `TPS.xlsx` |
| `ceess_model_confusion_matrix.tsv` | multi-class confusion matrix over TPS types |
| `ceess_group_confusion_matrix.tsv` | binary CeeSs vs non-CeeSs confusion matrix |
| `type_score_hits/` | per-type TSV of high-confidence candidates |
| `type_score_fastas/` | per-type FASTA of high-confidence candidates |

`classification.tsv` also carries the CeeSs columns for every coral-like candidate: `is_coral_like`, `esm_type_prediction`, `esm_ceess_label`, `esm_ceess_probability`, `is_ceess_candidate`, and the per-type `esm_type_probability_*` columns.

## Suggested Reading Order

1. `pipeline_summary.tsv`
2. `03_classification/classification.tsv`
3. `03_classification/nearest_neighbors.tsv`
4. `03_classification/embedding.svg`
5. `03_classification/embedding_3d_sections.svg`
6. optional `03_classification/ceess_candidates.tsv`
