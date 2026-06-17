# Outputs

A standard `ariadne run` produces one directory per stage plus a top-level summary:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
├── 04_phylogeny/
└── pipeline_summary.tsv
```

## `01_discovery`

| File | Contents |
|---|---|
| `all_predicted_proteins.faa` | all input/ORF-derived proteins searched |
| `candidates.protein.faa` | proteins that passed HMM discovery — the handoff into filtering |
| `candidates.orf.fna` | nucleotide ORFs (transcriptome mode) |
| `candidates.hits.tsv` | per-record HMM hit table (score, E-value) |

**Inspect first:** `candidates.hits.tsv` to confirm which records passed discovery.

## `02_filtering`

| File | Contents |
|---|---|
| `candidates.filtered.faa` | the canonical handoff into classification and phylogeny |
| `filter_report.tsv` | status of every candidate (`kept`, `removed`, `deduplicated_against`) |
| `dedupe_clusters.tsv` | representative ↔ member relationships from near-duplicate collapsing |
| `reference_matches.tsv` | candidates matching a reference sequence — **kept** in the filtered output; logged for traceability only |
| `manual_review.tsv` | borderline cases flagged for manual inspection |

## `03_classification`

| File | Contents |
|---|---|
| `classification.tsv` | best single-file summary of candidate predictions |
| `nearest_neighbors.tsv` | reference evidence for each predicted assignment |
| `tps_features.tsv` | per-sequence profile-HMM feature matrix |
| `embedding.tsv` | embedding coordinates |
| `embedding.svg` | 2-D visual inspection of candidate placement |
| `embedding_3d_sections.svg` | publication-style multi-view embedding figure |
| `candidate_cluster_context.tsv` | local cluster context per candidate |
| `global_context_tree.nwk` | global context tree (Newick) |

**Inspect first:** `classification.tsv` for predictions, then `nearest_neighbors.tsv` for the supporting references.

### Optional CeeSs outputs

When `--ceess-xlsx` points to a valid workbook and the ESM stack is installed, the following are added:

| File | Contents |
|---|---|
| `ceess_predictions.tsv` | per candidate: `esm_type_prediction`, aggregated `P(CeeSs)`, one `esm_type_probability_*` column per TPS type |
| `ceess_candidates.tsv` | shortlist with `P(CeeSs) ≥ --ceess-threshold` (default 0.9) |
| `ceess_candidates.fasta` | FASTA of the shortlist |
| `ceess_projection.tsv` | 2-D LDA/PCA coordinates of the ESM embedding (references + candidates) |
| `ceess_embedding.svg` | combined labeled-reference + candidate projection figure |
| `ceess_model_metrics.tsv` | cross-validated metrics on `TPS.xlsx` (accuracy, macro-F1, per-type precision/recall) |
| `ceess_model_confusion_matrix.tsv` | multi-class confusion matrix over TPS types |
| `ceess_group_confusion_matrix.tsv` | binary CeeSs vs non-CeeSs confusion matrix |
| `type_score_hits/` | per-type TSV of candidates scoring above 0.95 |
| `type_score_fastas/` | per-type FASTA of the same high-confidence candidates |

`classification.tsv` also carries the CeeSs columns for every coral-like candidate: `is_coral_like`, `esm_type_prediction`, `esm_ceess_label`, `esm_ceess_probability`, `is_ceess_candidate`, and the per-type `esm_type_probability_*` columns.

## `04_phylogeny`

| File | Contents |
|---|---|
| `phylogeny_input.fasta` | merged, deduplicated candidates + references |
| `phylogeny_alignment.fasta` | MAFFT alignment |
| `phylogeny_sequence_map.tsv` | maps tree-safe identifiers back to original headers |
| `iqtree.treefile` | final maximum-likelihood phylogeny (Newick) |
| `iqtree.iqtree` | IQ-TREE model and inference summary |
| `iqtree.log` | IQ-TREE run log |
| `phylogeny_preview.svg` | compact SVG preview rendered from the final tree |

## Suggested reading order

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
    Figure. Classification embedding from a representative run: 100 candidates discovered → 36 retained after filtering → 36 classified as coral-like → 5 CeeSs candidates shortlisted at <code>P(CeeSs) ≥ 0.9</code> (ESM2-650M, MLP head; cross-validated accuracy 77.6%).
  </figcaption>
</figure>
