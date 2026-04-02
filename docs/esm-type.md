# CeeSs Classifier (`ariadne.model`)

## Motivation

The `TPS/TPS.xlsx` spreadsheet provides a small but highly relevant supervised dataset for coral TPS representation learning:

- column 1: `Name`
- column 2: `Protein`
- column 3: `Type`
- column 4: `Species`

Because the protein sequence and the type annotation are already paired in the same table, this dataset is a natural entry point for an ESM-based CeeSs classification workflow.

In the current Ariadne design, this module (`ariadne.model`) is used as an integrated second-stage CeeSs head inside `ariadne classify` and `ariadne run`.

## Design

### Multi-class CeeSs head (integrated mode)

The integrated CeeSs head uses a deliberately simple small-sample strategy:

1. compute mean-pooled ESM2 embeddings for the labeled coral TPS proteins
2. keep the ESM backbone frozen and feed the embeddings into a lightweight classifier head
3. train either:
   - a small MLP with configurable epochs, learning rate, dropout, and batch size
   - or the legacy class-balanced logistic regression head for reproduction
4. use the resulting multi-class probabilities in two ways:
   - `esm_type_prediction` reports the top predicted TPS type directly
   - `esm_ceess_probability` is computed by summing the probabilities of all labels marked as CeeSs-positive in the workbook

### Workbook annotation for non-canonical CeeSs

If your dataset contains experimentally validated CeeSs that are not `cembrene A` or `cembrene B`, add an explicit workbook column such as `CeeSs_group` or `CeeSs`.

- rows marked `yes` / `true` / `1` / `CeeSs` are treated as CeeSs-positive
- rows marked `no` / `false` / `0` / `non-CeeSs` are treated as CeeSs-negative
- if the column is absent, Ariadne falls back to the legacy rule that only `cembrene A` and `cembrene B` are positive

Because the Cembrene A / B probabilities now come from the same multi-class head, there is no separate second-stage subtype model in integrated mode.

## Integrated CeeSs Mode

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_classification/
```

When `TPS/TPS.xlsx` is present and the optional ESM dependencies are installed, Ariadne automatically:

1. keeps only `coral-like` candidates from the normal HMM-based classification stage
2. trains a multi-class ESM-based classifier on the labeled coral TPS workbook
3. predicts a detailed TPS type for every coral-like candidate
4. sums the probabilities of all workbook-defined CeeSs-positive labels into `P(CeeSs)` and applies `--ceess-threshold`
5. exports final CeeSs candidates to `ceess_candidates.tsv` and `ceess_candidates.fasta`

## Outputs

### Integrated classification mode

Inside `ariadne classify` / `ariadne run`, the CeeSs head additionally writes:

- `ceess_predictions.tsv` — per coral-like candidate: `esm_type_prediction`, `esm_ceess_label`, `esm_ceess_probability`, `is_ceess_candidate`, plus one `esm_type_probability_*` column per TPS type
- `ceess_candidates.tsv` — shortlisted candidates with `P(CeeSs) ≥ --ceess-threshold` (default 0.9)
- `ceess_candidates.fasta` — FASTA export of the shortlisted candidates
- `ceess_projection.tsv` — 2D LDA/PCA projected coordinates (references + candidates)
- `ceess_embedding.svg` — combined labeled-reference plus candidate projection figure
- `ceess_model_metrics.tsv` — cross-validated metrics on `TPS.xlsx`: accuracy, macro-F1, per-type precision/recall/F1
- `ceess_model_confusion_matrix.tsv` — multi-class confusion matrix over TPS types
- `ceess_group_confusion_matrix.tsv` — binary CeeSs vs non-CeeSs confusion matrix
- `type_score_hits/` — per-type TSV of candidates above 0.95 probability
- `type_score_fastas/` — per-type FASTA of the same candidates

The `classification.tsv` produced by stage 3 also carries the CeeSs fields for every coral-like candidate:

| column | description |
| --- | --- |
| `is_coral_like` | `yes` / `no` — did the HMM feature space assign this candidate near coral references |
| `esm_type_prediction` | top predicted TPS type from the multi-class ESM head |
| `esm_ceess_label` | compatibility column carrying the same top predicted TPS type |
| `esm_ceess_probability` | aggregated P(CeeSs) computed from all workbook-defined positive labels |
| `esm_type_probability_*` | one probability column per TPS type, for example `esm_type_probability_cembrene_a` or `esm_type_probability_klysimplexin_r` |
| `is_ceess_candidate` | `yes` if `esm_ceess_probability ≥ --ceess-threshold` |

## Python API

The core integrated workflow is importable from `ariadne.model`:

```python
from ariadne.model import (
    classify_ceess_candidates_with_esm,
    compute_esm_embeddings,
    load_tps_xlsx,
    DEFAULT_ESM_MODEL_NAME,
)
```
