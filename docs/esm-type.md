# ESM Type Model

## Motivation

The `TPS/TPS.xlsx` spreadsheet provides a small but highly relevant supervised dataset for coral TPS representation learning:

- column 1: `Name`
- column 2: `Protein`
- column 3: `Type`
- column 4: `Species`

Because the protein sequence and the type annotation are already paired in the same table, this dataset is a natural entry point for an ESM-based type-separation workflow.

In the current Ariadne design, this model has two roles:

- standalone analysis through `ariadne esm-type`
- an integrated second-stage CeeSs head inside `ariadne classify` and `ariadne run`

## Design choice

For the current dataset size, Ariadne uses a deliberately stable small-sample strategy:

1. extract mean-pooled sequence embeddings with ESM2
2. standardize the embedding vectors
3. use supervised `LDA` for 2D visualization
4. use multinomial `Logistic Regression` as the type classifier
5. evaluate performance with stratified cross-validation

This is preferable to end-to-end fine-tuning when the labeled dataset is still relatively small.

## Command

### Standalone mode

```bash
ariadne esm-type \
  --xlsx TPS/TPS.xlsx \
  --output-dir esm_results/
```

The default model is:

```text
facebook/esm2_t33_650M_UR50D
```

You can override it with either a full Hugging Face model id or a short preset such as `150M`, `650M`, `3B`, or `15B`:

```bash
ariadne esm-type \
  --xlsx TPS/TPS.xlsx \
  --model-name 150M \
  --output-dir esm_results_large/
```

### Integrated CeeSs mode

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_classification/
```

When `TPS/TPS.xlsx` is present and the optional ESM dependencies are installed, Ariadne automatically:

1. keeps only `coral-like` candidates from the normal HMM-based classification stage
2. trains the ESM classifier from the labeled coral TPS workbook
3. ranks those coral-like candidates by `P(cembrene A) + P(cembrene B)`
4. exports final CeeSs candidates to `ceess_candidates.tsv` and `ceess_candidates.fasta`

## Outputs

The command writes:

- `esm_embeddings.npz`
  compressed raw ESM embeddings
- `esm_projection.tsv`
  2D projected coordinates for visualization
- `esm_predictions.tsv`
  per-sequence true labels, cross-validated predictions, and confidences
- `esm_confusion_matrix.tsv`
  type-level confusion matrix
- `esm_metrics.tsv`
  aggregate and per-class metrics
- `esm_type_centroids.tsv`
  type-level centroid summary
- `esm_embedding.svg`
  2D SVG visualization of type separation

Inside the integrated classification workflow, Ariadne additionally writes:

- `ceess_predictions.tsv`
  ESM predictions for each coral-like candidate
- `ceess_candidates.tsv`
  final shortlisted CeeSs candidates
- `ceess_candidates.fasta`
  FASTA export of the shortlisted candidates
- `ceess_embedding.svg`
  combined labeled-reference plus candidate projection
- `ceess_model_metrics.tsv`
  cross-validated metrics on the `TPS.xlsx` training set

## Current baseline on `TPS/TPS.xlsx`

Using the current `TPS/TPS.xlsx` spreadsheet and the default ESM2 model, Ariadne produced the following baseline:

- number of proteins: `50`
- number of types: `5`
- cross-validation folds: `4`
- cross-validated accuracy: `0.74`
- cross-validated macro-F1: `0.639939`

The present baseline suggests that:

- `cembrene B` is highly separable
- `klysimplexin R` is also strongly separable
- the main confusion currently lies among `cembrene A`, `Sesquiterpenes`, and `xeniaphyllene`

This behavior is consistent with a realistic biological setting in which some product classes are much more compact in sequence space than others.

## Why this baseline matters

This ESM workflow is useful for two complementary reasons:

- it provides a representation-level baseline for product-type separation
- it offers an interpretable path toward future refinement, such as prototype-based ranking, metric learning, or fine-tuned sequence encoders

## Practical extensions

Strong next steps would include:

- class-balanced reweighting with stronger calibration
- hierarchical classification, for example first separating broad TPS groups and then resolving subtype labels
- prototype-distance ranking against curated CeeSs subsets
- replacing logistic regression with a metric-learning or nearest-prototype objective
