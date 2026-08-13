# Tutorials

These tutorials mirror the current README workflow: run the pipeline, read the classification layer, and use PCA/LDA visualization for candidate triage.

## 1 · Run The Local Smoke Test

The repository ships a small script that runs the README workflow when you provide protein inputs and reference FASTAs:

```bash
bash tutorial/run_example.sh results/ my_proteins/ reference_fastas/
```

Expected stage directories:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
└── pipeline_summary.tsv
```

## 2 · Run Filtering And Classification Separately

When you already have a candidate TPS FASTA, you can inspect filtering and classification without rerunning the full pipeline.

```bash
# Filtering
ariadne filter \
  --input-fasta candidate_tps.faa \
  --reference-dir reference_fastas/ \
  --output-dir tmp_manual/02_filtering

# Classification + PCA/LDA visualization
ariadne classify \
  --candidates tmp_manual/02_filtering/candidates.filtered.faa \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir tmp_manual/03_classification
```

## 3 · Read Classification As The Screening Layer

Start with:

- `classification.tsv` - predicted source assignment for each candidate;
- `nearest_neighbors.tsv` - references that support each assignment;
- `embedding.svg` - global two-dimensional PCA/LDA view;
- `embedding_3d_sections.svg` - three-panel LDA/PCA summary (same layout as Figure S2);
- `ceess_predictions.tsv`, `ceess_candidates.tsv`, `ceess_candidates.fasta` - optional CeeSs shortlist.

A productive reading order:

1. inspect `classification.tsv`;
2. check `nearest_neighbors.tsv`;
3. inspect `embedding.svg`, then `embedding_3d_sections.svg`;
4. read `ceess_predictions.tsv` when ESM scoring is enabled;
5. open `ceess_candidates.tsv` / `.fasta` for the final shortlist.

## 4 · Run The Standalone Coral TPS ESM Model

To inspect the supervised coral TPS type model independently of discovery, call the Python API directly:

```python
import ariadne as ad

outputs = ad.analyze_tps_types_with_esm(
    "TPS.xlsx",
    output_dir="tmp_esm_results",
    classifier_kind="mlp",
)
print(outputs)
```

Files worth inspecting:

- `esm_embedding.svg` - labeled projection of the ESM embedding space;
- `esm_metrics.tsv` - cross-validated accuracy and macro-F1;
- `esm_confusion_matrix.tsv` - multi-class confusion matrix;
- `esm_predictions.tsv` - per-record cross-validated predictions.

## 5 · Iterate Quickly

For parameter tuning, keep the run focused on classification and visualization:

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir tmp_fast_iteration/
```

This gives repeated access to `classification.tsv`, `nearest_neighbors.tsv`, and `embedding.svg`.

## Suggested Workflow For Manuscript Preparation

1. run `ariadne run`;
2. screen candidates using `classification.tsv` and `embedding.svg`;
3. inspect the optional CeeSs shortlist;
4. archive `pipeline_summary.tsv`, `classification.tsv`, neighbor evidence, and visualization SVGs (`embedding.svg`, `embedding_3d_sections.svg`).
