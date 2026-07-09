# Advanced Usage

This page collects the knobs that are useful after the first successful run. Keep the defaults for routine screening, then tune one layer at a time so changes remain interpretable.

## Manuscript-grade Run Template

Use a fixed output directory, explicit thresholds, an explicit workbook path, and a reproducible ESM seed:

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --skip-phylogeny \
  --discovery-min-score 40 \
  --discovery-max-evalue 1e-5 \
  --min-coverage 10 \
  --min-length 300 \
  --identity-threshold 0.95 \
  --top-k 5 \
  --ceess-classifier mlp \
  --ceess-threshold 0.9 \
  --ceess-random-state 0 \
  --output-dir results_manuscript/
```

Archive these files with the candidate FASTA:

- `pipeline_summary.tsv`
- `02_filtering/filter_report.tsv`
- `02_filtering/dedupe_clusters.tsv`
- `03_classification/classification.tsv`
- `03_classification/nearest_neighbors.tsv`
- `03_classification/embedding.svg`
- `03_classification/ceess_predictions.tsv` when ESM2 scoring is enabled

## Discovery Sensitivity

Discovery controls the size of the candidate universe. Loosen thresholds when screening distant homologs; tighten them when the output is dominated by weak profile hits.

| Parameter | When to tune |
| --- | --- |
| `--discovery-min-score` | raise it to remove weak bit-score hits |
| `--discovery-max-evalue` | lower it for stricter profile-HMM evidence |
| `--query-hmm` | provide a curated HMM instead of letting Ariadne build one from references |

Example:

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --query-hmm curated_tps_query.hmm \
  --discovery-min-score 55 \
  --discovery-max-evalue 1e-8 \
  --output-dir results_strict_discovery/
```

## Filtering And De-duplication

Filtering should match the quality of the upstream assembly or protein prediction. Low-coverage transcriptomes may need a lower coverage cutoff; isoform-heavy assemblies may need stricter deduplication.

```bash
ariadne filter \
  --input-fasta results/01_discovery/candidates.protein.faa \
  --reference-dir reference_fastas/ \
  --min-coverage 8 \
  --min-length 300 \
  --identity-threshold 0.98 \
  --output-dir results_filter_tuned/
```

Read `filter_report.tsv` before changing multiple thresholds at once.

## Reference-space Classification

`--top-k` controls how much neighbor evidence is reported per candidate. A smaller value highlights the nearest references; a larger value gives more context for ambiguous candidates.

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir reference_fastas/ \
  --top-k 10 \
  --tree-neighbors 20 \
  --output-dir results/03_classification_top10/
```

When you want full control over the profile library, prebuild it and pass it explicitly:

```bash
ariadne build-tps-hmm-library \
  --alignment coral=reference_fastas/coral.fasta insect=reference_fastas/insect.fasta \
  --output-dir curated_tps_hmms/

ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir reference_fastas/ \
  --tps-hmm-dir curated_tps_hmms/ \
  --output-dir results/03_classification_curated_hmms/
```

## ESM2 CeeSs Classifier Heads

The ESM2 layer supports three classifier modes:

| Mode | Use case |
| --- | --- |
| `mlp` | default supervised head for most analyses |
| `logreg` | lightweight baseline for reproduction and small runs |
| `contrastive` | representation-learning experiment using the Barlow Twins path |

Baseline logistic regression:

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --ceess-classifier logreg \
  --ceess-threshold 0.9 \
  --output-dir results/03_classification_logreg/
```

Contrastive representation path:

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --ceess-classifier contrastive \
  --ceess-barlow-representation-dim 128 \
  --ceess-barlow-projection-dim 64 \
  --ceess-barlow-redundancy-weight 0.005 \
  --ceess-threshold 0.9 \
  --output-dir results/03_classification_contrastive/
```

## Reusing An MLP Checkpoint

For repeated candidate screens with the same `TPS.xlsx` label space, reuse a saved MLP checkpoint:

```bash
ariadne classify \
  --candidates new_candidates.filtered.faa \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --ceess-classifier mlp \
  --ceess-mlp-checkpoint previous_results/03_classification/ceess_classifier_checkpoint.pt \
  --output-dir results_reused_classifier/
```

Checkpoint reuse is currently limited to `--ceess-classifier mlp`.

## Reproducibility Checklist

- record the exact `ariadne --version` output;
- keep the reference FASTA directory used for the run;
- keep `TPS.xlsx` with the exact labels used for ESM2 scoring;
- archive `pipeline_summary.tsv` and all stage-level TSV outputs;
- report any non-default thresholds in Methods or Supplementary Methods;
- inspect both `classification.tsv` and `nearest_neighbors.tsv` before selecting final candidates.
