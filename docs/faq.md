# FAQ

## What Should I Use For `--reference-dir`?

Use a directory containing your curated multi-clade TPS reference FASTAs. The optional root-level `TPS.xlsx` workbook is used for CeeSs scoring when ESM dependencies are installed.

The reference FASTA directory is input data. Users do not need to generate a separate HMM directory from it for the documented workflow.

## What Inputs Does Ariadne Accept?

Either predicted proteins (`--protein-folder`, recursive `.faa` / `.fa` / `.fasta`) or transcriptome assemblies (`--transcriptomes`), in which case ORFs are predicted with Pyrodigal.

## What If I Only Want Classification And Visualization?

```bash
ariadne classify \
  --candidates my_candidates.faa \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir results_classification/
```

## Do I Need A GPU For CeeSs Scoring?

No. The ESM2 backbone can run on CPU, just more slowly. Select a device explicitly with `--ceess-device cuda:0` or `--ceess-device cpu`; by default Ariadne auto-detects.

## Which Outputs Should I Inspect First?

For most users: `classification.tsv`, then `nearest_neighbors.tsv`, then `embedding.svg` and `embedding_3d_sections.svg`. The three-panel SVG uses the same visual grammar as [Figure S2](method.md#stage-iv-pcalda-visualization). If ESM scoring is enabled, inspect `ceess_candidates.tsv` after that.

## Does Ariadne Still Use `Alignment.fasta`, Motif Analysis, Or Benchmark Mode?

No. The current release focuses on candidate discovery, filtering, classification, PCA/LDA visualization, and optional CeeSs prioritization.
