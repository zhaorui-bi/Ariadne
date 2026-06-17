<p align="center">
  <img src="docs/logo.png" alt="Ariadne logo" width="520">
</p>

<h1 align="center">Ariadne</h1>

<p align="center">
  <strong>Terpene synthase discovery platform for coral TPS mining<br>and cembrene-class synthase (CeeSs) prioritization</strong>
</p>

<p align="center">
  <a href="./README_ZH.md">中文</a> &nbsp;·&nbsp;
  <a href="./docs/index.md">Documentation</a> &nbsp;·&nbsp;
  <a href="#quick-start">Quick Start</a> &nbsp;·&nbsp;
  <a href="#citation">Citation</a>
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/Python-%E2%89%A53.9-0b132b?style=flat-square&logo=python&logoColor=white">
  <img alt="Version" src="https://img.shields.io/badge/Version-1.1.0-0f766e?style=flat-square">
  <img alt="License" src="https://img.shields.io/badge/License-MIT-c2410c?style=flat-square">
  <img alt="Visualization" src="https://img.shields.io/badge/Stage%204-PCA%2FLDA%20visualization-2563eb?style=flat-square">
  <img alt="ESM2" src="https://img.shields.io/badge/Optional-ESM2%20CeeSs-7c3aed?style=flat-square">
  <a href="./docs/index.md"><img alt="Docs" src="https://img.shields.io/badge/Docs-MkDocs-1d4ed8?style=flat-square"></a>
</p>

---

## Abstract

Terpene synthases (TPSs) generate one of the largest and most chemically diverse families of natural products. Coral genomes encode many candidate TPS proteins, but assigning candidates to specific cembrene-class diterpenoid products remains difficult.

**Ariadne** is a four-stage platform for coral TPS discovery and **cembrene-class synthase (CeeSs)** prioritization. The current README workflow ends with candidate classification, PCA/LDA visualization, and optional ESM2 CeeSs scoring. It no longer treats a final phylogenetic tree as the last stage.

The current repository no longer ships a generated `data/` or HMM bundle. Provide a reference FASTA directory with `--reference-dir`, and use the optional root-level `TPS.xlsx` workbook for the ESM2 scoring layer. Users do **not** need to generate or maintain a separate HMM directory for the workflow described here.

## Contributions

- **One reference directory, four analysis stages.** A user-supplied reference FASTA directory supports discovery, filtering, classification, and visualization.
- **Interpretable classification.** Candidates are placed in a multi-clade TPS feature space, assigned nearest-reference labels, and reported with supporting neighbors rather than only a single opaque score.
- **PCA/LDA visual interpretation.** Stage 4 focuses on `embedding.svg`, `embedding_3d_sections.svg`, `embedding.tsv`, `embedding_variance.tsv`, and cluster-context tables for candidate triage.
- **CeeSs scoring with protein language models.** When `TPS.xlsx` and the optional ESM stack are available, Ariadne ranks coral-like candidates by `P(CeeSs)`.
- **Auditable outputs.** Each run emits TSV tables and SVG figures that track candidates from raw input through final prioritization.

## Method

Ariadne runs as a sequential four-stage workflow. The visualization stage is generated from the classification output, so it lives under `03_classification/` rather than a separate tree directory.

| Stage | Command | Input | Key output |
|---|---|---|---|
| 1 · Discovery | `ariadne run` | Transcriptomes / protein FASTAs | `candidates.protein.faa` |
| 2 · Filtering | `ariadne filter` or `ariadne run` | Candidate FASTA | `candidates.filtered.faa` |
| 3 · Classification | `ariadne classify` or `ariadne run` | Filtered FASTA + references | `classification.tsv`, `nearest_neighbors.tsv` |
| 4 · PCA/LDA visualization | `ariadne classify` or `ariadne run` | Classification feature space | `embedding.svg`, `embedding_3d_sections.svg`, `embedding_variance.tsv` |

**Stage 1 - Discovery.** Ariadne starts from transcriptome assemblies or predicted protein FASTAs. Transcriptome mode predicts ORFs with Pyrodigal; protein mode uses the supplied protein FASTAs directly.

**Stage 2 - Filtering.** Candidates are filtered by coverage and minimum length, and near-duplicates are collapsed at 95% identity. Reference matches are logged in `reference_matches.tsv` but retained in `candidates.filtered.faa`, which keeps known-like coral TPS alleles visible for downstream interpretation.

**Stage 3 - Classification.** Filtered candidates and reference sequences are represented in a multi-clade TPS feature space. Ariadne normalizes the feature matrix, assigns nearest-reference labels by *k*-nearest-neighbor voting, and writes the main prediction table.

**Stage 4 - PCA/LDA visualization.** The same feature space is projected for visual inspection. Ariadne prefers supervised LDA when the reference labels support it and falls back to PCA when needed. The primary artifacts are the 2-D embedding, the three-section 3-D view, variance table, and candidate cluster context.

## Installation

**Recommended** - Python 3.11 with the bundled Conda environment:

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

**Minimal** - core dependencies only:

```bash
pip install -e .
```

**With ESM2 CeeSs scoring:**

```bash
pip install -e '.[esm]'   # adds torch, transformers, tqdm
```

Core dependencies: `numpy >= 1.24`, `pyhmmer >= 0.12.0`, `pyrodigal >= 3.7.0`, `scikit-learn >= 1.4`, `openpyxl >= 3.1`.

## Quick Start

### Full README workflow from protein FASTAs

Place your predicted protein FASTA files in a directory such as `my_proteins/`, then run:

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --skip-phylogeny \
  --output-dir results/
```

`--skip-phylogeny` keeps `ariadne run` aligned with the current README workflow, where Stage 4 is PCA/LDA visualization rather than tree reconstruction.

Output layout:

```text
results/
├── 01_discovery/          # discovery hits and per-sample protein FASTAs
├── 02_filtering/          # filtered FASTA, filter_report.tsv, dedupe_clusters.tsv
├── 03_classification/     # classification, nearest neighbors, PCA/LDA SVGs
│                          #   (+ ceess_* outputs when TPS.xlsx is used)
└── pipeline_summary.tsv
```

### From transcriptome assemblies

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --skip-phylogeny \
  --output-dir results_from_transcriptomes/
```

### Single-stage classification and visualization

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir results/03_classification/
```

### From Python or a notebook

The pipeline stages are re-exported at the package top level, so Ariadne can be driven directly from a script or notebook. A bare `import ariadne` stays light because heavy optional dependencies are imported lazily on first use.

```python
import ariadne as ad

filtering = ad.filter_candidates(
    "candidate_tps.faa",
    "results/02_filtering",
    reference_dir="reference_fastas/",
)
classification = ad.classify_candidates(
    filtering["filtered_fasta"],
    reference_dir="reference_fastas/",
    output_dir="results/03_classification",
    ceess_xlsx="TPS.xlsx",
)
```

See [`tutorial/tutorial.ipynb`](tutorial/tutorial.ipynb) for a notebook walkthrough and [`tutorial/run_example.sh`](tutorial/run_example.sh) for a command-line smoke test.

## CeeSs Scoring (Optional)

When `TPS.xlsx` is present and the ESM stack is installed, `ariadne run` and `ariadne classify` attach an ESM2 scoring pass after feature-space classification:

1. Load labeled TPS sequences from `TPS.xlsx` (`Name` / `Protein` / `Type` / `Species`).
2. Compute frozen, mean-pooled ESM2 embeddings for training sequences and coral-like candidates.
3. Train a lightweight classifier head on the training embeddings.
4. Score each candidate; `P(CeeSs)` is the summed probability over CeeSs-positive type labels.
5. Write candidates above `--ceess-threshold` (default `0.9`) to `ceess_candidates.{tsv,fasta}`.

| `--ceess-classifier` | Head |
|---|---|
| `mlp` *(default)* | Torch MLP, cross-entropy + AdamW, class weighting |
| `logreg` | Scikit-learn logistic regression with standard scaling |
| `contrastive` | Barlow Twins projection network + MLP head |

Key files: `ceess_predictions.tsv`, `ceess_candidates.{tsv,fasta}`, `ceess_embedding.svg`, `ceess_model_metrics.tsv`. See [docs/esm-type.md](./docs/esm-type.md) for the full output schema.

## Command-line Interface

| Command | Purpose |
|---|---|
| `ariadne run` | End-to-end workflow through classification and PCA/LDA visualization |
| `ariadne filter` | Stage 2 - coverage / length filtering and deduplication |
| `ariadne classify` | Stages 3-4 - feature-space classification plus visualization |
| `ariadne prepare-references` | Build clean reference FASTA files from source data |

> Global flags (`--verbose`, `--log-file`) must precede the subcommand: `ariadne --verbose run ...`

The main workflow parameter tables live in the **[CLI Reference](./docs/cli-reference.md)**.

## Repository Layout

```text
Ariadne/
├── ariadne/           # core package
│   ├── search.py      # Stage 1 - candidate discovery
│   ├── filter.py      # Stage 2 - coverage, length, and near-duplicate filtering
│   ├── embed.py       # Stages 3-4 - classification and PCA/LDA visualization
│   ├── model.py       # optional ESM2 CeeSs scoring
│   ├── data.py        # reference data management
│   ├── utils.py       # logging, FASTA/TSV I/O, sequence utilities
│   └── cli.py         # command-line interface
├── TPS.xlsx           # optional CeeSs training workbook
├── docs/              # documentation site
├── docs/logo.png      # README and documentation logo
├── tutorial/          # runnable script and notebook tutorial
├── environment.yml
└── pyproject.toml
```

## Documentation

The full documentation site is in [`docs/`](./docs/index.md):

- [Getting Started](./docs/getting-started.md) - installation and your first run
- [Method](./docs/method.md) - the four-stage design
- [Tutorials](./docs/tutorials.md) - practical analysis pathways
- [CLI Reference](./docs/cli-reference.md) - every command and parameter
- [Outputs](./docs/outputs.md) - every artifact the pipeline produces
- [CeeSs Classifier](./docs/esm-type.md) - the ESM2 scoring layer
- [Citation](./docs/citation.md)

## Example

A minimal local smoke test can be run with:

```bash
bash tutorial/run_example.sh results/ my_proteins/ reference_fastas/
RUN_CEESS=1 bash tutorial/run_example.sh results/ my_proteins/ reference_fastas/
```

The script requires an input protein directory and a reference FASTA directory, and keeps the workflow focused on classification plus PCA/LDA visualization.

## Development

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
ruff check ariadne
```

## Citation

If Ariadne is useful in your work, please cite it:

```bibtex
@software{jiang2026ariadne,
  author  = {Jiang, Zhaorui},
  title   = {Ariadne: A Coral-Centered Terpene Synthase Discovery and CeeSs Prioritization Platform},
  year    = {2026},
  url      = {https://github.com/zhaorui-bi/Ariadne},
  version = {1.1.0}
}
```

## License

Released under the [MIT License](./LICENSE).
