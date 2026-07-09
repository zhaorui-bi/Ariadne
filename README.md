<p align="center">
  <img src="docs/logo.png" alt="Ariadne logo" width="520">
</p>

<h1 align="center">Ariadne</h1>

<p align="center">
  <strong>Terpene synthase discovery, CeeSs prioritization, and PCA/LDA visualization</strong>
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
  <img alt="Docs" src="https://img.shields.io/badge/Docs-MkDocs-1d4ed8?style=flat-square">
</p>

---

## Overview

**Ariadne** is a research-grade Python package and CLI for terpene synthase (TPS) candidate discovery. It turns transcriptome or protein FASTA inputs into a traceable workflow:

```text
Discovery -> Filtering -> Classification -> PCA/LDA visualization
```

Optional ESM2 scoring can prioritize coral-like **cembrene-class synthase (CeeSs)** candidates when `TPS.xlsx` and the `[esm]` dependencies are available.

## Key Features

- HMM-based discovery from transcriptomes or predicted protein FASTAs.
- Coverage, length, and near-duplicate filtering with TSV audit reports.
- Reference-space TPS classification with nearest-neighbor evidence.
- PCA/LDA SVG projections for candidate triage.
- Optional ESM2 CeeSs scoring and FASTA export.

## Installation

Recommended:

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

Minimal:

```bash
pip install -e .
```

With ESM2 scoring:

```bash
pip install -e '.[esm]'
```

## Quick Start

Prepare a reference FASTA directory, then run:

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --skip-phylogeny \
  --output-dir results/
```

Transcriptome mode:

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir reference_fastas/ \
  --skip-phylogeny \
  --output-dir results_from_transcriptomes/
```

Core output layout:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
└── pipeline_summary.tsv
```

## Documentation

- [Getting Started](./docs/getting-started.md)
- [Method](./docs/method.md)
- [Advanced Usage](./docs/advanced-usage.md)
- [Software Architecture](./docs/software-architecture.md)
- [CLI Reference](./docs/cli-reference.md)
- [Outputs](./docs/outputs.md)
- [CeeSs Classifier](./docs/esm-type.md)

## Development

```bash
pip install -e '.[dev]'
python -m ariadne --help
mkdocs build --strict
```

## Citation

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
