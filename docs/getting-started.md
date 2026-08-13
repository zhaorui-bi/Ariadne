# Getting Started

Ariadne is designed for a reproducible research workflow: one environment, one reference FASTA directory, and one command that takes protein or transcriptome inputs through classification and PCA/LDA visualization.

## Requirements

Core dependencies installed with the package:

- Python `>= 3.9`
- `numpy >= 1.24`
- `openpyxl >= 3.1`
- `pyhmmer >= 0.12.0`
- `pyrodigal >= 3.7.0`
- `scikit-learn >= 1.4`

Optional, for the ESM2 CeeSs classifier (`pip install -e '.[esm]'`):

- `torch >= 2.2`
- `transformers >= 4.44`

No separate HMM preparation step is required for the documented workflow.

## Installation

### Conda (Recommended)

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

### Local Virtual Environment

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

To enable the ESM-based CeeSs head during `classify` and `run`:

```bash
python -m pip install -e '.[esm]'
```

## Project Layout

| Path | Role |
|---|---|
| `reference_fastas/` | user-provided reference FASTA directory |
| `TPS.xlsx` | optional CeeSs training workbook |
| `tutorial/` | runnable script and notebook tutorial |
| `docs/logo.png` | current logo |
| `ariadne/` | Python package and CLI implementation |

The reference FASTA directory is input data, not a generated HMM output directory.

## Your First Run

Place predicted protein FASTA files in `my_proteins/`, then run:

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir results/
```

This executes the README-aligned workflow:

1. candidate discovery;
2. coverage / length filtering and near-duplicate removal;
3. TPS feature-space classification;
4. PCA/LDA visualization under `03_classification/` (`embedding.svg` and the three-panel `embedding_3d_sections.svg`).

If `TPS.xlsx` is present and the ESM dependencies are installed, the classification stage also writes CeeSs predictions and candidate shortlists.

Expected output layout:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
└── pipeline_summary.tsv
```

## Transcriptome Mode

If your starting point is transcriptome FASTA rather than predicted proteins, Ariadne predicts ORFs with Pyrodigal before candidate screening:

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir results_from_transcriptomes/
```

## Sanity Checks

```bash
ariadne --help
ariadne run --help
ariadne classify --help
```

Inside the local project virtual environment:

```bash
.venv/bin/python -m ariadne --help
.venv/bin/python -m ariadne run --help
```

## Defaults Worth Remembering

- `--reference-dir reference_fastas/` points Ariadne at your reference FASTA directory.
- When `TPS.xlsx` is present and the ESM stack is installed, classification also writes `ceess_predictions.tsv`, `ceess_candidates.tsv`, and `ceess_candidates.fasta`.

## Suggested Reading Order

1. read [Method](method.md) to understand the four-stage design;
2. run the smoke test from [Tutorials](tutorials.md);
3. inspect `classification.tsv`, `nearest_neighbors.tsv`, and `embedding.svg`;
4. use [Advanced Usage](advanced-usage.md) when tuning thresholds or ESM2 classifier heads;
5. keep the [CLI Reference](cli-reference.md) open while tuning parameters.
