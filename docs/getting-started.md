# Getting Started

## Environment

Ariadne is a Python command-line tool designed to work with:

- Python `>= 3.9`
- `mafft`
- `iqtree` or `iqtree2`
- `numpy >= 1.24`
- `openpyxl >= 3.1`
- `pyhmmer >= 0.12.0`
- `pyrodigal >= 3.7.0`
- `scikit-learn >= 1.4`

The recommended setup is the bundled conda environment.

## Installation

### Option 1. Conda

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

### Option 2. Local virtual environment

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## Repository inputs

The current Ariadne release assumes a tree-native project structure:

- `input/`
  protein FASTA inputs for standard discovery
- `tree/`
  multi-clade TPS reference FASTA collection
- `output/`
  legacy example outputs, not required by the active workflow

## First end-to-end run

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

This run will create:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
├── 04_phylogeny/
└── pipeline_summary.tsv
```

## Transcriptome mode

If your inputs are transcriptomes instead of predicted proteins:

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results_from_transcriptomes/
```

In this mode, Ariadne uses `pyrodigal` to predict ORFs before HMM screening.

## Local sanity checks

After installation, these commands should work:

```bash
ariadne --help
ariadne run --help
ariadne classify --help
```

If you installed Ariadne inside the project virtual environment, the equivalent commands are:

```bash
.venv/bin/python -m ariadne --help
.venv/bin/python -m ariadne run --help
```

## Design principles to keep in mind

- `tree/` is the default reference source across the whole workflow.
- A discovery query HMM is automatically built from the coral reference when `--query-hmm` is omitted.
- A TPS HMM library is automatically built from all FASTA files under `tree/` when `--tps-hmm-dir` is omitted.
- After classification, the workflow directly proceeds to alignment and phylogeny.

## Next steps

- Read [Method](method.md) for a deeper explanation of the four stages.
- Use [Tutorials](tutorials.md) to reproduce a practical end-to-end analysis.
- Keep [CLI Reference](cli-reference.md) open when tuning parameters.
