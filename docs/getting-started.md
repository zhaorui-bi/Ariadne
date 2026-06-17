# Getting Started

Ariadne is designed to be reproducible in a research setting: one environment, one reference root, and one command that takes you from raw protein inputs to a phylogeny-ready result directory.

## Requirements

Core dependencies (installed automatically with the package):

- Python `≥ 3.9`
- `numpy ≥ 1.24`
- `openpyxl ≥ 3.1`
- `pyhmmer ≥ 0.12.0`
- `pyrodigal ≥ 3.7.0`
- `scikit-learn ≥ 1.4`

External tools required only for the phylogeny stage:

- `mafft`
- `iqtree` or `iqtree2`

Optional, for the ESM2 CeeSs classifier (`pip install -e '.[esm]'`):

- `torch ≥ 2.2`
- `transformers ≥ 4.44`

## Installation

### Conda (recommended)

The bundled Conda environment keeps the bioinformatics dependencies reproducible across platforms.

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

### Local virtual environment

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

## Project layout

The Ariadne workflow assumes a tree-native repository structure:

| Path | Role |
|---|---|
| `input/` | example protein FASTA inputs for discovery |
| `tree/` | the multi-clade TPS reference collection — the reference backbone of the whole pipeline |
| `TPS/` | optional labeled coral TPS workbook (`TPS.xlsx`) used to train the CeeSs classifier |
| `ariadne/hmm/` | bundled discovery and TPS-library HMMs built from the current `tree/` dataset |

The key conceptual point is that `tree/` is **not** just a phylogeny folder — it is the reference backbone reused by every stage.

## Your first run

The simplest complete run is:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

This executes the full four-stage workflow:

1. HMM-guided discovery;
2. coverage / length filtering and near-duplicate removal;
3. TPS feature-space classification;
4. MAFFT alignment and IQ-TREE reconstruction.

If `TPS/TPS.xlsx` is present and the ESM dependencies are installed, the classification stage additionally:

1. keeps the `coral-like` candidates from the HMM classification layer;
2. trains a small ESM2 type classifier on the labeled coral TPS workbook;
3. shortlists `cembrene A / cembrene B` candidates as the final CeeSs set.

The expected output layout is:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
├── 04_phylogeny/
└── pipeline_summary.tsv
```

## Transcriptome mode

If your starting point is transcriptome FASTA rather than predicted proteins, Ariadne predicts ORFs with Pyrodigal before HMM search:

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results_from_transcriptomes/
```

## Sanity checks

Once the environment is ready, these commands should succeed:

```bash
ariadne --help
ariadne run --help
ariadne classify --help
```

Inside the local project virtual environment, the equivalents are:

```bash
.venv/bin/python -m ariadne --help
.venv/bin/python -m ariadne run --help
```

## Defaults worth remembering

- If `--query-hmm` is omitted, Ariadne uses the bundled `ariadne/hmm/query.hmm`, falling back to building a discovery HMM from the coral reference in `tree/`.
- If `--tps-hmm-dir` is omitted, Ariadne uses the bundled `ariadne/hmm/` library, falling back to building a fresh TPS HMM library from `tree/`.
- Classification flows directly into alignment and phylogeny; use `--skip-phylogeny` to stop after classification.
- When `TPS/TPS.xlsx` is present and the ESM stack is installed, classification also writes `ceess_predictions.tsv`, `ceess_candidates.tsv`, and `ceess_candidates.fasta`.

## Suggested reading order

1. read [Method](method.md) to understand the four-stage design;
2. run the bundled example from [Tutorials](tutorials.md);
3. inspect `classification.tsv`, `embedding.svg`, and `iqtree.treefile`;
4. keep the [CLI Reference](cli-reference.md) open while tuning parameters.
