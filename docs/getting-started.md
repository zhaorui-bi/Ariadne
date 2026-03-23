# Getting Started

## Installation philosophy

Ariadne is intended to be easy to reproduce in a research setting: one environment, one reference root, and one command that can take you from raw protein inputs to a phylogeny-ready result directory.

The current release depends on:

- Python `>= 3.9`
- `mafft`
- `iqtree` or `iqtree2`
- `numpy >= 1.24`
- `openpyxl >= 3.1`
- `pyhmmer >= 0.12.0`
- `pyrodigal >= 3.7.0`
- `scikit-learn >= 1.4`

The recommended setup is the bundled conda environment, because it keeps the bioinformatics dependencies reproducible across platforms.

## Recommended installation

### Conda

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

### Local virtual environment

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## Minimal project assumptions

The active Ariadne workflow assumes a tree-native repository structure:

- `input/`
  protein FASTA inputs for standard discovery
- `tree/`
  the multi-clade TPS reference collection used across discovery, classification, and phylogeny
- `output/`
  a historical example-output directory that is no longer required by the current workflow

The key conceptual shift is that `tree/` is not just a phylogeny folder. It is the reference backbone of the whole pipeline.

## First end-to-end run

The simplest complete run is:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

This command executes the full four-stage workflow:

1. discovery by HMM-guided screening
2. filtering and near-duplicate removal
3. TPS feature-space classification
4. MAFFT alignment and IQ-TREE reconstruction

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

If your starting point is transcriptome FASTA rather than predicted proteins, Ariadne can infer ORFs before HMM search:

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results_from_transcriptomes/
```

In this mode, `pyrodigal` is used upstream of the discovery step.

## Sanity checks after installation

These commands should succeed once the environment is ready:

```bash
ariadne --help
ariadne run --help
ariadne classify --help
```

If you are working inside the local project virtual environment:

```bash
.venv/bin/python -m ariadne --help
.venv/bin/python -m ariadne run --help
.venv/bin/python -m ariadne classify --help
```

## Practical assumptions worth remembering

- If `--query-hmm` is omitted, Ariadne automatically builds a discovery HMM from the coral reference in `tree/`.
- If `--tps-hmm-dir` is omitted, Ariadne builds a TPS HMM library from all FASTA files under `tree/`.
- The active workflow proceeds directly from classification to alignment and phylogeny.
- The software is currently framed around coral TPS mining and CeeSs prioritization, but the underlying reference logic is still cross-clade.

## Suggested reading order for new users

If this is your first time using Ariadne, the most productive next sequence is:

1. read [Method](method.md) to understand the four-stage design
2. run the bundled example from [Tutorials](tutorials.md)
3. inspect `classification.tsv`, `embedding.svg`, and `iqtree.treefile`
4. keep [CLI Reference](cli-reference.md) open while tuning parameters
