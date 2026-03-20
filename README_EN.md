<!-- # Ariadne -->

<p align="center">
  <img src="fig/ariadne_cvpr_logo.svg" alt="Ariadne CVPR-style logo" width="980">
</p>

<!-- <p align="center">
  <img src="fig/ariadne_icon_minimal.svg" alt="Ariadne minimal icon" width="92">
</p> -->

<p align="center">
  <a href="./README.md">中文</a>
</p>

<p align="center">
  🧬 <strong>Coral TPS / CeSS Discovery Platform</strong><br>
  🔬 HMM-guided mining · 🧠 motif calibration · 🌊 coral-aware embedding · 🌳 phylogeny-ready outputs
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/Python-3.11%2B-0F172A?style=flat-square&logo=python&logoColor=white">
  <img alt="MAFFT" src="https://img.shields.io/badge/MAFFT-required-0F766E?style=flat-square">
  <img alt="IQ-TREE" src="https://img.shields.io/badge/IQ--TREE-required-1D4ED8?style=flat-square">
  <img alt="CeSS" src="https://img.shields.io/badge/Target-CeSS%20Mining-F59E0B?style=flat-square">
  <img alt="Benchmark" src="https://img.shields.io/badge/Benchmark-Opt--in-7C3AED?style=flat-square">
</p>

> 🧵 Ariadne is designed not merely to collect generic TPS hits, but to narrow coral-centered datasets into a smaller, experiment-ready set of **CeSS-like / cembrene synthase** candidates.

## ✨ What Ariadne Does

Ariadne is a command-line platform for coral-centered terpene synthase discovery. The workflow includes:

1. `discovery`
   Automatically builds HMMs from `tree/` references and discovers candidate TPS sequences from proteins or transcriptomes
2. `filtering`
   Performs coverage, length, and deduplication filtering
3. `motif`
   Applies the TPS motif gate and combines cembrene-family plus validated CeSS references to decide `CeSS-like` candidates
4. `classification`
   Runs nearest-neighbor classification, clade-aware embedding, and context-tree generation in TPS HMM feature space
5. `phylogeny`
   Builds phylogenetic trees with MAFFT + IQ-TREE
6. `benchmark`
   Compares predictions against an expected FASTA set

## 🚀 Current Defaults

The current defaults are tuned for routine mining instead of mandatory benchmarking:

- ✅ Benchmark is **off by default**
- ✅ Center fallback is **on by default**
- ✅ `tree/` is used to automatically build:
  - `query.hmm`
  - TPS HMM library
  - motif reference inputs
  - phylogeny reference inputs
- ✅ CeSS-oriented calibrated motif thresholds are preserved by default
- ✅ Benchmarking is now an explicit **opt-in**

In practice, most users only need:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

## 🗂️ Repository Layout

```text
Ariadne/
├── ariadne/                # core source code
├── input/                  # analysis input
├── output/                 # expected output / benchmark examples
├── tree/                   # multi-clade TPS references and default HMM source
├── fig/                    # logo / figures
├── environment.yml         # conda environment
├── install.sh
└── pyproject.toml
```

### 📁 Data Convention

- `input/`
  Main analysis input directory. The current example is `input/all_tps.fasta`
- `tree/`
  The most important reference directory. It is used for HMM building, motif analysis, classification, and phylogeny
- `output/`
  Expected output example directory. The repository keeps the current filename `output/outpu.fasta`

## 🧩 Core Modules

| Module                      | Role                                                                                        |
| --------------------------- | ------------------------------------------------------------------------------------------- |
| `ariadne/cli.py`            | CLI entrypoint for `run`, `discover`, `motif`, `classify`, `phylogeny`, and `compare-fasta` |
| `ariadne/discovery.py`      | HMM construction, candidate discovery, ORF prediction                                       |
| `ariadne/filtering.py`      | Coverage, length, and deduplication filters                                                 |
| `ariadne/motif.py`          | TPS motif gate, CeSS assignment, priority ranking                                           |
| `ariadne/classification.py` | Feature matrix, embedding, nearest-neighbor classification, context tree                    |
| `ariadne/phylogeny.py`      | MAFFT + IQ-TREE phylogeny                                                                   |
| `ariadne/benchmark.py`      | FASTA benchmarking                                                                          |
| `ariadne/fasta_utils.py`    | FASTA / TSV I/O and cleanup                                                                 |

## 🛠️ Installation

### 1. Recommended: conda

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
```

### 2. venv

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
bash install.sh
```

### 3. Manual

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

### 📦 Runtime Dependencies

- Python `>= 3.9`
- `mafft`
- `iqtree`
- `numpy >= 1.24`
- `openpyxl >= 3.1`
- `pyhmmer >= 0.12.0`
- `pyrodigal >= 3.7.0`
- `scikit-learn >= 1.4`

Recommended Python version: `3.11`.

## ⚙️ Default Parameter Profile

These defaults are currently tuned for routine coral / CeSS mining:

| Parameter                           | Current Default | Meaning                                                                                |
| ----------------------------------- | --------------- | -------------------------------------------------------------------------------------- |
| Benchmark                           | `off`           | `05_benchmark/` is only produced when `--enable-benchmark` is passed                   |
| Center fallback                     | `on`            | If a direct anchor is missing, Ariadne falls back to a center window for CeSS judgment |
| `validated_cess_identity_threshold` | `0.95`          | High-confidence validated CeSS support threshold                                       |
| `cembrene_identity_threshold`       | `0.75`          | Family-level cembrene support threshold                                                |
| `cembrene_margin_threshold`         | `0.10`          | Minimum identity margin over non-cembrene references                                   |
| IQ-TREE fast mode                   | `on`            | Runs the full pipeline with `LG + --fast` by default                                   |

## 🧪 Quick Start

### 1. Daily Mining Run

This is the recommended daily entrypoint: no benchmark by default, center fallback on, and automatic tree-native model building.

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

### 2. CeSS-Calibrated Run

If you already have validated CeSS references, this is the recommended calibrated mode:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --reference-alignment tree/coral.fasta \
  --validated-cess-fasta output/outpu.fasta \
  --output-dir results_calibrated/
```

### 3. Benchmark Run

Only enable benchmark when you explicitly want method-to-reference comparison:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --reference-alignment tree/coral.fasta \
  --validated-cess-fasta output/outpu.fasta \
  --expected-fasta output/outpu.fasta \
  --enable-benchmark \
  --output-dir results_benchmark/
```

### 4. Transcriptome Input

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results_from_transcriptome/
```

## 🔍 Stage-by-Stage Usage

### `discover`

```bash
ariadne discover \
  --protein-folder input/ \
  --hmm query.hmm \
  --output-dir 01_discovery/
```

### `filter`

```bash
ariadne filter \
  --input-fasta 01_discovery/candidates.protein.faa \
  --output-dir 02_filtering/
```

### `motif`

```bash
ariadne motif \
  --candidates 02_filtering/candidates.filtered.faa \
  --reference-alignment tree/coral.fasta \
  --validated-cess-fasta output/outpu.fasta \
  --output-dir 04_motif/
```

Common tuning knobs:

- `--disable-center-fallback`
  Switches back to strict mode by disabling center fallback
- `--validated-cess-identity-threshold`
  Adjusts validated CeSS support strength
- `--cembrene-identity-threshold`
  Adjusts the family-level identity floor
- `--cembrene-margin-threshold`
  Adjusts the separation margin against non-cembrene references

### `classify`

```bash
ariadne classify \
  --candidates 02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir 03_classification/
```

### `phylogeny`

```bash
ariadne phylogeny \
  --candidates 02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir 06_phylogeny/
```

### `compare-fasta`

```bash
ariadne compare-fasta \
  --predicted-fasta 04_motif/cess_candidates.fasta \
  --expected-fasta output/outpu.fasta \
  --output-dir 05_benchmark/
```

## 📤 Output Map

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
│   ├── classification.tsv
│   ├── cess_priority_ranking.tsv
│   ├── embedding.tsv
│   ├── embedding.svg
│   ├── embedding_3d_sections.svg
│   ├── global_context_tree.nwk
│   └── _auto_tps_hmms/
├── 04_motif/
│   ├── motif_summary.tsv
│   ├── targeted_mining_summary.tsv
│   ├── cess_candidates.tsv
│   ├── cess_candidates.fasta
│   ├── cess_priority_ranking.tsv
│   └── motif_windows.svg
├── 05_benchmark/          # only created with --enable-benchmark
└── 06_phylogeny/
```

### 🎯 Most Important Outputs

- `04_motif/targeted_mining_summary.tsv`
  Overall candidate counts, TPS+/TPS-, and final CeSS-like counts
- `04_motif/cess_candidates.fasta`
  Final CeSS-like candidate sequences
- `03_classification/cess_priority_ranking.tsv`
  Best table for experimental prioritization
- `03_classification/embedding.svg`
  Overview of the current clade-aware embedding
- `06_phylogeny/iqtree.treefile`
  Final phylogeny tree output

## 🧠 Embedding Design

The current classification view is not a crowded vanilla PCA plot. It prioritizes:

- `lda_coral_subclade_spread`
  Supervised separation of reference clades and candidate groups
- An additional split of `coral` into:
  - `cembrene_a`
  - `cembrene_b`
  - `cembrene_c`
  - multiple `tps_subclade_*`

As a result, the newer `embedding.svg` is better suited for reading:

- whether reference clades are separated
- whether CeSS / TPS subclusters within coral are separated
- where candidate TPS+/TPS-/CeSS-like groups are located

## 🧵 CeSS Calibration Logic

The current CeSS calling logic is:

- If motif-window identity to a validated CeSS reference is `>= 0.95`
  the candidate is directly supported as `CeSS-like`
- If a direct anchor is found and:
  - identity to cembrene references is `>= 0.75`
  - and the margin over non-cembrene references is `>= 0.10`
    then the candidate is also supported as `CeSS-like`
- If only center fallback is available without validated CeSS support
  the sequence is not directly called CeSS-like, but it is retained in `cess_priority_ranking.tsv`

## 📊 Real Data Run in This Repository

This repository already includes a full real-data run:

- 📂 [tmp_run_calibrated_v3/](tmp_run_calibrated_v3)
- 🖼️ [tmp_run_calibrated_v3/03_classification/embedding.svg](tmp_run_calibrated_v3/03_classification/embedding.svg)
- 🏁 [tmp_run_calibrated_v3/05_benchmark/benchmark_summary.tsv](tmp_run_calibrated_v3/05_benchmark/benchmark_summary.tsv)

Key summary:

- `36` input candidates
- `32` TPS-positive
- `4` TPS-negative
- `4` CeSS-like
- `3` exact matches in benchmark mode
- plus `1` high-similarity near hit

## 🌳 HMM / Tree Notes

- `tree/` is now the default reference entrypoint instead of `Alignment.fasta`
- If `tree/*.fasta` are not MSAs, Ariadne automatically runs MAFFT before `hmmbuild`
- `run` and `phylogeny` default to IQ-TREE with `LG + --fast`
- If you only want mining and not benchmarking, you do not need `--expected-fasta`

## 🧰 Legacy Compatibility

These legacy scripts are still present, but the main workflow should now use the CLI:

- `ariadne/filter_contigs.py`
- `ariadne/filter_contigs_long.py`
- `ariadne/filter_coverage.py`
- `ariadne/DupRemover.py`
- `ariadne/visualization.py`

Preferred entrypoints:

- `ariadne run`
- `ariadne discover`
- `ariadne filter`
- `ariadne motif`
- `ariadne classify`
- `ariadne phylogeny`
- `ariadne compare-fasta`

## 📎 Citation / Positioning

If you position Ariadne as a paper-style platform, it is best understood as:

> 🪸 A coral TPS / CeSS-centered workflow for discovery, prioritization, and phylogeny that is especially friendly for experimental screening.

## 🎨 Logo Assets

- `fig/ariadne_cvpr_logo.svg`
  Horizontal banner logo for README headers, slides, and project pages
- `fig/ariadne_icon_minimal.svg`
  Minimal icon for paper figure corners, GitHub avatars, and favicon-style usage

## License

MIT License
