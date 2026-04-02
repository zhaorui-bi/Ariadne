<p align="center">
  <img src="fig/ariadne_cvpr_logo.svg" alt="Ariadne logo" width="980">
</p>

<p align="center">
  <a href="./README_ZH.md"><strong>中文</strong></a>
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/Python-%3E%3D3.9-0F172A?style=flat-square&logo=python&logoColor=white">
  <img alt="Version" src="https://img.shields.io/badge/Version-1.0.0-0F766E?style=flat-square">
  <img alt="License" src="https://img.shields.io/badge/License-MIT-C2410C?style=flat-square">
  <img alt="ESM2" src="https://img.shields.io/badge/Optional-ESM2%20CeeSs-7C3AED?style=flat-square">
  <img alt="MAFFT" src="https://img.shields.io/badge/Requires-MAFFT%20%2B%20IQ--TREE-14532D?style=flat-square">
  <a href="./docs/index.md"><img alt="Docs" src="https://img.shields.io/badge/Docs-Read%20the%20Docs-2563EB?style=flat-square&logo=readthedocs&logoColor=white"></a>
</p>

---

## Abstract

Coral terpene synthases (TPSs) represent an underexplored frontier in natural product biosynthesis. Identifying which of the hundreds of predicted coral TPS proteins is responsible for a specific terpenoid product — in particular the cembrene-class (CeeSs) compounds — requires more than sequence homology: it demands systematic embedding in a curated reference space.

**Ariadne** is a tree-native, four-stage computational platform for coral TPS mining and CeeSs prioritization. Starting from raw transcriptomes or predicted proteomes, it combines profile HMM-guided discovery, coverage- and length-aware filtering, TPS feature-space embedding with supervised dimensionality reduction, and an optional ESM2-based scoring layer for CeeSs candidates. The same `tree/` reference backbone drives all four stages, ensuring biological consistency from candidate discovery through final phylogenetic placement.

---

## Highlights

- **Tree-native design** — a single `tree/` reference directory drives discovery, classification, and phylogeny without manual bookkeeping between stages.
- **HMM feature space** — candidates are scored against a multi-source TPS HMM library, embedded with LDA/PCA, and assigned nearest reference neighbors before any phylogenetic tree is built.
- **ESM2 CeeSs scoring** — when `TPS/TPS.xlsx` and the optional ESM stack are available, a frozen ESM2 backbone with a trainable MLP head scores each coral-like candidate for cembrene A / cembrene B probability. A Barlow Twins contrastive variant is also supported.
- **Publication-quality outputs** — the pipeline emits SVG embedding figures, per-candidate UPGMA trees, a MAFFT alignment, an IQ-TREE phylogeny, and a complete TSV audit trail.

---

## Method Overview

<p align="center">
  <img src="./docs/assets/overview_pipeline.svg" alt="Ariadne pipeline overview" width="100%">
</p>

Ariadne runs as a sequential four-stage pipeline. Each stage can also be invoked independently.

| Stage | Command | Input | Key Output |
|---|---|---|---|
| 1. Discovery | `ariadne discover` | Transcriptomes / protein FASTAs | `candidates.protein.faa` |
| 2. Filtering | `ariadne filter` | Candidate FASTA | `candidates.filtered.faa` |
| 3. Classification | `ariadne classify` | Filtered FASTA + `tree/` | `classification.tsv`, `embedding.svg` |
| 4. Phylogeny | `ariadne phylogeny` | Filtered FASTA + `tree/` | `iqtree.treefile`, `phylogeny_preview.svg` |

**Stage 1 — Discovery.**  
ORFs are predicted from transcriptome assemblies with Pyrodigal (meta mode) and translated proteins are searched with a query profile HMM built from the coral reference alignment. When protein FASTAs are provided directly, the ORF step is skipped.

**Stage 2 — Filtering.**  
Candidates are filtered by coverage (`≥ 10×` default), minimum length (`≥ 300 aa`), and near-duplicate collapsing at 95% identity using a bounded edit-distance algorithm. Optionally, sequences already present in the reference set are removed.

**Stage 3 — Classification.**  
All sequences (references + candidates) are scored against the TPS HMM library, producing a per-sequence feature vector. The matrix is z-scored and reduced to 3D via LDA (supervised, with KMeans subclustering of the large coral reference set) or PCA as fallback. Candidates are assigned their nearest-reference labels by k-NN voting. An optional ESM2 CeeSs sub-stage runs a frozen ESM2 backbone + MLP head trained on `TPS.xlsx` labeled sequences, producing per-candidate `P(CeeSs)` scores.

**Stage 4 — Phylogeny.**  
References and filtered candidates are merged into a deduplicated FASTA, aligned with MAFFT, and a maximum-likelihood tree is inferred with IQ-TREE. A compact SVG preview is rendered from the resulting Newick file.

---

## Results Preview

<p align="center">
  <img src="./docs/assets/latest_embedding.svg" alt="TPS feature-space embedding" width="100%">
</p>

<p align="center">
  <em>TPS HMM feature-space embedding from a representative run. 100 candidates discovered → 36 retained after filtering → 36 classified as coral-like → 5 CeeSs candidates shortlisted (P(CeeSs) ≥ 0.9).</em>
</p>

<p align="center">
  <img src="./docs/assets/latest_embedding_3d_sections.svg" alt="3D embedding sections" width="100%">
</p>

---

## Installation

**Recommended** — Python 3.11 with the bundled Conda environment:

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

**Minimal** — core dependencies only:

```bash
pip install -e .
# Runtime: mafft and iqtree (or iqtree2) must be on PATH
```

**With ESM2 CeeSs scoring:**

```bash
pip install -e '.[esm]'
# Requires: torch, transformers, tqdm
```

Core Python dependencies: `numpy >= 1.24`, `pyhmmer >= 0.12.0`, `pyrodigal >= 3.7.0`, `scikit-learn >= 1.4`, `openpyxl >= 3.1`

---

## Quick Start

### Full pipeline from protein FASTAs

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

Output layout:

```
results/
├── 01_discovery/          # HMM hits, per-sample protein FASTAs
├── 02_filtering/          # filtered FASTA, filter_report.tsv, dedupe_clusters.tsv
├── 03_classification/     # classification.tsv, embedding.svg, per-candidate trees
│   └── (ceess_*/          # optional CeeSs outputs when TPS/TPS.xlsx is present)
├── 04_phylogeny/          # iqtree.treefile, phylogeny_preview.svg
└── pipeline_summary.tsv
```

### From transcriptome assemblies

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results/
```

### Classification only

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results/03_classification/
```

### Phylogeny only

```bash
ariadne phylogeny \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results/04_phylogeny/
```

---

## CeeSs Scoring (Optional)

When `TPS/TPS.xlsx` is present and the ESM stack is installed, stages `run` and `classify` automatically attach an ESM2 CeeSs scoring pass after the HMM classification step.

**Scoring pipeline:**

1. Load labeled TPS sequences from `TPS.xlsx` (Name / Protein / Type / Species columns).
2. Compute frozen ESM2 mean-pooled embeddings for all training sequences and coral-like candidates.
3. Train a small MLP head (default) or logistic regression on the training embeddings.
4. Score each coral-like candidate; report `P(CeeSs)` as the summed probability over all CeeSs-positive type labels.
5. Candidates above `--ceess-threshold` (default 0.9) are written to `ceess_candidates.tsv` and `ceess_candidates.fasta`.

**Classifier options** (`--ceess-classifier`):

| Option | Description |
|---|---|
| `mlp` (default) | Torch MLP with cross-entropy loss, AdamW, class weighting |
| `logreg` | Sklearn logistic regression with standard scaler |
| `contrastive` | Barlow Twins projection network + MLP head |

**Key CeeSs outputs:**

| File | Description |
|---|---|
| `ceess_predictions.tsv` | Per-candidate scores: `esm_type_prediction`, `esm_ceess_probability`, one probability column per TPS type |
| `ceess_candidates.tsv` | Shortlisted candidates above threshold |
| `ceess_candidates.fasta` | FASTA of shortlisted candidates |
| `ceess_embedding.svg` | 2D LDA/PCA projection of training references and candidates |
| `ceess_model_metrics.tsv` | Cross-validated accuracy, F1, confusion matrix |

---

## Module Layout

```
ariadne/
├── utils.py       # logging, terminal output, FASTA I/O, sequence utilities
├── data.py        # reference data management (coral, insect, plant, fungi, bacteria)
├── search.py      # Stage 1: HMM construction and candidate discovery
├── filter.py      # Stage 2: coverage, length, and near-duplicate filtering
├── embed.py       # Stage 3: HMM feature matrix, embedding, classification
├── model.py       # ESM2 CeeSs scoring (MLP, logistic regression, Barlow Twins)
├── tree.py        # Stage 4: MAFFT alignment, IQ-TREE phylogeny, SVG preview
├── demo.py        # demo workspace generator
├── cli.py         # command-line interface
├── __init__.py
├── __main__.py
└── hmm/           # bundled query HMM and TPS HMM library
```

---

## Repository Layout

```
Ariadne/
├── ariadne/           # core package
├── docs/              # project documentation
├── fig/               # figures and logos
├── input/             # example protein inputs
├── TPS/               # optional labeled coral TPS workbook (TPS.xlsx)
├── tree/              # default reference FASTA collection
├── environment.yml
└── pyproject.toml
```

---

## CLI Reference

```
ariadne --help

Commands:
  prepare-references    Prepare reference FASTA files from coral / insect / extra sources
  prepare-demo          Create a small reproducible demo workspace
  build-hmm             Build a profile HMM from a reference alignment
  build-tps-hmm-library Build a TPS HMM library from multiple reference FASTAs
  discover              Stage 1: HMM-guided candidate discovery
  filter                Stage 2: quality filtering and deduplication
  classify              Stage 3: feature-space embedding and classification
  phylogeny             Stage 4: MAFFT alignment + IQ-TREE phylogeny
  run                   Full pipeline (stages 1–4)
```

Full CLI documentation: [docs/cli-reference.md](./docs/cli-reference.md)

---

## Documentation

- [Getting Started](./docs/index.md)
- [CLI Reference](./docs/cli-reference.md)
- [Output Files](./docs/outputs.md)
- [CeeSs / ESM Type](./docs/esm-type.md)
- [Citation](./docs/citation.md)

---

## Citation

```bibtex
@software{jiang2026ariadne,
  author    = {Jiang, Zhaorui},
  title     = {Ariadne: A Coral-Centered Terpene Synthase Discovery and CeeSs Prioritization Platform},
  year      = {2026},
  url       = {https://github.com/zhaoruijiang26/Ariadne},
  version   = {1.0.0}
}
```
