<p align="center">
  <img src="fig/logo.png" alt="Ariadne logo" width="980">
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

## 🪸 Abstract

Coral terpene synthases (TPSs) represent an underexplored frontier in natural product biosynthesis. Identifying which of the hundreds of predicted coral TPS proteins is responsible for a specific terpenoid product — in particular the cembrene-class (CeeSs) compounds — requires more than sequence homology: it demands systematic embedding in a curated reference space.

**Ariadne** is a tree-native, four-stage computational platform for coral TPS mining and CeeSs prioritization. Starting from raw transcriptomes or predicted proteomes, it combines profile HMM-guided discovery, coverage- and length-aware filtering, TPS feature-space embedding with supervised dimensionality reduction, and an optional ESM2-based scoring layer for CeeSs candidates. The same `tree/` reference backbone drives all four stages, ensuring biological consistency from candidate discovery through final phylogenetic placement.

---

## ✨ Contributions

- **Tree-native design** — a single `tree/` reference directory drives discovery, classification, and phylogeny without manual bookkeeping between stages.
- **HMM feature space** — candidates are scored against a multi-source TPS HMM library, embedded with LDA/PCA, and assigned nearest reference neighbors before any phylogenetic tree is built.
- **ESM2 CeeSs scoring** — when `TPS/TPS.xlsx` and the optional ESM stack are available, a frozen ESM2 backbone with a trainable MLP head scores each coral-like candidate for cembrene A / cembrene B probability. A Barlow Twins contrastive variant is also supported.
- **Publication-quality outputs** — the pipeline emits SVG embedding figures, per-candidate UPGMA trees, a MAFFT alignment, an IQ-TREE phylogeny, and a complete TSV audit trail.

---

## 🔬 Method

<p align="center">
  <img src="fig/framework.png" alt="Ariadne pipeline overview" width="100%">
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

## ⚙️ Installation

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

## 🚀 Quick Start

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

## 🧬 CeeSs Scoring (Optional)

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

## 📦 Module Layout

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

## 🗂️ Repository Layout

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

## ⌨️ CLI Reference

| Command | Purpose |
|---|---|
| `ariadne run` | Full end-to-end pipeline (stages 1–4) |
| `ariadne discover` | Stage 1: HMM-guided candidate discovery |
| `ariadne filter` | Stage 2: quality filtering and deduplication |
| `ariadne classify` | Stage 3: feature-space embedding and classification |
| `ariadne phylogeny` | Stage 4: MAFFT alignment + IQ-TREE phylogeny |
| `ariadne prepare-references` | Prepare reference FASTA files from source data |
| `ariadne build-hmm` | Build a single profile HMM from an alignment |
| `ariadne build-tps-hmm-library` | Build a TPS HMM library from multiple alignments |

### `ariadne run` — full pipeline

| Parameter | Default | Description |
|---|---|---|
| `--protein-folder PATH` | `None` | directory of protein FASTA files |
| `--transcriptomes PATH …` | `None` | transcriptome FASTAs; ORFs predicted with Pyrodigal |
| `--protein-glob GLOB …` | auto | override recursive glob patterns under `--protein-folder` |
| `--query-hmm PATH` | bundled | profile HMM for discovery; auto-built from `tree/` as fallback |
| `--reference-dir PATH` | **required** | tree-native reference FASTA directory |
| `--output-dir PATH` | **required** | output root directory |
| `--hmm-name NAME` | `ariadne_query` | name for any auto-built discovery HMM |
| `--discovery-min-score FLOAT` | `None` | minimum HMM bitscore for discovery |
| `--discovery-max-evalue FLOAT` | `None` | maximum E-value for discovery |
| `--min-coverage FLOAT` | `10.0` | minimum sequencing coverage (filtering) |
| `--min-length INT` | `300` | minimum protein length in aa (filtering) |
| `--identity-threshold FLOAT` | `0.95` | near-duplicate collapsing threshold |
| `--tps-hmm-dir PATH` | auto | TPS HMM library; uses bundled or auto-builds from `tree/` |
| `--top-k INT` | `5` | nearest-reference voting size |
| `--tree-neighbors INT` | `12` | neighbors for local candidate context trees |
| `--ceess-xlsx PATH` | `TPS/TPS.xlsx` | labeled coral TPS workbook for ESM CeeSs scoring |
| `--skip-ceess-model` | `False` | skip ESM-based CeeSs scoring |
| `--ceess-threshold FLOAT` | `0.9` | minimum P(CeeSs) for `ceess_candidates.tsv` |
| `--ceess-classifier` | `mlp` | classifier head: `mlp`, `logreg`, or `contrastive` |
| `--ceess-model-name NAME` | `facebook/esm2_t33_650M_UR50D` | ESM2 preset or Hugging Face model ID |
| `--ceess-batch-size INT` | `4` | ESM2 inference batch size |
| `--ceess-max-length INT` | `2048` | maximum tokenized length for ESM2 |
| `--ceess-device DEVICE` | auto | torch device, e.g. `cuda:0` or `cpu` |
| `--ceess-cv-folds INT` | `5` | cross-validation folds for classifier evaluation |
| `--ceess-random-state INT` | `0` | random seed |
| `--ceess-epochs INT` | `200` | MLP training epochs |
| `--ceess-hidden-dim INT` | `128` | MLP hidden layer width |
| `--ceess-dropout FLOAT` | `0.1` | MLP dropout rate |
| `--ceess-learning-rate FLOAT` | `1e-3` | MLP learning rate |
| `--ceess-weight-decay FLOAT` | `1e-4` | MLP weight decay |
| `--ceess-train-batch-size INT` | `8` | MLP training batch size |
| `--ceess-barlow-representation-dim INT` | `None` | Barlow Twins encoder width (`contrastive` only) |
| `--ceess-barlow-projection-dim INT` | `None` | Barlow Twins projection width (`contrastive` only) |
| `--ceess-barlow-redundancy-weight FLOAT` | `0.005` | Barlow Twins off-diagonal penalty (`contrastive` only) |
| `--ceess-mlp-checkpoint PATH` | `None` | pretrained MLP `.pt` checkpoint; skips training |
| `--skip-phylogeny` | `False` | skip MAFFT + IQ-TREE |
| `--mafft-bin PATH` | auto | explicit MAFFT binary |
| `--mafft-mode FLAG` | `--auto` | MAFFT alignment mode |
| `--iqtree-bin PATH` | auto | explicit IQ-TREE binary |
| `--iqtree-model MODEL` | `LG` | IQ-TREE substitution model |
| `--iqtree-threads INT\|AUTO` | `AUTO` | IQ-TREE thread count |
| `--iqtree-bootstrap INT` | `None` | ultrafast bootstrap replicates |
| `--no-iqtree-fast` | `False` | disable IQ-TREE `--fast` mode |

### `ariadne discover`

| Parameter | Default | Description |
|---|---|---|
| `--protein-folder PATH` | `None` | protein FASTA directory |
| `--transcriptomes PATH …` | `None` | transcriptome FASTA inputs |
| `--protein-glob GLOB …` | auto | recursive search patterns |
| `--hmm PATH` | **required** | discovery HMM |
| `--output-dir PATH` | **required** | discovery output directory |
| `--min-score FLOAT` | `None` | minimum HMM bitscore |
| `--max-evalue FLOAT` | `None` | maximum E-value |

### `ariadne filter`

| Parameter | Default | Description |
|---|---|---|
| `--input-fasta PATH` | **required** | candidate protein FASTA |
| `--output-dir PATH` | **required** | filter output directory |
| `--min-coverage FLOAT` | `10.0` | minimum sequencing coverage |
| `--min-length INT` | `300` | minimum protein length (aa) |
| `--identity-threshold FLOAT` | `0.95` | near-duplicate threshold |
| `--reference-dir PATH` | `None` | reference directory; matches are logged in `reference_matches.tsv` and **kept** |

### `ariadne classify`

Accepts the same `--ceess-*` flags as `ariadne run` (see full table above).

| Parameter | Default | Description |
|---|---|---|
| `--candidates PATH` | **required** | filtered candidate FASTA |
| `--reference-dir PATH` | **required** | reference FASTA directory |
| `--output-dir PATH` | **required** | classification output directory |
| `--tps-hmm-dir PATH` | auto | TPS HMM library directory |
| `--top-k INT` | `5` | voting neighbors |
| `--tree-neighbors INT` | `12` | local context-tree neighbors |

### `ariadne phylogeny`

| Parameter | Default | Description |
|---|---|---|
| `--candidates PATH` | **required** | filtered candidate FASTA |
| `--reference-dir PATH` | **required** | reference FASTA directory |
| `--output-dir PATH` | **required** | phylogeny output directory |
| `--mafft-bin PATH` | auto | explicit MAFFT binary |
| `--mafft-mode FLAG` | `--auto` | MAFFT alignment mode |
| `--iqtree-bin PATH` | auto | explicit IQ-TREE binary |
| `--iqtree-model MODEL` | `LG` | substitution model |
| `--iqtree-threads INT\|AUTO` | `AUTO` | threads |
| `--iqtree-bootstrap INT` | `None` | optional bootstrap replicates |
| `--no-iqtree-fast` | `False` | disable fast mode |

### `ariadne prepare-references`

| Parameter | Default | Description |
|---|---|---|
| `--coral PATH` | `None` | coral reference FASTA |
| `--coral-limit INT` | `None` | maximum coral sequences |
| `--insect-xlsx PATH` | `None` | insect TPS Excel workbook |
| `--insect-limit INT` | `None` | maximum insect sequences |
| `--bacteria-fasta PATH` | `None` | bacterial TPS FASTA |
| `--fungal-fasta PATH` | `None` | fungal TPS FASTA |
| `--plant-fasta PATH` | `None` | plant TPS FASTA |
| `--extra-fasta PATH …` | `None` | additional FASTA files |
| `--output-dir PATH` | **required** | output directory |

### `ariadne build-hmm` / `ariadne build-tps-hmm-library`

| Command | Parameter | Description |
|---|---|---|
| `build-hmm` | `--alignment PATH` (required) | input alignment or FASTA |
| `build-hmm` | `--output PATH` (required) | output `.hmm` file |
| `build-hmm` | `--name NAME` | profile name |
| `build-tps-hmm-library` | `--alignment NAME=PATH …` (required) | named alignment files |
| `build-tps-hmm-library` | `--output-dir PATH` (required) | output HMM library directory |

> **Note:** Global flags (`--verbose`, `--log-file`) must precede the subcommand: `ariadne --verbose run ...`

Full CLI documentation: [docs/cli-reference.md](./docs/cli-reference.md)

---

## 📖 Documentation

- [Getting Started](./docs/index.md)
- [CLI Reference](./docs/cli-reference.md)
- [Output Files](./docs/outputs.md)
- [CeeSs / ESM Type](./docs/esm-type.md)
- [Citation](./docs/citation.md)

---

## 📄 Citation

```bibtex
@software{jiang2026ariadne,
  author    = {Jiang, Zhaorui},
  title     = {Ariadne: A Coral-Centered Terpene Synthase Discovery and CeeSs Prioritization Platform},
  year      = {2026},
  url       = {https://github.com/zhaoruijiang26/Ariadne},
  version   = {1.0.0}
}
```
