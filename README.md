<p align="center">
  <img src="fig/ariadne_cvpr_logo.svg" alt="Ariadne logo" width="980">
</p>

<p align="center">
  <a href="./README_ZH.md"><strong>中文</strong></a>
</p>

<p align="center">
  🧬 <strong>Ariadne</strong><br>
  A tree-native terpene synthase discovery and phylogeny platform for coral-centered and cross-clade TPS mining
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/Python-3.11%2B-0F172A?style=flat-square&logo=python&logoColor=white">
  <img alt="Workflow" src="https://img.shields.io/badge/Workflow-4%20Stages-0F766E?style=flat-square">
  <img alt="References" src="https://img.shields.io/badge/Default%20Reference-tree%2F-C2410C?style=flat-square">
  <img alt="MAFFT" src="https://img.shields.io/badge/MAFFT-required-14532D?style=flat-square">
  <img alt="IQ-TREE" src="https://img.shields.io/badge/IQ--TREE-required-1D4ED8?style=flat-square">
  <a href="./docs/index.md"><img alt="Documentation" src="https://img.shields.io/badge/Documentation-English%20Docs-2563EB?style=flat-square&logo=readthedocs&logoColor=white"></a>
</p>

<p align="center">
  <a href="./docs/index.md">
    <img alt="Open Documentation" src="https://img.shields.io/badge/Open-Documentation-0B132B?style=for-the-badge&logo=readthedocs&logoColor=white">
  </a>
</p>

> Ariadne is a research-oriented command-line platform for terpene synthase discovery, feature-space classification, multiple sequence alignment, and phylogenetic reconstruction.  
> The current release uses a streamlined `discovery -> filtering -> classification -> phylogeny` workflow with `tree/` as the default reference source.

## 📢 News

- `2026-03-21`: GitHub landing page switched to English by default, with a dedicated Chinese documentation page at [README_ZH.md](./README_ZH.md).
- `2026-03-21`: The pipeline was simplified to a cleaner four-stage workflow. All `motif` and `benchmark` functionality was removed.
- `2026-03-21`: `tree/` is now the default source for query-HMM construction, TPS HMM library generation, classification background, and phylogeny references.

## ✨ Highlights

- 🌊 Tree-native workflow: one curated `tree/` directory drives discovery, classification, and phylogeny.
- 🧭 Feature-space screening: candidates are projected into a TPS HMM embedding for nearest-reference assignment and visual inspection.
- 🌳 Phylogeny-ready by default: after classification, Ariadne directly builds a MAFFT alignment and IQ-TREE phylogeny.
- 🧪 Practical CLI design: full end-to-end `run`, plus modular `discover`, `filter`, `classify`, and `phylogeny` commands.
- 📄 Figure-friendly outputs: the pipeline exports `embedding.svg`, `embedding_3d_sections.svg`, local context trees, and final Newick trees.

## 🧠 Overview

The current Ariadne workflow is intentionally compact:

1. `discovery`
   Build a query HMM from the reference FASTA files in `tree/`, then search protein FASTA files or transcriptome-derived ORFs for candidate TPS sequences.
2. `filtering`
   Remove low-coverage, too-short, and near-duplicate candidates.
3. `classification`
   Score all references and candidates against a TPS HMM library, embed them in feature space, and assign each candidate to its nearest reference source.
4. `phylogeny`
   Merge filtered candidates with the reference collection, run MAFFT, and reconstruct a phylogeny with IQ-TREE.

## 🗂️ Repository Layout

```text
Ariadne/
├── ariadne/                # core package
├── input/                  # example protein inputs
├── tree/                   # default reference FASTA collection
├── output/                 # historical example outputs, kept only for reference
├── fig/                    # logos and figures
├── environment.yml         # conda environment
├── install.sh
└── pyproject.toml
```

## 📁 Data Convention

- `input/`
  Example input folder for candidate mining. This is the standard source for `--protein-folder`.
- `tree/`
  The primary reference directory. Ariadne reads all multi-clade TPS FASTA files here and uses them to build the query HMM, TPS HMM library, classification background, and phylogeny background.
- `output/`
  A legacy example-output folder. It is no longer part of the default pipeline logic.

## 🛠️ Installation

### Conda

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

### venv

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
bash install.sh
```

### Manual

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## 📦 Dependencies

- Python `>= 3.9`
- `mafft`
- `iqtree` or `iqtree2`
- `numpy >= 1.24`
- `openpyxl >= 3.1`
- `pyhmmer >= 0.12.0`
- `pyrodigal >= 3.7.0`
- `scikit-learn >= 1.4`

Recommended environment: Python `3.11` via [environment.yml](./environment.yml).

## 🚀 Quick Start

### End-to-end run

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

This command will automatically:

- build the discovery query HMM from `tree/`
- build the TPS HMM library from `tree/`
- discover candidate proteins
- filter and deduplicate candidates
- classify them in TPS feature space
- generate the final alignment and phylogeny

### Transcriptome mode

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results_from_transcriptomes/
```

### Classification only

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_classification/
```

### Alignment and phylogeny only

```bash
ariadne phylogeny \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_phylogeny/
```

## 🎓 Tutorial

### Tutorial 1. Run the repository example

Use the included `input/` and `tree/` folders:

```bash
.venv/bin/python -m ariadne run \
  --protein-folder input \
  --reference-dir tree \
  --output-dir tmp_run_example
```

Expected outputs:

- `tmp_run_example/01_discovery/`
- `tmp_run_example/02_filtering/`
- `tmp_run_example/03_classification/`
- `tmp_run_example/04_phylogeny/`
- `tmp_run_example/pipeline_summary.tsv`

### Tutorial 2. Inspect classification outputs

Start with:

- `03_classification/classification.tsv`
- `03_classification/nearest_neighbors.tsv`
- `03_classification/embedding.svg`
- `03_classification/embedding_3d_sections.svg`

These files tell you which reference source each candidate is closest to, how confident that assignment is, and where the candidate lies in TPS feature space.

### Tutorial 3. Inspect phylogeny outputs

Start with:

- `04_phylogeny/phylogeny_input.fasta`
- `04_phylogeny/phylogeny_alignment.fasta`
- `04_phylogeny/iqtree.treefile`
- `04_phylogeny/iqtree.iqtree`

These files give you the exact alignment and tree used for downstream interpretation, figure generation, or manual curation.

## 🧪 CLI Reference

### Main commands

| Command | Purpose |
| --- | --- |
| `ariadne run` | full end-to-end workflow |
| `ariadne discover` | HMM-based candidate discovery from proteins or transcriptomes |
| `ariadne filter` | coverage, length, and near-duplicate filtering |
| `ariadne classify` | TPS feature-space classification |
| `ariadne phylogeny` | MAFFT alignment plus IQ-TREE reconstruction |
| `ariadne build-hmm` | build a single query HMM from a FASTA/MSA source |
| `ariadne build-tps-hmm-library` | build a TPS HMM library directory from reference FASTA/MSA files |

### `ariadne run`

| Parameter | Default | Description |
| --- | --- | --- |
| `--protein-folder` | `None` | folder of protein FASTA files |
| `--transcriptomes` | `None` | transcriptome FASTA files for ORF prediction mode |
| `--reference-dir` | required | tree-native reference directory |
| `--output-dir` | required | output root |
| `--query-hmm` | auto | use a prebuilt discovery HMM instead of auto-building from `tree/` |
| `--tps-hmm-dir` | auto | use a prebuilt TPS HMM library instead of auto-building from `tree/` |
| `--min-coverage` | `10.0` | filtering coverage threshold |
| `--min-length` | `300` | filtering minimum amino-acid length |
| `--identity-threshold` | `0.95` | near-duplicate collapse threshold |
| `--top-k` | `5` | nearest-reference voting size |
| `--tree-neighbors` | `12` | reference neighbors per local context tree |
| `--skip-phylogeny` | `False` | skip MAFFT + IQ-TREE |
| `--mafft-mode` | `--auto` | MAFFT mode |
| `--iqtree-model` | `LG` | IQ-TREE model |
| `--iqtree-threads` | `AUTO` | IQ-TREE threads |
| `--iqtree-bootstrap` | `None` | optional ultrafast bootstrap replicates |
| `--no-iqtree-fast` | `False` | disable default `--fast` mode |

### `ariadne discover`

| Parameter | Default | Description |
| --- | --- | --- |
| `--protein-folder` | `None` | folder of predicted protein FASTA files |
| `--transcriptomes` | `None` | transcript FASTA inputs |
| `--protein-glob` | common FASTA patterns | file globs under `--protein-folder` |
| `--hmm` | required | discovery HMM |
| `--output-dir` | required | discovery output directory |
| `--min-score` | `None` | optional minimum HMM score |
| `--max-evalue` | `None` | optional maximum E-value |

### `ariadne filter`

| Parameter | Default | Description |
| --- | --- | --- |
| `--input-fasta` | required | candidate protein FASTA |
| `--output-dir` | required | filter output directory |
| `--min-coverage` | `10.0` | remove low-coverage records |
| `--min-length` | `300` | remove short proteins |
| `--identity-threshold` | `0.95` | near-duplicate threshold |

### `ariadne classify`

| Parameter | Default | Description |
| --- | --- | --- |
| `--candidates` | required | filtered candidate FASTA |
| `--reference-dir` | required | reference FASTA directory |
| `--output-dir` | required | classification output directory |
| `--tps-hmm-dir` | auto | TPS HMM library directory |
| `--top-k` | `5` | nearest-reference voting size |
| `--tree-neighbors` | `12` | neighbors used for local context trees |

### `ariadne phylogeny`

| Parameter | Default | Description |
| --- | --- | --- |
| `--candidates` | required | filtered candidate FASTA |
| `--reference-dir` | required | reference FASTA directory |
| `--output-dir` | required | phylogeny output directory |
| `--mafft-bin` | auto | explicit MAFFT binary |
| `--mafft-mode` | `--auto` | MAFFT mode |
| `--iqtree-bin` | auto | explicit IQ-TREE binary |
| `--iqtree-model` | `LG` | substitution model |
| `--iqtree-threads` | `AUTO` | thread setting |
| `--iqtree-bootstrap` | `None` | optional bootstrap |
| `--no-iqtree-fast` | `False` | disable default fast mode |

## 📤 Output Structure

```text
results/
├── 01_discovery/
│   ├── candidates.protein.faa
│   ├── candidates.orf.fna
│   └── candidates.hits.tsv
├── 02_filtering/
│   ├── candidates.filtered.faa
│   ├── filter_report.tsv
│   ├── dedupe_clusters.tsv
│   └── manual_review.tsv
├── 03_classification/
│   ├── tps_features.tsv
│   ├── embedding.tsv
│   ├── embedding.svg
│   ├── embedding_3d_sections.svg
│   ├── classification.tsv
│   ├── nearest_neighbors.tsv
│   ├── candidate_cluster_context.tsv
│   ├── global_context_tree.nwk
│   └── trees/
├── 04_phylogeny/
│   ├── phylogeny_input.fasta
│   ├── phylogeny_alignment.fasta
│   ├── phylogeny_sequence_map.tsv
│   ├── iqtree.treefile
│   └── iqtree.iqtree
└── pipeline_summary.tsv
```

## 🔬 Recommended Reading Order for Results

1. Open `pipeline_summary.tsv` to verify all stages completed.
2. Read `03_classification/classification.tsv` for per-candidate assignments.
3. Inspect `03_classification/embedding.svg` for global candidate placement.
4. Read `04_phylogeny/iqtree.treefile` and `04_phylogeny/iqtree.iqtree` for phylogenetic interpretation.

## ⚠️ Notes

- `tree/` is now the canonical source of reference FASTA files.
- `Alignment.fasta`-based entrypoints are no longer part of the main workflow.
- `motif` and `benchmark` functionality were intentionally removed to keep the current release focused and stable.
- The repository still contains historical example outputs under `output/`, but they are not required for the current pipeline.

## 📚 Citation

If Ariadne is useful in your work, please cite the software entry below and update it with your preferred version tag or DOI when available.

```bibtex
@software{jiang2026ariadne,
  author       = {Jiang, Zhaorui},
  title        = {Ariadne: A Tree-Native Terpene Synthase Discovery and Phylogeny Platform},
  year         = {2026},
  url          = {https://github.com/zhaoruijiang26/Ariadne},
  version      = {0.1.0}
}
```

## 🤝 Acknowledgement

Ariadne is designed as a practical bridge between TPS candidate mining and downstream phylogenetic interpretation, especially for coral-centered discovery settings where curated cross-clade references matter.
