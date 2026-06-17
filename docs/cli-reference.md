# CLI Reference

Ariadne exposes one command for each main workflow stage plus `run`, which executes the README-aligned workflow through classification and PCA/LDA visualization. Required parameters are marked **required**.

## Command Map

| Command | Purpose |
|---|---|
| `ariadne run` | end-to-end workflow through classification and visualization |
| `ariadne filter` | Stage 2: quality control and near-duplicate collapsing |
| `ariadne classify` | Stages 3-4: classification plus PCA/LDA visualization |
| `ariadne prepare-references` | prepare reference FASTA files from source data |

The installed CLI may still expose legacy or developer helper commands, including lower-level discovery helpers. They are not part of the current no-prep README workflow; use `ariadne run` for Stage 1.

---

## `ariadne run`

End-to-end workflow through discovery, filtering, classification, and visualization. Requires either `--protein-folder` or `--transcriptomes`.

### Input

| Parameter | Default | Description |
|---|---|---|
| `--protein-folder PATH` | `None` | directory of protein FASTA files |
| `--transcriptomes PATH [PATH ...]` | `None` | one or more transcriptome FASTA files; ORFs are predicted with Pyrodigal |
| `--protein-glob GLOB [GLOB ...]` | common FASTA extensions | override recursive glob patterns under `--protein-folder` |
| `--reference-dir PATH` | **required** | reference FASTA directory |
| `--output-dir PATH` | **required** | root output directory |

### Stage 1 - Discovery

| Parameter | Default | Description |
|---|---|---|
| `--discovery-min-score FLOAT` | `None` | minimum discovery score |
| `--discovery-max-evalue FLOAT` | `None` | maximum E-value |

### Stage 2 - Filtering

| Parameter | Default | Description |
|---|---|---|
| `--min-coverage FLOAT` | `10.0` | minimum per-base sequencing coverage parsed from the FASTA header |
| `--min-length INT` | `300` | minimum protein length in amino acids |
| `--identity-threshold FLOAT` | `0.95` | near-duplicate collapsing threshold and reference-match logging threshold |

### Stages 3-4 - Classification And Visualization

| Parameter | Default | Description |
|---|---|---|
| `--top-k INT` | `5` | number of nearest reference neighbors used for label voting |
| `--tree-neighbors INT` | `12` | number of nearest references used for local candidate context summaries |
| `--skip-phylogeny` | `False` | keep `run` aligned with the current README workflow by skipping the legacy tree step |

### CeeSs Scoring (Optional)

| Parameter | Default | Description |
|---|---|---|
| `--ceess-xlsx PATH` | `TPS.xlsx` when supplied | labeled coral TPS workbook (`Name` / `Protein` / `Type` / `Species`) |
| `--skip-ceess-model` | `False` | skip the ESM-based CeeSs scoring stage |
| `--ceess-threshold FLOAT` | `0.9` | minimum aggregated `P(CeeSs)` for `ceess_candidates.tsv` |
| `--ceess-classifier {mlp,logreg,contrastive}` | `mlp` | classifier head on top of frozen ESM2 embeddings |
| `--ceess-model-name NAME` | `facebook/esm2_t33_650M_UR50D` | ESM2 preset name or Hugging Face model ID |
| `--ceess-batch-size INT` | `4` | batch size for ESM2 inference |
| `--ceess-max-length INT` | `2048` | maximum tokenized sequence length |
| `--ceess-device DEVICE` | auto | explicit torch device, for example `cuda:0` or `cpu` |
| `--ceess-cv-folds INT` | `5` | target cross-validation folds |
| `--ceess-random-state INT` | `0` | random seed |
| `--ceess-epochs INT` | `200` | MLP training epochs |
| `--ceess-hidden-dim INT` | `128` | MLP hidden layer width |
| `--ceess-dropout FLOAT` | `0.1` | MLP dropout rate |
| `--ceess-learning-rate FLOAT` | `1e-3` | MLP AdamW learning rate |
| `--ceess-weight-decay FLOAT` | `1e-4` | MLP AdamW weight decay |
| `--ceess-train-batch-size INT` | `8` | MLP training mini-batch size |
| `--ceess-barlow-representation-dim INT` | `None` | encoder output width for Barlow Twins (`contrastive` only) |
| `--ceess-barlow-projection-dim INT` | `None` | projection head output width for Barlow Twins (`contrastive` only) |
| `--ceess-barlow-redundancy-weight FLOAT` | `0.005` | off-diagonal redundancy penalty (`contrastive` only) |
| `--ceess-mlp-checkpoint PATH` | `None` | pretrained MLP `.pt` checkpoint; skips training |

---

## `ariadne filter`

Stage 2 standalone. Requires `--input-fasta` and `--output-dir`.

| Parameter | Default | Description |
|---|---|---|
| `--input-fasta PATH` | **required** | candidate protein FASTA from stage 1 |
| `--output-dir PATH` | **required** | filter output directory |
| `--min-coverage FLOAT` | `10.0` | remove candidates below this sequencing coverage |
| `--min-length INT` | `300` | remove proteins shorter than this length |
| `--identity-threshold FLOAT` | `0.95` | near-duplicate collapsing threshold |
| `--reference-dir PATH` | `None` | optional reference FASTA directory; reference matches are logged but kept |

---

## `ariadne classify`

Stages 3-4 standalone. Requires `--candidates`, `--reference-dir`, and `--output-dir`.

### Core

| Parameter | Default | Description |
|---|---|---|
| `--candidates PATH` | **required** | filtered candidate FASTA |
| `--reference-dir PATH` | **required** | reference FASTA directory |
| `--output-dir PATH` | **required** | classification and visualization output directory |
| `--top-k INT` | `5` | nearest-reference voting neighbors |
| `--tree-neighbors INT` | `12` | local context neighbor count |

### CeeSs Scoring (Optional)

The `classify` command accepts the same `--ceess-*` parameters as `run`.

---

## `ariadne prepare-references`

Prepare reference FASTA files from coral, insect, fungal, plant, bacterial, or extra sources.

| Parameter | Default | Description |
|---|---|---|
| `--coral PATH` | `None` | coral reference FASTA |
| `--coral-limit INT` | `None` | maximum number of coral reference sequences to include |
| `--insect-xlsx PATH` | `None` | insect TPS Excel workbook |
| `--insect-limit INT` | `None` | maximum number of insect sequences to include |
| `--bacteria-fasta PATH` | `None` | bacterial TPS FASTA |
| `--fungal-fasta PATH` | `None` | fungal TPS FASTA |
| `--fungi-fasta PATH` | `None` | alias of `--fungal-fasta` |
| `--plant-fasta PATH` | `None` | plant TPS FASTA |
| `--extra-fasta PATH [PATH ...]` | `None` | additional FASTA files to include |
| `--output-dir PATH` | **required** | output directory for prepared reference files |

---

## Practical Notes

- Use `--reference-dir reference_fastas/` or another directory containing your reference FASTAs.
- No separate HMM generation step is required for the documented workflow.
- Global flags (`--verbose`, `--log-file`) must be placed before the subcommand: `ariadne --verbose run ...`
- Use `--skip-phylogeny` with `ariadne run` to keep outputs focused on `03_classification/`.
- The CeeSs head requires optional ESM dependencies. Install with `pip install 'ariadne-tps[esm]'` or `pip install torch transformers`.
- `--ceess-classifier contrastive` activates the Barlow Twins variant.
- `--ceess-mlp-checkpoint` loads a previously saved `.pt` file and skips training.
