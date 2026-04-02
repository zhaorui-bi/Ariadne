# CLI Reference

## Command map

| Command | Purpose |
| --- | --- |
| `ariadne run` | full end-to-end workflow (stages 1–4) |
| `ariadne discover` | Stage 1: HMM-based candidate discovery |
| `ariadne filter` | Stage 2: quality control and near-duplicate collapsing |
| `ariadne classify` | Stage 3: TPS feature-space classification |
| `ariadne phylogeny` | Stage 4: MAFFT alignment and IQ-TREE reconstruction |
| `ariadne prepare-references` | prepare reference FASTA files from source data |
| `ariadne build-hmm` | build a single HMM from a FASTA/MSA source |
| `ariadne build-tps-hmm-library` | build a directory of TPS HMMs from references |

---

## `ariadne run`

Full end-to-end pipeline (stages 1–4). Requires either `--protein-folder` or `--transcriptomes`.

### Input

| Parameter | Default | Description |
| --- | --- | --- |
| `--protein-folder PATH` | `None` | directory of protein FASTA files; all `.faa`, `.fa`, `.fasta` files are discovered recursively |
| `--transcriptomes PATH [PATH ...]` | `None` | one or more transcriptome FASTA files; ORFs are predicted with Pyrodigal before HMM search |
| `--protein-glob GLOB [GLOB ...]` | common FASTA extensions | override the recursive glob patterns used under `--protein-folder` |
| `--query-hmm PATH` | bundled `ariadne/hmm/query.hmm` | explicit profile HMM for candidate discovery; auto-built from `tree/` if bundled HMM is absent |
| `--reference-dir PATH` | **required** | tree-native reference FASTA directory used by all four stages |
| `--output-dir PATH` | **required** | root output directory; stage subdirectories are created automatically |
| `--hmm-name NAME` | `ariadne_query` | name written into any auto-built discovery HMM |

### Stage 1 — Discovery

| Parameter | Default | Description |
| --- | --- | --- |
| `--discovery-min-score FLOAT` | `None` | minimum HMM bitscore; sequences below this value are excluded |
| `--discovery-max-evalue FLOAT` | `None` | maximum E-value; sequences above this value are excluded |

### Stage 2 — Filtering

| Parameter | Default | Description |
| --- | --- | --- |
| `--min-coverage FLOAT` | `10.0` | minimum per-base sequencing coverage; parsed from the FASTA header |
| `--min-length INT` | `300` | minimum protein length in amino acids |
| `--identity-threshold FLOAT` | `0.95` | near-duplicate collapsing threshold; also used for reference-match logging |

### Stage 3 — Classification

| Parameter | Default | Description |
| --- | --- | --- |
| `--tps-hmm-dir PATH` | auto | explicit TPS HMM library directory; uses bundled `ariadne/hmm/` or auto-builds from `tree/` if absent |
| `--top-k INT` | `5` | number of nearest reference neighbors used for label voting |
| `--tree-neighbors INT` | `12` | number of nearest neighbors used to build per-candidate local context trees |

### Stage 3 — CeeSs scoring (optional)

| Parameter | Default | Description |
| --- | --- | --- |
| `--ceess-xlsx PATH` | `TPS/TPS.xlsx` when present | labeled coral TPS workbook (Name / Protein / Type / Species columns) |
| `--skip-ceess-model` | `False` | skip the ESM-based CeeSs scoring stage entirely |
| `--ceess-threshold FLOAT` | `0.9` | minimum aggregated P(CeeSs) for a candidate to appear in `ceess_candidates.tsv` |
| `--ceess-classifier {mlp,logreg,contrastive}` | `mlp` | classifier head on top of frozen ESM2 embeddings |
| `--ceess-model-name NAME` | `facebook/esm2_t33_650M_UR50D` | ESM2 preset name or Hugging Face model ID |
| `--ceess-batch-size INT` | `4` | batch size for ESM2 inference |
| `--ceess-max-length INT` | `2048` | maximum tokenized sequence length for ESM2 |
| `--ceess-device DEVICE` | auto | explicit torch device, e.g. `cuda:0` or `cpu` |
| `--ceess-cv-folds INT` | `5` | target cross-validation folds for classifier evaluation on `TPS.xlsx` |
| `--ceess-random-state INT` | `0` | random seed for reproducible training |
| `--ceess-epochs INT` | `200` | MLP training epochs |
| `--ceess-hidden-dim INT` | `128` | MLP hidden layer width |
| `--ceess-dropout FLOAT` | `0.1` | MLP dropout rate |
| `--ceess-learning-rate FLOAT` | `1e-3` | MLP AdamW learning rate |
| `--ceess-weight-decay FLOAT` | `1e-4` | MLP AdamW weight decay |
| `--ceess-train-batch-size INT` | `8` | MLP training mini-batch size |
| `--ceess-barlow-representation-dim INT` | `None` | encoder output width for Barlow Twins (`contrastive` only) |
| `--ceess-barlow-projection-dim INT` | `None` | projection head output width for Barlow Twins (`contrastive` only) |
| `--ceess-barlow-redundancy-weight FLOAT` | `0.005` | off-diagonal redundancy penalty in Barlow Twins loss (`contrastive` only) |
| `--ceess-mlp-checkpoint PATH` | `None` | path to a pretrained MLP `.pt` checkpoint; skips training and loads directly |

### Stage 4 — Phylogeny

| Parameter | Default | Description |
| --- | --- | --- |
| `--skip-phylogeny` | `False` | skip MAFFT alignment and IQ-TREE inference |
| `--mafft-bin PATH` | auto | explicit path to the MAFFT binary |
| `--mafft-mode FLAG` | `--auto` | MAFFT alignment mode flag, e.g. `--linsi` |
| `--iqtree-bin PATH` | auto | explicit path to the IQ-TREE binary |
| `--iqtree-model MODEL` | `LG` | substitution model passed to IQ-TREE `-m` |
| `--iqtree-threads INT\|AUTO` | `AUTO` | number of threads for IQ-TREE |
| `--iqtree-bootstrap INT` | `None` | number of ultrafast bootstrap replicates; omit for no bootstrap |
| `--no-iqtree-fast` | `False` | disable the default `--fast` mode in IQ-TREE |

---

## `ariadne discover`

Stage 1 standalone. Requires either `--protein-folder` or `--transcriptomes`, plus `--hmm` and `--output-dir`.

| Parameter | Default | Description |
| --- | --- | --- |
| `--protein-folder PATH` | `None` | directory of protein FASTA files |
| `--transcriptomes PATH [PATH ...]` | `None` | transcriptome FASTA inputs; ORFs predicted with Pyrodigal |
| `--protein-glob GLOB [GLOB ...]` | common FASTA extensions | override recursive glob patterns under `--protein-folder` |
| `--hmm PATH` | **required** | profile HMM used for candidate discovery |
| `--output-dir PATH` | **required** | discovery output directory |
| `--min-score FLOAT` | `None` | minimum HMM bitscore cutoff |
| `--max-evalue FLOAT` | `None` | maximum E-value cutoff |

---

## `ariadne filter`

Stage 2 standalone. Requires `--input-fasta` and `--output-dir`.

| Parameter | Default | Description |
| --- | --- | --- |
| `--input-fasta PATH` | **required** | candidate protein FASTA from stage 1 |
| `--output-dir PATH` | **required** | filter output directory |
| `--min-coverage FLOAT` | `10.0` | remove candidates below this sequencing coverage |
| `--min-length INT` | `300` | remove proteins shorter than this length (aa) |
| `--identity-threshold FLOAT` | `0.95` | near-duplicate collapsing threshold |
| `--reference-dir PATH` | `None` | optional reference FASTA directory; candidates matching reference sequences at `--identity-threshold` are logged in `reference_matches.tsv` but are **kept** in the filtered output |

---

## `ariadne classify`

Stage 3 standalone. Requires `--candidates`, `--reference-dir`, and `--output-dir`.

### Core

| Parameter | Default | Description |
| --- | --- | --- |
| `--candidates PATH` | **required** | filtered candidate FASTA (output of stage 2) |
| `--reference-dir PATH` | **required** | reference FASTA directory |
| `--output-dir PATH` | **required** | classification output directory |
| `--tps-hmm-dir PATH` | auto | explicit TPS HMM library directory |
| `--top-k INT` | `5` | nearest-reference voting neighbors |
| `--tree-neighbors INT` | `12` | neighbors for local candidate context trees |

### CeeSs scoring (optional)

| Parameter | Default | Description |
| --- | --- | --- |
| `--ceess-xlsx PATH` | `TPS/TPS.xlsx` when present | labeled coral TPS workbook |
| `--skip-ceess-model` | `False` | skip ESM-based CeeSs scoring |
| `--ceess-threshold FLOAT` | `0.9` | minimum P(CeeSs) threshold |
| `--ceess-classifier {mlp,logreg,contrastive}` | `mlp` | classifier head |
| `--ceess-model-name NAME` | `facebook/esm2_t33_650M_UR50D` | ESM2 model name or Hugging Face ID |
| `--ceess-batch-size INT` | `4` | ESM2 inference batch size |
| `--ceess-max-length INT` | `2048` | maximum tokenized length for ESM2 |
| `--ceess-device DEVICE` | auto | torch device |
| `--ceess-cv-folds INT` | `5` | cross-validation folds |
| `--ceess-random-state INT` | `0` | random seed |
| `--ceess-epochs INT` | `200` | MLP training epochs |
| `--ceess-hidden-dim INT` | `128` | MLP hidden layer width |
| `--ceess-dropout FLOAT` | `0.1` | MLP dropout rate |
| `--ceess-learning-rate FLOAT` | `1e-3` | MLP learning rate |
| `--ceess-weight-decay FLOAT` | `1e-4` | MLP weight decay |
| `--ceess-train-batch-size INT` | `8` | MLP training batch size |
| `--ceess-barlow-representation-dim INT` | `None` | Barlow Twins encoder output width (`contrastive` only) |
| `--ceess-barlow-projection-dim INT` | `None` | Barlow Twins projection head width (`contrastive` only) |
| `--ceess-barlow-redundancy-weight FLOAT` | `0.005` | Barlow Twins redundancy penalty (`contrastive` only) |
| `--ceess-mlp-checkpoint PATH` | `None` | pretrained MLP checkpoint `.pt`; skips training |

---

## `ariadne phylogeny`

Stage 4 standalone. Requires `--candidates`, `--reference-dir`, and `--output-dir`.

| Parameter | Default | Description |
| --- | --- | --- |
| `--candidates PATH` | **required** | filtered candidate FASTA |
| `--reference-dir PATH` | **required** | reference FASTA directory |
| `--output-dir PATH` | **required** | phylogeny output directory |
| `--mafft-bin PATH` | auto | explicit MAFFT binary path |
| `--mafft-mode FLAG` | `--auto` | MAFFT alignment mode |
| `--iqtree-bin PATH` | auto | explicit IQ-TREE binary path |
| `--iqtree-model MODEL` | `LG` | IQ-TREE substitution model |
| `--iqtree-threads INT\|AUTO` | `AUTO` | IQ-TREE thread count |
| `--iqtree-bootstrap INT` | `None` | ultrafast bootstrap replicates |
| `--no-iqtree-fast` | `False` | disable IQ-TREE `--fast` mode |

---

## `ariadne prepare-references`

Prepare reference FASTA files from coral, insect, fungal, plant, and bacterial sources.

| Parameter | Default | Description |
| --- | --- | --- |
| `--coral PATH` | `None` | coral reference FASTA |
| `--coral-limit INT` | `None` | maximum number of coral reference sequences to include |
| `--insect-xlsx PATH` | `None` | insect TPS Excel workbook |
| `--insect-limit INT` | `None` | maximum number of insect sequences to include |
| `--bacteria-fasta PATH` | `None` | bacterial TPS FASTA |
| `--fungal-fasta PATH` | `None` | fungal TPS FASTA (alias: `--fungi-fasta`) |
| `--fungi-fasta PATH` | `None` | fungal TPS FASTA (alias: `--fungal-fasta`) |
| `--plant-fasta PATH` | `None` | plant TPS FASTA |
| `--extra-fasta PATH [PATH ...]` | `None` | additional FASTA files to include |
| `--output-dir PATH` | **required** | output directory for prepared reference files |

---

## `ariadne build-hmm`

Build a single profile HMM from an alignment or FASTA file.

| Parameter | Default | Description |
| --- | --- | --- |
| `--alignment PATH` | **required** | input alignment or FASTA file |
| `--output PATH` | **required** | output `.hmm` file path |
| `--name NAME` | `None` | name written into the HMM profile |

---

## `ariadne build-tps-hmm-library`

Build a directory of per-clade TPS HMMs from named alignment files.

| Parameter | Default | Description |
| --- | --- | --- |
| `--alignment NAME=PATH [NAME=PATH ...]` | **required** | one or more `NAME=PATH` pairs, where `NAME` becomes the HMM profile name and `PATH` is the input alignment |
| `--output-dir PATH` | **required** | output directory for the HMM library |

---

## Practical notes

- If `--query-hmm` is omitted, Ariadne prefers the bundled `ariadne/hmm/query.hmm` and only auto-builds from `tree/` as a fallback.
- If `--tps-hmm-dir` is omitted, Ariadne first uses the bundled `ariadne/hmm/` library and only auto-builds one HMM per FASTA file in `tree/` as a fallback.
- Global flags (`--verbose`, `--log-file`) must be placed **before** the subcommand: `ariadne --verbose run ...`
- If `TPS/TPS.xlsx` is present and the optional ESM stack is installed, `ariadne classify` and `ariadne run` automatically add a CeeSs scoring pass after HMM classification. Use `--skip-ceess-model` to suppress it.
- `--skip-phylogeny` is useful for fast iteration when only candidate ranking and embeddings are needed.
- The CeeSs head requires optional ESM dependencies. Install with `pip install 'ariadne-tps[esm]'` or `pip install torch transformers`.
- `--ceess-classifier contrastive` activates the Barlow Twins variant; the three `--ceess-barlow-*` flags are only used in this mode.
- `--ceess-mlp-checkpoint` loads a previously saved `.pt` file and skips all training; useful for reproducible inference.
