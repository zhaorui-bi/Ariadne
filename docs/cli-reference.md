# CLI Reference

## Command map

| Command | Purpose |
| --- | --- |
| `ariadne run` | full end-to-end workflow |
| `ariadne discover` | HMM-based candidate discovery |
| `ariadne filter` | quality control and near-duplicate collapsing |
| `ariadne classify` | TPS feature-space classification |
| `ariadne phylogeny` | MAFFT alignment and IQ-TREE reconstruction |
| `ariadne prepare-references` | prepare reference FASTA files from source data |
| `ariadne build-hmm` | build a single HMM from a FASTA/MSA source |
| `ariadne build-tps-hmm-library` | build a directory of TPS HMMs from references |

## `ariadne run`

| Parameter | Default | Description |
| --- | --- | --- |
| `--protein-folder` | `None` | directory of protein FASTA files |
| `--transcriptomes` | `None` | transcriptome FASTA inputs |
| `--protein-glob` | common FASTA globs | recursive file patterns under `--protein-folder` |
| `--query-hmm` | bundled `ariadne/hmm/query.hmm` | explicit discovery HMM |
| `--reference-dir` | required | tree-native reference FASTA directory |
| `--output-dir` | required | output root |
| `--hmm-name` | `ariadne_query` | name for an auto-built discovery HMM |
| `--discovery-min-score` | `None` | optional HMM score cutoff |
| `--discovery-max-evalue` | `None` | optional E-value cutoff |
| `--min-coverage` | `10.0` | filtering coverage threshold |
| `--min-length` | `300` | filtering minimum protein length |
| `--identity-threshold` | `0.95` | near-duplicate collapsing threshold, also used for optional reference-match removal |
| `--tps-hmm-dir` | auto | explicit TPS HMM library directory |
| `--top-k` | `5` | nearest-reference voting size |
| `--tree-neighbors` | `12` | neighbors used for local candidate context trees |
| `--ceess-xlsx` | `TPS/TPS.xlsx` when present | labeled coral TPS workbook for the optional ESM CeeSs head |
| `--skip-ceess-model` | `False` | skip the optional ESM-based CeeSs scoring stage |
| `--ceess-threshold` | `0.9` | minimum aggregated `P(CeeSs)` from the multi-class ESM head for a final CeeSs candidate |
| `--ceess-classifier` | `mlp` | classifier head on top of frozen ESM embeddings (`mlp` or `logreg`) |
| `--ceess-model-name` | `facebook/esm2_t33_650M_UR50D` | ESM preset or Hugging Face model id for CeeSs scoring |
| `--ceess-batch-size` | `4` | ESM inference batch size for the CeeSs head |
| `--ceess-max-length` | `2048` | maximum tokenized sequence length for the CeeSs head |
| `--ceess-device` | auto | explicit torch device for the CeeSs head |
| `--ceess-cv-folds` | `5` | target cross-validation folds for training-set evaluation |
| `--ceess-random-state` | `0` | random seed for the CeeSs head |
| `--ceess-epochs` | `200` | training epochs for the MLP CeeSs head |
| `--ceess-hidden-dim` | `128` | hidden layer width for the MLP CeeSs head |
| `--ceess-dropout` | `0.1` | dropout rate for the MLP CeeSs head |
| `--ceess-learning-rate` | `1e-3` | learning rate for the MLP CeeSs head |
| `--ceess-weight-decay` | `1e-4` | weight decay for the MLP CeeSs head |
| `--ceess-train-batch-size` | `8` | classifier-training batch size for the MLP CeeSs head |
| `--ceess-barlow-representation-dim` | `None` | encoded representation width for the Barlow Twins pipeline (contrastive only) |
| `--ceess-barlow-projection-dim` | `None` | projection-head output width for the Barlow Twins pipeline (contrastive only) |
| `--ceess-barlow-redundancy-weight` | `0.005` | off-diagonal redundancy penalty for Barlow Twins loss (contrastive only) |
| `--ceess-mlp-checkpoint` | `None` | path to a pretrained MLP checkpoint (`.pt`) — skips training and loads directly |
| `--skip-phylogeny` | `False` | skip MAFFT + IQ-TREE |
| `--mafft-bin` | auto | explicit MAFFT binary |
| `--mafft-mode` | `--auto` | MAFFT mode flag |
| `--iqtree-bin` | auto | explicit IQ-TREE binary |
| `--iqtree-model` | `LG` | substitution model |
| `--iqtree-threads` | `AUTO` | thread setting |
| `--iqtree-bootstrap` | `None` | optional bootstrap replicates |
| `--no-iqtree-fast` | `False` | disable default `--fast` mode |

## `ariadne discover`

| Parameter | Default | Description |
| --- | --- | --- |
| `--protein-folder` | `None` | protein FASTA directory |
| `--transcriptomes` | `None` | transcript FASTA inputs |
| `--protein-glob` | common FASTA globs | recursive search patterns |
| `--hmm` | required | discovery HMM |
| `--output-dir` | required | discovery output directory |
| `--min-score` | `None` | optional score cutoff |
| `--max-evalue` | `None` | optional E-value cutoff |

## `ariadne filter`

| Parameter | Default | Description |
| --- | --- | --- |
| `--input-fasta` | required | candidate protein FASTA |
| `--output-dir` | required | filter output directory |
| `--min-coverage` | `10.0` | remove low-coverage candidates |
| `--min-length` | `300` | remove short proteins |
| `--identity-threshold` | `0.95` | near-duplicate threshold |
| `--reference-dir` | `None` | optional reference FASTA directory — candidates matching reference sequences are logged in `reference_matches.tsv` but are **kept** in the filtered output |

## `ariadne classify`

| Parameter | Default | Description |
| --- | --- | --- |
| `--candidates` | required | filtered candidate FASTA |
| `--reference-dir` | required | reference FASTA directory |
| `--output-dir` | required | classification output directory |
| `--tps-hmm-dir` | auto | TPS HMM library directory |
| `--top-k` | `5` | voting neighbors |
| `--tree-neighbors` | `12` | local context-tree neighbors |
| `--ceess-xlsx` | `TPS/TPS.xlsx` when present | labeled coral TPS workbook for the optional ESM CeeSs head |
| `--skip-ceess-model` | `False` | skip the optional ESM-based CeeSs scoring stage |
| `--ceess-threshold` | `0.9` | minimum aggregated `P(CeeSs)` from the multi-class ESM head for a final CeeSs candidate |
| `--ceess-classifier` | `mlp` | classifier head on top of frozen ESM embeddings (`mlp` or `logreg`) |
| `--ceess-model-name` | `facebook/esm2_t33_650M_UR50D` | ESM preset or Hugging Face model id for CeeSs scoring |
| `--ceess-batch-size` | `4` | ESM inference batch size for the CeeSs head |
| `--ceess-max-length` | `2048` | maximum tokenized sequence length for the CeeSs head |
| `--ceess-device` | auto | explicit torch device for the CeeSs head |
| `--ceess-cv-folds` | `5` | target cross-validation folds for training-set evaluation |
| `--ceess-random-state` | `0` | random seed for the CeeSs head |
| `--ceess-epochs` | `200` | training epochs for the MLP CeeSs head |
| `--ceess-hidden-dim` | `128` | hidden layer width for the MLP CeeSs head |
| `--ceess-dropout` | `0.1` | dropout rate for the MLP CeeSs head |
| `--ceess-learning-rate` | `1e-3` | learning rate for the MLP CeeSs head |
| `--ceess-weight-decay` | `1e-4` | weight decay for the MLP CeeSs head |
| `--ceess-train-batch-size` | `8` | classifier-training batch size for the MLP CeeSs head |
| `--ceess-barlow-representation-dim` | `None` | encoded representation width for the Barlow Twins pipeline (contrastive only) |
| `--ceess-barlow-projection-dim` | `None` | projection-head output width for the Barlow Twins pipeline (contrastive only) |
| `--ceess-barlow-redundancy-weight` | `0.005` | off-diagonal redundancy penalty for Barlow Twins loss (contrastive only) |
| `--ceess-mlp-checkpoint` | `None` | path to a pretrained MLP checkpoint (`.pt`) — skips training and loads directly |

## `ariadne phylogeny`

| Parameter | Default | Description |
| --- | --- | --- |
| `--candidates` | required | filtered candidate FASTA |
| `--reference-dir` | required | reference FASTA directory |
| `--output-dir` | required | phylogeny output directory |
| `--mafft-bin` | auto | explicit MAFFT binary |
| `--mafft-mode` | `--auto` | MAFFT mode |
| `--iqtree-bin` | auto | explicit IQ-TREE binary |
| `--iqtree-model` | `LG` | substitution model |
| `--iqtree-threads` | `AUTO` | threads |
| `--iqtree-bootstrap` | `None` | optional bootstrap |
| `--no-iqtree-fast` | `False` | disable fast mode |

## Practical notes

- If `--query-hmm` is omitted, Ariadne prefers the bundled `ariadne/hmm/query.hmm` and only auto-builds from `tree/` as a fallback.
- If `--tps-hmm-dir` is omitted, Ariadne first uses the bundled `ariadne/hmm/` library and only auto-builds one HMM per FASTA file in `tree/` as a fallback.
- If `TPS/TPS.xlsx` is present and the optional ESM stack is installed, `ariadne classify` and `ariadne run` automatically add a second-stage CeeSs scorer on top of the normal coral classification step.
- `--skip-phylogeny` is useful for fast iteration when you only need candidate ranking and embeddings.
- The integrated CeeSs head requires the optional ESM dependencies. Install them with `pip install 'ariadne-tps[esm]'` or `pip install torch transformers`.
