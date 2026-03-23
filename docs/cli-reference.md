# CLI Reference

## Command map

| Command | Purpose |
| --- | --- |
| `ariadne run` | full end-to-end workflow |
| `ariadne discover` | HMM-based candidate discovery |
| `ariadne filter` | quality control and near-duplicate collapsing |
| `ariadne classify` | TPS feature-space classification |
| `ariadne phylogeny` | MAFFT alignment and IQ-TREE reconstruction |
| `ariadne esm-type` | ESM2 embedding and supervised coral TPS type classification from `TPS.xlsx` |
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
| `--identity-threshold` | `0.95` | near-duplicate collapsing threshold |
| `--tps-hmm-dir` | auto | explicit TPS HMM library directory |
| `--top-k` | `5` | nearest-reference voting size |
| `--tree-neighbors` | `12` | neighbors used for local candidate context trees |
| `--ceess-xlsx` | `TPS/TPS.xlsx` when present | labeled coral TPS workbook for the optional ESM CeeSs head |
| `--skip-ceess-model` | `False` | skip the optional ESM-based CeeSs scoring stage |
| `--ceess-threshold` | `0.5` | minimum combined `cembrene A + cembrene B` probability for a final CeeSs candidate |
| `--ceess-model-name` | `facebook/esm2_t33_650M_UR50D` | ESM preset or Hugging Face model id for coral TPS type scoring |
| `--ceess-batch-size` | `4` | ESM inference batch size for the CeeSs head |
| `--ceess-max-length` | `2048` | maximum tokenized sequence length for the CeeSs head |
| `--ceess-device` | auto | explicit torch device for the CeeSs head |
| `--ceess-cv-folds` | `5` | target cross-validation folds for training-set evaluation |
| `--ceess-random-state` | `0` | random seed for the CeeSs head |
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
| `--ceess-threshold` | `0.5` | minimum combined `cembrene A + cembrene B` probability for a final CeeSs candidate |
| `--ceess-model-name` | `facebook/esm2_t33_650M_UR50D` | ESM preset or Hugging Face model id for coral TPS type scoring |
| `--ceess-batch-size` | `4` | ESM inference batch size for the CeeSs head |
| `--ceess-max-length` | `2048` | maximum tokenized sequence length for the CeeSs head |
| `--ceess-device` | auto | explicit torch device for the CeeSs head |
| `--ceess-cv-folds` | `5` | target cross-validation folds for training-set evaluation |
| `--ceess-random-state` | `0` | random seed for the CeeSs head |

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

## `ariadne esm-type`

| Parameter | Default | Description |
| --- | --- | --- |
| `--xlsx` | `TPS/TPS.xlsx` when present | labeled coral TPS spreadsheet |
| `--output-dir` | required | ESM analysis output directory |
| `--sheet-name` | first sheet | optional worksheet selection |
| `--model-name` | `facebook/esm2_t33_650M_UR50D` | ESM preset or Hugging Face model id |
| `--batch-size` | `4` | ESM inference batch size |
| `--max-length` | `2048` | maximum tokenized sequence length |
| `--device` | auto | explicit torch device |
| `--cv-folds` | `5` | target cross-validation folds |
| `--random-state` | `0` | random seed for evaluation splits |

## Practical notes

- If `--query-hmm` is omitted, Ariadne prefers the bundled `ariadne/hmm/query.hmm` and only auto-builds from `tree/` as a fallback.
- If `--tps-hmm-dir` is omitted, Ariadne first uses the bundled `ariadne/hmm/` library and only auto-builds one HMM per FASTA file in `tree/` as a fallback.
- If `TPS/TPS.xlsx` is present and the optional ESM stack is installed, `ariadne classify` and `ariadne run` automatically add a second-stage CeeSs scorer on top of the normal coral classification step.
- `--skip-phylogeny` is useful for fast iteration when you only need candidate ranking and embeddings.
- `ariadne esm-type` and the integrated CeeSs head require the optional ESM dependencies. Install them with `pip install 'ariadne-tps[esm]'` or `pip install torch transformers`.
