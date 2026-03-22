# CLI Reference

## Command map

| Command | Purpose |
| --- | --- |
| `ariadne run` | full end-to-end workflow |
| `ariadne discover` | HMM-based candidate discovery |
| `ariadne filter` | quality control and near-duplicate collapsing |
| `ariadne classify` | TPS feature-space classification |
| `ariadne phylogeny` | MAFFT alignment and IQ-TREE reconstruction |
| `ariadne build-hmm` | build a single HMM from a FASTA/MSA source |
| `ariadne build-tps-hmm-library` | build a directory of TPS HMMs from references |

## `ariadne run`

| Parameter | Default | Description |
| --- | --- | --- |
| `--protein-folder` | `None` | directory of protein FASTA files |
| `--transcriptomes` | `None` | transcriptome FASTA inputs |
| `--protein-glob` | common FASTA globs | recursive file patterns under `--protein-folder` |
| `--query-hmm` | auto | explicit discovery HMM |
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

- If `--query-hmm` is omitted, Ariadne auto-builds the discovery HMM from the coral reference in `tree/`.
- If `--tps-hmm-dir` is omitted, Ariadne auto-builds one HMM per FASTA file in `tree/`.
- `--skip-phylogeny` is useful for fast iteration when you only need candidate ranking and embeddings.
