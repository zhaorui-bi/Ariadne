# FAQ

## Why is `tree/` so important?

`tree/` is the canonical reference directory and the backbone of the whole pipeline. The same collection is reused for:

- discovery query-HMM generation;
- TPS HMM library generation;
- classification background sequences;
- phylogeny background sequences.

Keeping every stage anchored to one reference universe is what makes Ariadne's outputs comparable across stages.

## What inputs does Ariadne accept?

Either predicted proteins (`--protein-folder`, recursive `.faa` / `.fa` / `.fasta`) or transcriptome assemblies (`--transcriptomes`), in which case ORFs are predicted with Pyrodigal before HMM search.

## What if I already have my own HMMs?

Pass `--query-hmm` for a prebuilt discovery HMM and/or `--tps-hmm-dir` for a prebuilt TPS HMM library. Otherwise Ariadne uses the bundled HMMs in `ariadne/hmm/`, falling back to building them from `tree/`.

## What if I only want classification?

```bash
ariadne classify \
  --candidates my_candidates.faa \
  --reference-dir tree/ \
  --output-dir results_classification/
```

## What if MAFFT or IQ-TREE is not installed?

Discovery, filtering, and classification still run; only the phylogeny stage requires `mafft` and `iqtree`/`iqtree2`. Skip the tree with `--skip-phylogeny`:

```bash
ariadne run --protein-folder input/ --reference-dir tree/ \
  --skip-phylogeny --output-dir results_no_tree/
```

## Do I need a GPU for CeeSs scoring?

No. The ESM2 backbone runs on CPU, just more slowly. Select a device explicitly with `--ceess-device cuda:0` or `--ceess-device cpu`; by default Ariadne auto-detects.

## Which outputs should I inspect first?

For most users: `classification.tsv`, then `embedding.svg`, then `iqtree.treefile` and `iqtree.iqtree`.

## Does Ariadne still use `Alignment.fasta`, motif analysis, or benchmark mode?

No. The current release is tree-native and no longer depends on the old `Alignment.fasta` entrypoints; motif-based post-processing and benchmark-versus-expected comparison were both removed to keep the workflow focused (see [Method → Scope](method.md#scope-of-the-current-release)).
