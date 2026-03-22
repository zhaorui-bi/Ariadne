# FAQ

## Why is `tree/` so important?

`tree/` is the canonical reference directory in the current Ariadne release. It is reused for:

- discovery query-HMM generation
- TPS HMM library generation
- classification background sequences
- phylogeny background sequences

## Does Ariadne still use `Alignment.fasta`?

No. The current workflow is tree-native and no longer depends on the old `Alignment.fasta` entrypoints.

## Does Ariadne still include motif analysis?

No. Motif-based post-processing was intentionally removed from the current release.

## Does Ariadne still include benchmark mode?

No. Benchmark-vs-expected FASTA comparison was also removed from the active workflow.

## Which outputs should I inspect first?

For most users:

1. `classification.tsv`
2. `embedding.svg`
3. `iqtree.treefile`
4. `iqtree.iqtree`

## What if I already have my own HMMs?

You can pass:

- `--query-hmm` for a prebuilt discovery HMM
- `--tps-hmm-dir` for a prebuilt TPS HMM library

## What if I only want classification?

Use:

```bash
ariadne classify \
  --candidates my_candidates.faa \
  --reference-dir tree/ \
  --output-dir results_classification/
```

## What if MAFFT or IQ-TREE is not installed?

Discovery, filtering, and classification can still run, but the phylogeny stage requires:

- `mafft`
- `iqtree` or `iqtree2`

You can temporarily skip the final tree stage with:

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --skip-phylogeny \
  --output-dir results_no_tree/
```
