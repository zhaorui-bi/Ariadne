# Outputs

## Output root

A standard `ariadne run` creates:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
├── 04_phylogeny/
└── pipeline_summary.tsv
```

## `01_discovery`

Key files:

- `all_predicted_proteins.faa`
- `candidates.protein.faa`
- `candidates.orf.fna`
- `candidates.hits.tsv`

What to inspect:

- use `candidates.hits.tsv` to check which records passed HMM discovery
- use `candidates.protein.faa` as the starting point for downstream filtering

## `02_filtering`

Key files:

- `candidates.filtered.faa`
- `filter_report.tsv`
- `dedupe_clusters.tsv`
- `manual_review.tsv`

What to inspect:

- `filter_report.tsv` explains why records were removed
- `dedupe_clusters.tsv` records representative-member relationships
- `candidates.filtered.faa` is the canonical handoff into classification and phylogeny

## `03_classification`

Key files:

- `tps_features.tsv`
- `embedding.tsv`
- `embedding.svg`
- `embedding_3d_sections.svg`
- `classification.tsv`
- `nearest_neighbors.tsv`
- `candidate_cluster_context.tsv`
- `global_context_tree.nwk`

What to inspect:

- `classification.tsv`
  best single-file summary of candidate predictions
- `nearest_neighbors.tsv`
  evidence for the predicted source assignment
- `embedding.svg`
  fast visual inspection of candidate placement
- `embedding_3d_sections.svg`
  publication-style multi-view embedding figure

## `04_phylogeny`

Key files:

- `phylogeny_input.fasta`
- `phylogeny_alignment.fasta`
- `phylogeny_sequence_map.tsv`
- `iqtree.treefile`
- `iqtree.iqtree`
- `iqtree.log`

What to inspect:

- `phylogeny_sequence_map.tsv`
  maps tree-safe identifiers back to original sequence headers
- `iqtree.treefile`
  final phylogeny in Newick format
- `iqtree.iqtree`
  model and inference summary from IQ-TREE

## Reading order

If you are new to Ariadne, this reading order usually works best:

1. `pipeline_summary.tsv`
2. `03_classification/classification.tsv`
3. `03_classification/embedding.svg`
4. `04_phylogeny/iqtree.treefile`
5. `04_phylogeny/iqtree.iqtree`

## Representative preview

<figure class="paper-figure">
  <img src="assets/results_preview.svg" alt="Representative outputs preview">
  <figcaption>
    Figure 2. The most user-facing artifacts are the classification embedding panels and the final phylogeny outputs.
  </figcaption>
</figure>
