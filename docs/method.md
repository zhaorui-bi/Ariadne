# Method

## Conceptual Framing

Ariadne is built around a user-provided reference FASTA directory. That directory should contain the multi-clade TPS references used for classification and visualization. The optional root-level `TPS.xlsx` workbook is used by the ESM2 CeeSs scoring layer.

The current main workflow has four stages:

1. candidate discovery;
2. candidate filtering;
3. feature-space classification;
4. PCA/LDA visualization and candidate triage.

The last stage is no longer a tree-building step. The documented workflow does not require users to generate HMM files from the reference FASTA directory.

## Stage I - Discovery

The first stage creates the candidate universe from either:

- **protein mode** (`--protein-folder`) - input proteins are searched directly;
- **transcriptome mode** (`--transcriptomes`) - ORFs are first predicted with Pyrodigal, then translated proteins are screened.

Optional `--discovery-min-score` and `--discovery-max-evalue` cutoffs can be used to tune sensitivity.

**Primary outputs:** `candidates.protein.faa`, `candidates.orf.fna`, `candidates.hits.tsv`.

## Stage II - Filtering

Filtering applies three transparent quality steps:

1. **coverage filtering** (default `>= 10x`, parsed from the FASTA header);
2. **minimum-length filtering** (default `>= 300 aa`);
3. **near-duplicate collapsing** at 95% identity.

Candidates that match a reference sequence are retained in `candidates.filtered.faa` and logged in `reference_matches.tsv`. This keeps known-like coral TPS alleles visible for later interpretation.

**Primary outputs:** `candidates.filtered.faa`, `filter_report.tsv`, `dedupe_clusters.tsv`, `reference_matches.tsv`, `manual_review.tsv`.

## Stage III - Classification

Classification places every filtered candidate into a shared TPS feature space defined by the supplied reference FASTA directory:

1. build a per-sequence feature matrix for candidates and references;
2. normalize the matrix so features are comparable;
3. assign candidate labels by nearest-reference voting;
4. write neighbor evidence and assignment summaries.

The result is a structured evidence layer rather than only a final label.

**Primary outputs:** `tps_features.tsv`, `classification.tsv`, `nearest_neighbors.tsv`, `candidate_cluster_context.tsv`, `assignment_summary.tsv`.

### Optional ESM2 CeeSs Scoring

When `TPS.xlsx` and the optional ESM dependencies are available, classification is followed by a protein-language-model scoring pass:

1. compute frozen, mean-pooled ESM2 embeddings for labeled training proteins and candidates;
2. train a lightweight classifier head;
3. predict a fine-grained TPS type;
4. aggregate `P(CeeSs)` and apply `--ceess-threshold`.

The full scoring design and output schema are documented in [CeeSs Classifier](esm-type.md).

## Stage IV - PCA/LDA Visualization

The fourth stage projects the same feature space for visual inspection. Ariadne prefers supervised LDA when labels support it and falls back to PCA when the supervised projection is not well-defined.

The visualization layer is produced under `03_classification/`:

- `embedding.tsv` - coordinates and labels;
- `embedding_variance.tsv` - explained variance / projection metadata;
- `embedding.svg` - 2-D visual inspection plot;
- `embedding_3d_sections.svg` - three orthogonal 3-D projection panels;
- `candidate_cluster_context.tsv` - per-candidate placement context.

## Why This Structure Works

<div class="card-grid">
  <div class="paper-card">
    <h3>Reference consistency</h3>
    <p>The same reference FASTAs support classification, visualization, and optional CeeSs scoring.</p>
  </div>
  <div class="paper-card">
    <h3>Fast triage</h3>
    <p>The main workflow stops at tables and SVGs that are quick to inspect, rerun, and archive.</p>
  </div>
</div>

## Scope Of The Current Release

The active release focuses on interpretable candidate discovery, classification, PCA/LDA visualization, and optional CeeSs scoring. The old final tree/phylogeny section is not part of the main documented workflow.
