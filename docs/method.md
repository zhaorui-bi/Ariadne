# Method

## Conceptual framing

Ariadne is built on a single premise: **the same curated reference collection should support candidate discovery, candidate interpretation, and final phylogenetic placement.** In the current release that shared backbone is the `tree/` directory — a multi-clade collection of TPS reference sequences (coral, insect, plant, fungi, bacteria).

Coupling all stages to one reference universe keeps the discovery query, the classification feature space, and the phylogenetic background biologically consistent, which is what makes downstream interpretation tractable for coral TPS mining and CeeSs prioritization.

<figure class="paper-figure">
  <img src="assets/overview_pipeline.svg" alt="Ariadne pipeline method figure">
  <figcaption>
    Figure 1. Ariadne is a four-stage workflow driven by a single multi-clade TPS reference directory.
  </figcaption>
</figure>

The tree-native design has three direct consequences:

- the discovery query HMM and the classification feature space are derived from the same reference universe;
- the final phylogeny is reconstructed against the same multi-clade background that informed screening;
- nearest-neighbor assignments and phylogenetic placement can be compared within one consistent context.

## Stage I — Discovery

The first stage is tuned for sensitivity. Unless a prebuilt `--query-hmm` is supplied, Ariadne builds a discovery profile HMM from the coral reference in `tree/`. Two entry modes are supported:

- **protein mode** (`--protein-folder`) — input proteins are searched directly;
- **transcriptome mode** (`--transcriptomes`) — ORFs are first predicted with Pyrodigal (meta mode), then the translated proteins are searched.

Profile-HMM search (via `pyhmmer`) yields the candidate universe refined downstream. Optional `--discovery-min-score` and `--discovery-max-evalue` cutoffs gate the hits.

**Primary outputs:** `candidates.protein.faa`, `candidates.orf.fna`, `candidates.hits.tsv`.

## Stage II — Filtering

The second stage enforces candidate quality without imposing biological interpretation. It applies, in order:

1. **coverage filtering** (default ≥ 10×, parsed from the FASTA header);
2. **minimum-length filtering** (default ≥ 300 aa);
3. **near-duplicate collapsing** at 95% identity, using a bounded edit-distance test so that long sequences are compared efficiently.

Candidates with ≥ 95% identity to any reference in `tree/` are **retained** in `candidates.filtered.faa`; their matches are recorded in `reference_matches.tsv` for traceability. This is deliberate — novel alleles and species-specific variants of known coral TPS genes proceed through classification and CeeSs scoring rather than being silently discarded.

The stage is conservative and transparent: rather than hiding decisions inside a single score, Ariadne exports per-candidate reports that make selection and deduplication explicit.

**Primary outputs:** `candidates.filtered.faa`, `filter_report.tsv`, `dedupe_clusters.tsv`, `reference_matches.tsv`, `manual_review.tsv`.

## Stage III — Classification

The third stage is the interpretive core of the workflow. It places every filtered candidate into a TPS feature space defined by the reference collection, in four steps:

1. **Feature construction.** All filtered candidates and all reference sequences are scored against a multi-clade TPS profile-HMM library derived from `tree/`. Each sequence becomes a feature vector of per-profile bit scores.
2. **Normalization.** The feature matrix is z-scored across profiles so that no single HMM dominates the geometry.
3. **Low-dimensional embedding.** The matrix is projected to 3-D by supervised linear discriminant analysis (LDA), using *k*-means subclustering of the large coral reference set to define discriminative classes; principal component analysis (PCA) is used as a fallback when supervision is degenerate.
4. **Nearest-reference label transfer.** Each candidate is assigned a label by *k*-nearest-neighbor voting (`--top-k`, default 5) over reference neighbors in the embedding, and a local context tree is built from its `--tree-neighbors` (default 12) nearest references.

The result is not a single label but a structured evidence layer: embedding coordinates, nearest neighbors, local and global context trees, and the per-sequence feature matrix.

**Primary outputs:** `tps_features.tsv`, `embedding.tsv`, `classification.tsv`, `nearest_neighbors.tsv`, `candidate_cluster_context.tsv`, `embedding.svg`, `embedding_3d_sections.svg`, `global_context_tree.nwk`.

### Optional ESM2 CeeSs scoring

When `TPS/TPS.xlsx` and the optional ESM dependencies are available, classification is followed by a protein-language-model scoring pass over the coral-like candidates:

1. compute frozen, mean-pooled ESM2 embeddings for the labeled training proteins and the candidates;
2. train a lightweight head — an MLP (default), logistic regression, or a Barlow Twins contrastive variant — on the training embeddings while keeping the ESM2 backbone frozen;
3. predict a fine-grained TPS type for each candidate;
4. aggregate `P(CeeSs)` as the summed probability over all workbook-defined CeeSs-positive labels and apply `--ceess-threshold`.

This places candidates inside the broader TPS landscape *before* phylogenetic reconstruction. The full scoring design and output schema are documented in [CeeSs Classifier](esm-type.md).

## Stage IV — Phylogeny

The fourth stage converts the screened candidate set into a phylogeny-ready analysis object. Filtered candidates are merged with the references loaded from `tree/`, deduplicated, and then:

1. aligned with **MAFFT** (`--mafft-mode`, default `--auto`);
2. reconstructed into a maximum-likelihood tree with **IQ-TREE** (`--iqtree-model`, default `LG`; optional ultrafast bootstrap via `--iqtree-bootstrap`).

A compact SVG preview is rendered directly from the resulting Newick tree, providing the bridge from candidate discovery to evolutionary interpretation.

**Primary outputs:** `phylogeny_input.fasta`, `phylogeny_alignment.fasta`, `phylogeny_sequence_map.tsv`, `iqtree.treefile`, `iqtree.iqtree`, `phylogeny_preview.svg`.

## Why this structure works

<div class="card-grid">
  <div class="paper-card">
    <h3>Reference consistency</h3>
    <p>One tree-native reference collection is reused across discovery, classification, and phylogeny, so evidence from each stage is directly comparable.</p>
  </div>
  <div class="paper-card">
    <h3>Screening-to-evolution continuity</h3>
    <p>Nearest-neighbor evidence, embeddings, and phylogenetic outputs live in one coherent result directory, end to end.</p>
  </div>
</div>

## Scope of the current release

The active release intentionally focuses on a stable, interpretable four-stage pipeline and excludes two earlier experimental paths: motif-centric post-filtering, and benchmark-versus-expected FASTA comparison. Removing them keeps the implementation focused and the outputs straightforward to reason about.
