<div class="hero-panel">
  <div class="hero-copy">
    <h1>Ariadne</h1>
    <p><strong>A tree-native terpene synthase discovery and phylogeny platform</strong> for coral-centered and cross-clade TPS mining.</p>
    <p>Ariadne turns a curated <code>tree/</code> reference directory into a practical four-stage workflow: HMM-guided discovery, filtering, TPS feature-space classification, and alignment-driven phylogeny.</p>
    <div class="hero-actions">
      <a class="md-button md-button--primary" href="getting-started/">Get Started</a>
      <a class="md-button" href="cli-reference/">CLI Reference</a>
      <a class="md-button" href="https://github.com/zhaoruijiang26/Ariadne">GitHub</a>
    </div>
    <div class="hero-meta">
      <span class="hero-pill">Python 3.11+</span>
      <span class="hero-pill">MkDocs + Material</span>
      <span class="hero-pill">MAFFT + IQ-TREE</span>
      <span class="hero-pill">4-stage pipeline</span>
    </div>
  </div>
  <div class="hero-visual">
    <img src="assets/overview_pipeline.svg" alt="Ariadne overview figure">
  </div>
</div>

## 📢 News

- `2026-03-22` Ariadne now ships with an English documentation site built with MkDocs + Material for Read the Docs deployment.
- `2026-03-21` The software workflow was simplified to a focused four-stage pipeline: `discovery -> filtering -> classification -> phylogeny`.
- `2026-03-21` `tree/` became the default source for query-HMM construction, TPS HMM library generation, classification references, and phylogeny references.

## ✨ Why Ariadne?

<div class="card-grid card-grid--three">
  <div class="paper-card">
    <h3>🌊 Tree-native by design</h3>
    <p>A single curated <code>tree/</code> directory drives discovery, feature-space classification, and phylogenetic reconstruction.</p>
  </div>
  <div class="paper-card">
    <h3>🧭 Feature-space aware</h3>
    <p>Candidates are embedded in a TPS HMM score space, allowing fast nearest-reference assignment and visual screening.</p>
  </div>
  <div class="paper-card">
    <h3>🌳 Phylogeny-ready outputs</h3>
    <p>After classification, Ariadne directly builds a MAFFT alignment and an IQ-TREE phylogeny without extra manual glue code.</p>
  </div>
</div>

<div class="mini-kpi">
  <div class="paper-card"><strong>4</strong><span>Stages</span></div>
  <div class="paper-card"><strong>1</strong><span>Reference Root</span></div>
  <div class="paper-card"><strong>3</strong><span>Core Output Types</span></div>
  <div class="paper-card"><strong>0</strong><span>Benchmark Dependency</span></div>
</div>

## 🧠 Method At A Glance

<figure class="paper-figure">
  <img src="assets/overview_pipeline.svg" alt="Ariadne method overview">
  <figcaption>
    Figure 1. Ariadne uses a four-stage, tree-native workflow. The same <code>tree/</code> reference collection is reused across discovery, classification, and phylogeny.
  </figcaption>
</figure>

<div class="overview-grid">
  <div class="paper-card">
    <h3>01. Discovery</h3>
    <p>Build a discovery HMM from the default coral reference under <code>tree/</code>, then search protein inputs or transcriptome-derived ORFs for TPS candidates.</p>
  </div>
  <div class="paper-card">
    <h3>02. Filtering</h3>
    <p>Apply coverage filtering, minimum-length checks, and near-duplicate collapsing to keep a clean candidate set.</p>
  </div>
  <div class="paper-card">
    <h3>03. Classification</h3>
    <p>Score all references and candidates against a TPS HMM library, embed them in feature space, and assign the most likely reference source.</p>
  </div>
  <div class="paper-card">
    <h3>04. Phylogeny</h3>
    <p>Merge filtered candidates with reference sequences, run MAFFT, and reconstruct the final phylogeny using IQ-TREE.</p>
  </div>
</div>

## 🖼️ Results Preview

<figure class="paper-figure">
  <img src="assets/results_preview.svg" alt="Representative Ariadne outputs">
  <figcaption>
    Figure 2. Conference-style preview of the two most important result families: TPS feature-space embedding and final alignment-driven phylogeny.
  </figcaption>
</figure>

<div class="preview-grid">
  <div class="paper-card">
    <h3>Embedding panel</h3>
    <p><code>embedding.svg</code> and <code>embedding_3d_sections.svg</code> summarize how candidates relate to the reference clades in TPS HMM feature space.</p>
  </div>
  <div class="paper-card">
    <h3>Tree panel</h3>
    <p><code>phylogeny_alignment.fasta</code>, <code>iqtree.treefile</code>, and <code>iqtree.iqtree</code> provide the alignment and the final tree for downstream interpretation.</p>
  </div>
</div>

## 🚀 Quick Start

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

This command will:

- auto-build the discovery query HMM from `tree/`
- auto-build the TPS HMM library from `tree/`
- discover and filter candidates
- classify them in TPS feature space
- directly generate the final alignment and phylogeny

## 🧭 Where To Go Next

- Start with [Getting Started](getting-started.md) for installation and your first run.
- Read [Method](method.md) for a stage-by-stage explanation of the pipeline.
- Open [Tutorials](tutorials.md) for practical command sequences and example analysis flow.
- Use [CLI Reference](cli-reference.md) when tuning parameters.
- Check [Outputs](outputs.md) to understand every key artifact produced by the pipeline.
