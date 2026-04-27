<div class="hero-panel">
  <div class="hero-copy">
    <h1>Ariadne</h1>
    <p><strong>A coral-centered terpene synthase discovery and CeeSs prioritization platform</strong> for TPS mining, feature-space interpretation, and phylogenetic analysis.</p>
    <p>Ariadne is a tree-native, four-stage computational platform. It combines profile HMM-guided discovery, coverage- and length-aware filtering, TPS feature-space embedding, and an optional ESM2-based scoring layer for CeeSs candidates.</p>
    <div class="hero-actions">
      <a class="md-button md-button--primary" href="getting-started/">Get Started</a>
      <a class="md-button" href="cli-reference/">CLI Reference</a>
      <a class="md-button" href="https://github.com/zhaorui-bi/Ariadne">GitHub</a>
    </div>
    <div class="hero-meta">
      <span class="hero-pill">Python >= 3.9</span>
      <span class="hero-pill">Optional ESM2 CeeSs</span>
      <span class="hero-pill">MAFFT + IQ-TREE</span>
      <span class="hero-pill">4-stage pipeline</span>
    </div>
  </div>
  <div class="hero-visual">
    <img src="assets/framework.png" alt="Ariadne overview figure">
  </div>
</div>

## 📢 News

- `2026-04-02` Architecture refactored to a clean 7-module layout: `utils`, `data`, `search`, `filter`, `embed`, `model`, `tree`. Old verbose names retired.
- `2026-04-02` Filtering updated: candidates matching reference sequences are now **retained** in `candidates.filtered.faa`; matches are still logged in `reference_matches.tsv` for traceability.
- `2026-04-02` Example run summary: 100 candidates discovered -> 36 retained after filtering -> 36 classified as coral-like -> **5 CeeSs candidates** shortlisted (P(CeeSs) >= 0.9).
- `2026-03-23` The repository introduction was updated around the new CeeSs framing, reflecting Ariadne as a platform for coral TPS mining and CeeSs prioritization.
- `2026-03-22` Ariadne now ships with an English documentation site built with MkDocs + Material for Read the Docs deployment.
- `2026-03-21` The software workflow was simplified to a focused four-stage pipeline: `discovery -> filtering -> classification -> phylogeny`.

## 🪸 Abstract

Coral terpene synthases (TPSs) represent an underexplored frontier in natural product biosynthesis. Identifying which of the hundreds of predicted coral TPS proteins is responsible for a specific terpenoid product, especially cembrene-class (CeeSs) compounds, requires more than sequence homology: it demands systematic embedding in a curated reference space.

Ariadne starts from raw transcriptomes or predicted proteomes and uses the same <code>tree/</code> reference backbone across discovery, classification, and final phylogenetic placement, keeping candidate interpretation biologically consistent from the first HMM search through the final IQ-TREE result.

## ✨ Why Ariadne?

<div class="card-grid card-grid--three">
  <div class="paper-card">
    <h3>🌊 Tree-native design</h3>
    <p>A single curated <code>tree/</code> directory drives discovery, classification, and phylogeny without manual bookkeeping between stages.</p>
  </div>
  <div class="paper-card">
    <h3>🧭 HMM feature space</h3>
    <p>Candidates are scored against a multi-source TPS HMM library, embedded with LDA/PCA, and assigned nearest reference neighbors.</p>
  </div>
  <div class="paper-card">
    <h3>🧬 ESM2 CeeSs scoring</h3>
    <p>When <code>TPS/TPS.xlsx</code> and the optional ESM stack are available, a frozen ESM2 backbone with a trainable head scores coral-like candidates.</p>
  </div>
</div>

<div class="mini-kpi">
  <div class="paper-card"><strong>4</strong><span>Stages</span></div>
  <div class="paper-card"><strong>1</strong><span>Reference Root</span></div>
  <div class="paper-card"><strong>2</strong><span>Input Modes</span></div>
  <div class="paper-card"><strong>SVG</strong><span>Publication Outputs</span></div>
</div>

## 🧠 Method At A Glance

<figure class="paper-figure">
  <img src="assets/framework.png" alt="Ariadne method overview">
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
  <img src="assets/latest_embedding.svg" alt="Latest Ariadne classification embedding output">
  <figcaption>
    Figure 2. Current bundled local result preview, synced from <code>output/03_classification/embedding.svg</code>. In this local run, <code>36</code> coral-like candidates were scored and <code>5</code> were retained as final CeeSs candidates (P(CeeSs) ≥ 0.9).
  </figcaption>
</figure>

## 🚀 Quick Start

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

This command will:

- use the bundled discovery query HMM from `ariadne/hmm/query.hmm`
- use the bundled TPS HMM library from `ariadne/hmm/`
- discover and filter candidates
- classify them in TPS feature space
- directly generate the final alignment and phylogeny

## 🧭 Where To Go Next

- Start with [Getting Started](getting-started.md) for installation and your first run.
- Read [Method](method.md) for a stage-by-stage explanation of the pipeline.
- Open [Tutorials](tutorials.md) for practical command sequences and example analysis flow.
- Use [CLI Reference](cli-reference.md) when tuning parameters.
- Check [Outputs](outputs.md) to understand every key artifact produced by the pipeline.
