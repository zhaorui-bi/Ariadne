<div class="hero-panel">
  <div class="hero-copy">
    <h1>Ariadne</h1>
    <p><strong>Terpene synthase discovery platform for coral TPS mining and cembrene-class synthase (CeeSs) prioritization.</strong></p>
    <p>Ariadne turns a user-provided reference FASTA directory into a four-stage workflow: candidate discovery, quality filtering, feature-space classification, and PCA/LDA visualization. The current main workflow ends with interpretable candidate triage, not a final tree-building stage.</p>
    <div class="hero-actions">
      <a class="md-button md-button--primary" href="getting-started/">Get Started</a>
      <a class="md-button" href="method/">Method</a>
      <a class="md-button" href="https://github.com/zhaorui-bi/Ariadne">GitHub</a>
    </div>
    <div class="hero-meta">
      <span class="hero-pill">Python 3.9+</span>
      <span class="hero-pill">reference FASTA directory</span>
      <span class="hero-pill">PCA/LDA visualization</span>
      <span class="hero-pill">Optional ESM2 scoring</span>
    </div>
  </div>
  <div class="hero-visual">
    <img src="fig/logo.png" alt="Ariadne logo">
  </div>
</div>

## Background

Terpene synthases (TPSs) produce a large and chemically diverse family of natural products. Coral genomes encode many candidate TPS proteins, while product-specific assignment remains difficult. Ariadne focuses on making this screening step traceable: every candidate is carried through filtering, nearest-reference classification, geometric visualization, and optional CeeSs scoring.

## Design Principles

<div class="card-grid card-grid--three">
  <div class="paper-card">
    <h3>Data-first workflow</h3>
    <p>Use <code>--reference-dir</code> to provide reference FASTAs. Users do not need to generate separate HMM resources for the documented workflow.</p>
  </div>
  <div class="paper-card">
    <h3>Feature-space aware</h3>
    <p>Candidates are embedded with reference sequences, assigned nearest-reference labels, and reported with supporting neighbors.</p>
  </div>
  <div class="paper-card">
    <h3>Visualization-led triage</h3>
    <p>Stage 4 is PCA/LDA visualization: <code>embedding.svg</code>, <code>embedding_3d_sections.svg</code>, variance, and cluster context.</p>
  </div>
</div>

<div class="mini-kpi">
  <div class="paper-card"><strong>4</strong><span>Pipeline Stages</span></div>
  <div class="paper-card"><strong>1</strong><span>Reference Directory</span></div>
  <div class="paper-card"><strong>5</strong><span>Reference Clades</span></div>
  <div class="paper-card"><strong>0</strong><span>Required HMM Prep Steps</span></div>
</div>

## Method At A Glance

<div class="overview-grid">
  <div class="paper-card">
    <h3>1 · Discovery</h3>
    <p>Start from protein FASTAs or transcriptomes. Transcriptome inputs are converted to protein candidates with Pyrodigal.</p>
  </div>
  <div class="paper-card">
    <h3>2 · Filtering</h3>
    <p>Apply coverage and minimum-length filters, collapse near-duplicates, and retain reference-like candidates with traceable logging.</p>
  </div>
  <div class="paper-card">
    <h3>3 · Classification</h3>
    <p>Place candidates in the multi-clade TPS feature space, transfer labels from nearest references, and optionally run ESM2 CeeSs scoring.</p>
  </div>
  <div class="paper-card">
    <h3>4 · Visualization</h3>
    <p>Render PCA/LDA views and supporting tables for candidate triage. No <code>04_phylogeny/</code> directory is part of the main workflow.</p>
  </div>
</div>

## Quick Start

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --skip-phylogeny \
  --output-dir results/
```

This command will:

- discover and filter candidate TPS proteins;
- classify retained candidates against the supplied reference FASTA directory;
- write PCA/LDA visualization artifacts under `03_classification/`;
- add CeeSs scoring when `TPS.xlsx` and the `[esm]` extra are available.

## Where To Go Next

- **[Getting Started](getting-started.md)** - installation and your first run.
- **[Method](method.md)** - a stage-by-stage explanation of the current workflow.
- **[Tutorials](tutorials.md)** - practical command sequences and analysis pathways.
- **[CLI Reference](cli-reference.md)** - main workflow commands and parameters.
- **[Outputs](outputs.md)** - how to read every main artifact.
