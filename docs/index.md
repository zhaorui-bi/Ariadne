# Ariadne

<section class="research-hero" markdown>
<div markdown>
<span class="research-hero__eyebrow">TPS Discovery Framework</span>

**Ariadne** is a four-stage computational framework for large-scale terpene synthase discovery, quality-controlled candidate reduction, subtype classification, and downstream CeeSs prioritization.

<div class="research-hero__subtitle" markdown>
The documentation follows the project model used in the OUC-HWU guide: concise navigation, research-oriented pages, stable figures, and a restrained Times New Roman visual system.
</div>

<div class="research-hero__meta" markdown>
**Scope.** Transcriptome or protein FASTA input -> discovery -> filtering -> feature-space classification -> PCA/LDA visualization and optional ESM2 CeeSs scoring.
</div>

<div class="research-hero__actions" markdown>
[Get Started](getting-started.md){ .md-button .md-button--primary }
[Read Method](method.md){ .md-button }
[Advanced Usage](advanced-usage.md){ .md-button }
</div>
</div>

<div class="research-hero__panel" markdown>
**Documentation Status**

- Version baseline: `ariadne-tps 1.1.0`
- Updated: 2026-07-09
- Maintainer: Zhaorui Jiang
- Primary outputs: TSV evidence tables and SVG projections
- Repository: [zhaorui-bi/Ariadne](https://github.com/zhaorui-bi/Ariadne)
</div>
</section>

<div class="metric-strip" markdown>
<div markdown>
**4**
<span>Analysis stages</span>
</div>
<div markdown>
**ESM2**
<span>Optional scoring</span>
</div>
<div markdown>
**PCA/LDA**
<span>Projection layer</span>
</div>
<div markdown>
**TSV/SVG**
<span>Auditable output</span>
</div>
</div>

## Algorithm Framework

<figure class="algorithm-figure">
  <img src="images/algorithm-framework.png" alt="Ariadne algorithm framework with discovery, filtering, classification, and visualization stages">
  <figcaption><strong>Figure 1.</strong> Ariadne starts from transcriptomic or protein resources, screens sequence hits with profile-HMM evidence, filters candidates by quality and redundancy, extracts ESM2 and profile-space representations, and prioritizes subtypes for visualization and wet-lab follow-up.</figcaption>
</figure>

## Workflow At A Glance

<div class="workflow-map workflow-map--four" markdown>

<div class="workflow-step" markdown>
<span class="stage-label">Stage I</span>
### Discovery
Search transcriptome-derived ORFs or predicted protein FASTAs with a query HMM. The result is a candidate sequence universe with traceable hit evidence.

Primary artifact: `candidates.protein.faa`
</div>

<div class="workflow-step" markdown>
<span class="stage-label">Stage II</span>
### Filtering
Apply coverage, minimum-length, and near-duplicate controls to build a compact, non-redundant candidate set.

Primary artifact: `candidates.filtered.faa`
</div>

<div class="workflow-step" markdown>
<span class="stage-label">Stage III</span>
### Classification
Represent candidates and references in the same TPS feature space, then report nearest-reference evidence and optional ESM2 CeeSs probabilities.

Primary artifact: `classification.tsv`
</div>

<div class="workflow-step" markdown>
<span class="stage-label">Stage IV</span>
### Visualization
Project the feature space with LDA or PCA so candidate placement can be inspected before experimental validation.

Primary artifact: `embedding.svg`
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

This run produces:

- discovery and hit evidence under `01_discovery/`;
- filtered, non-redundant sequences under `02_filtering/`;
- subtype assignments, nearest-neighbor evidence, ESM2 scoring outputs, and projections under `03_classification/`;
- a machine-readable `pipeline_summary.tsv` at the output root.

## Documentation Map

<div class="evidence-grid" markdown>

<div class="evidence-card" markdown>
### [Getting Started](getting-started.md)
Install Ariadne, prepare a reference FASTA directory, and run the first complete workflow.
</div>

<div class="evidence-card" markdown>
### [Method](method.md)
Read the staged algorithmic design behind discovery, filtering, classification, and visualization.
</div>

<div class="evidence-card" markdown>
### [Advanced Usage](advanced-usage.md)
Tune thresholds, switch CeeSs classifier heads, reuse checkpoints, and design manuscript-grade runs.
</div>

<div class="evidence-card" markdown>
### [Software Architecture](software-architecture.md)
Review package structure, public API contracts, dependency boundaries, and release checks.
</div>

<div class="evidence-card" markdown>
### [Outputs](outputs.md)
Interpret every table and figure emitted by the pipeline, including CeeSs and embedding artifacts.
</div>

</div>

!!! note "Current workflow boundary"
    Ariadne prioritizes computational triage and visualization. Wet-lab validation is the downstream experimental step informed by the software outputs, not a step executed by the package.
