<p align="center">
  <img src="docs/logo.png" alt="Ariadne, a terpene synthase discovery platform" width="380">
</p>

<p align="center">
  <strong>Genome-wide targeted mining of coral terpene synthases</strong><br>
  <em>Journal of Natural Products</em> (2026) &nbsp;·&nbsp; DOI: <a href="https://doi.org/10.1021/acs.jnatprod.6c00686">10.1021/acs.jnatprod.6c00686</a>
</p>

<p align="center">
  <a href="./README_ZH.md">中文</a> &nbsp;·&nbsp;
  <a href="https://ariadne-platform.readthedocs.io/en/latest/">Read the Docs</a> &nbsp;·&nbsp;
  <a href="https://doi.org/10.1021/acs.jnatprod.6c00686">Paper</a> &nbsp;·&nbsp;
  <a href="#citation">Citation</a>
</p>

<p align="center">
  <a href="https://ariadne-platform.readthedocs.io/en/latest/"><img alt="Read the Docs" src="https://img.shields.io/readthedocs/ariadne-platform/latest?style=flat-square&label=Read%20the%20Docs&color=1d4ed8"></a>
  <a href="https://doi.org/10.1021/acs.jnatprod.6c00686"><img alt="Journal of Natural Products 2026" src="https://img.shields.io/badge/J.%20Nat.%20Prod.-2026-084c61?style=flat-square"></a>
  <a href="https://doi.org/10.1021/acs.jnatprod.6c00686"><img alt="DOI" src="https://img.shields.io/badge/DOI-10.1021%2Facs.jnatprod.6c00686-084c61?style=flat-square"></a>
  <img alt="Python" src="https://img.shields.io/badge/Python-%E2%89%A53.9-0b132b?style=flat-square&logo=python&logoColor=white">
  <img alt="Version" src="https://img.shields.io/badge/Version-1.1.0-0f766e?style=flat-square">
  <img alt="License" src="https://img.shields.io/badge/License-MIT-c2410c?style=flat-square">
</p>

---

**Ariadne** is the computational platform described in Wang et al., *J. Nat. Prod.* **2026**, for genome-wide targeted mining of coral terpene synthases (TPS) and prioritization of **cembrene-class synthases (CeeSs)**. It converts transcriptome assemblies or predicted protein FASTAs into an auditable four-stage workflow: discovery, filtering, classification, and PCA/LDA visualization, with optional ESM2 scoring of coral-like CeeSs.

> Genome mining is an efficient route to natural-product enzymes, yet TPS searches often rediscover the same products. Ariadne was built to target the 14-membered carbocyclic cembranoid pathway in corals and to surface experimentally tractable CeeSs rather than a long list of unprioritized homologs.

## Published findings

Results from the accompanying article, not claims about a local software run:

<table>
  <tr>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/flask.svg" width="28" height="28" alt=""><br>
      <strong>5 CeeSs</strong><br>
      cembrene synthases identified from coral TPS homologs
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/target.svg" width="28" height="28" alt=""><br>
      <strong>80% accuracy</strong><br>
      experimental validation rate across 90 TPS homologs
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/ring.svg" width="28" height="28" alt=""><br>
      <strong>14-membered rings</strong><br>
      cembranoid scaffolds prioritized for biosynthesis
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/genome.svg" width="28" height="28" alt=""><br>
      <strong>Polyphyletic origins</strong><br>
      evolutionary trajectory used to guide enzyme engineering
    </td>
  </tr>
</table>

The published set includes two enzymes with high similarity to known cembrene B synthases and three low-similarity enzymes confirmed by heterologous expression in yeast. Phylogeny and subsequent engineering of an early-diverging coral enzyme, which produced an unreported cembrene scaffold, are scientific analyses enabled by these candidates. They are not additional CLI stages in this repository.

## Scientific context

Cembranoids are coral-derived diterpenoids built on a 14-membered carbocyclic skeleton. Chemical synthesis of these scaffolds remains difficult, and only two biosynthetic precursors had been obtained from fifteen previously reported coral enzymes across six studies. Conventional TPS genome mining tends to recover the same product families. Ariadne addresses that bottleneck by combining profile-HMM discovery, quality filtering, reference-space classification, and CeeSs-focused visualization so that wet-lab effort can be spent on novel or poorly characterized coral enzymes.

## Platform

<table>
  <tr>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/discovery.svg" width="28" height="28" alt=""><br>
      <strong>I. Discovery</strong><br>
      HMM search of transcriptome ORFs or protein FASTAs
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/filter.svg" width="28" height="28" alt=""><br>
      <strong>II. Filtering</strong><br>
      coverage, length, and near-duplicate controls
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/classify.svg" width="28" height="28" alt=""><br>
      <strong>III. Classification</strong><br>
      reference-space labels with nearest-neighbor evidence
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/visualize.svg" width="28" height="28" alt=""><br>
      <strong>IV. Visualization</strong><br>
      PCA/LDA projections for candidate triage
    </td>
  </tr>
</table>

<p align="center">
  <img src="docs/images/algorithm-framework.png" alt="Ariadne four-stage algorithm framework" width="880">
</p>

## Capabilities

<table>
  <tr>
    <td valign="top" width="50%">
      <img src="docs/images/icons/genome.svg" width="22" height="22" alt="">
      <strong> Targeted coral TPS mining</strong><br>
      Genome- or transcriptome-wide screening against a user-supplied reference FASTA directory, with HMM resources built when prebuilt profiles are absent.
    </td>
    <td valign="top" width="50%">
      <img src="docs/images/icons/output.svg" width="22" height="22" alt="">
      <strong> Auditable intermediates</strong><br>
      Every stage writes TSV evidence so candidates can be kept, discarded, or re-ranked without a black-box endpoint.
    </td>
  </tr>
  <tr>
    <td valign="top">
      <img src="docs/images/icons/target.svg" width="22" height="22" alt="">
      <strong> CeeSs prioritization</strong><br>
      Optional ESM2 scoring of coral-like candidates when <code>TPS.xlsx</code> and the <code>[esm]</code> extras are available.
    </td>
    <td valign="top">
      <img src="docs/images/icons/visualize.svg" width="22" height="22" alt="">
      <strong> Manuscript-style figures</strong><br>
      Three-panel LDA/PCA SVGs (LD1–LD2, LD1–LD3, LD2–LD3) for inspection before heterologous expression.
    </td>
  </tr>
</table>

## Installation

Recommended:

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

Minimal:

```bash
pip install -e .
```

With ESM2 CeeSs scoring:

```bash
pip install -e '.[esm]'
```

## Quick Start

Prepare a reference FASTA directory, then run:

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir results/
```

Transcriptome mode predicts ORFs with Pyrodigal before discovery:

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir reference_fastas/ \
  --output-dir results_from_transcriptomes/
```

Core output layout:

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
└── pipeline_summary.tsv
```

Wet-lab validation remains a downstream experimental step. Ariadne reports computational evidence; it does not execute enzyme assays.

## Documentation

- [Read the Docs](https://ariadne-platform.readthedocs.io/en/latest/)
- [Getting Started](./docs/getting-started.md)
- [Method](./docs/method.md)
- [CeeSs Classifier](./docs/esm-type.md)
- [Outputs](./docs/outputs.md)
- [CLI Reference](./docs/cli-reference.md)
- [Citation](./docs/citation.md)

## Development

```bash
pip install -e '.[dev]'
python -m ariadne --help
mkdocs build --strict
```

## Citation

If you use Ariadne in your work, please cite:

Wang, Y.; Yu, M.; Jiang, Z.; Zhou, C.; Feng, W.; Yu, K.; Ju, J.; Li, F. Identification of Cembrene Synthases Uncovers Polyphyletic Origins of 14-Membered Carbocyclic Cembranoids Biosynthesis in Corals. *J. Nat. Prod.* **2026**. DOI: [10.1021/acs.jnatprod.6c00686](https://doi.org/10.1021/acs.jnatprod.6c00686)

```bibtex
@article{wang2026cembrene,
  author    = {Wang, Yinhao and Yu, Mengmeng and Jiang, Zhaorui and Zhou, Chengyu and Feng, Wei and Yu, Ke and Ju, Jianhua and Li, Feng},
  title     = {Identification of Cembrene Synthases Uncovers Polyphyletic Origins of 14-Membered Carbocyclic Cembranoids Biosynthesis in Corals},
  journal   = {Journal of Natural Products},
  year      = {2026},
  month     = {sep},
  publisher = {American Chemical Society},
  doi       = {10.1021/acs.jnatprod.6c00686},
  url       = {https://doi.org/10.1021/acs.jnatprod.6c00686},
  note      = {Published online September 8, 2026}
}
```

## License

Released under the [MIT License](./LICENSE).
