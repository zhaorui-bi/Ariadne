# Documentation Site

The documentation is built with MkDocs + Material.

## Local Preview

```bash
python -m pip install -r docs/requirements.txt
mkdocs serve
```

## Current Content Model

The site follows a research-documentation model inspired by the OUC-HWU guide:

- Times New Roman is the primary typography for the documentation site;
- the landing page is a method-oriented project overview rather than a marketing page;
- users provide reference FASTA directories with `--reference-dir`;
- Ariadne can build required HMM resources from references when prebuilt profiles are not supplied;
- Stage 4 is PCA/LDA visualization and candidate triage, not the wet-lab validation itself;
- generated outputs are described under `03_classification/`.

## Asset Policy

Stable documentation images should live under `docs/images/` or next to the page that uses them. The current logo lives at `docs/logo.png`; the algorithm framework figure lives at `docs/images/algorithm-framework.png`. Do not reference temporary chat or WeChat image paths in committed docs.
