# Documentation Site

The documentation is built with MkDocs + Material.

## Local Preview

```bash
python -m pip install -r docs/requirements.txt
mkdocs serve
```

## Current Content Model

The site is framed around the *Journal of Natural Products* (2026) article and a restrained academic visual system:

- Times New Roman is the primary typography;
- the landing page leads with the published scientific context, then the four computational stages;
- line-art icons live under `docs/images/icons/` and are reused by the README;
- users provide reference FASTA directories with `--reference-dir`;
- Ariadne can build required HMM resources from references when prebuilt profiles are not supplied;
- Stage 4 is PCA/LDA visualization and candidate triage, not wet-lab validation;
- generated outputs are described under `03_classification/`.

## Asset Policy

Stable documentation images should live under `docs/images/` or next to the page that uses them. The current logo lives at `docs/logo.png`; the algorithm framework figure lives at `docs/images/algorithm-framework.png`; the LDA example (Figure S2) lives at `docs/images/published_embedding_3d_sections.svg`; line-art icons live at `docs/images/icons/`. Working copies, recovery scripts, and unpublished drafts stay in the local `fig/` directory and are not committed. Do not reference temporary chat or WeChat image paths in committed docs.
