# Documentation Site

The documentation is built with MkDocs + Material.

## Local Preview

```bash
python -m pip install -r docs/requirements.txt
mkdocs serve
```

## Current Content Model

The site follows the README workflow:

- users provide reference FASTA directories with `--reference-dir`;
- no separate HMM generation step is documented for users;
- Stage 4 is PCA/LDA visualization, not a final tree section;
- generated outputs are described under `03_classification/`.

## Asset Policy

The current logo lives at `docs/fig/logo.png`. Do not reference the old `docs/assets/` SVG files unless they are restored.
