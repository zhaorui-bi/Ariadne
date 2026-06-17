# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.1.0] - 2026-06-17

### Added
- Top-level Python API: the main pipeline functions (`discover_candidates`,
  `filter_candidates`, `classify_candidates`, `build_phylogeny`, …) are now
  re-exported from `ariadne` and resolved lazily (PEP 562), so `import ariadne`
  stays light while `import ariadne as ad; ad.classify_candidates(...)` works
  from any script or notebook.
- `examples/tutorial.ipynb`: a runnable, end-to-end notebook walkthrough that
  uses the bundled HMMs and renders the embedding inline (no external tools).
- Test suite (`tests/`) covering FASTA/TSV I/O, sequence utilities, filtering,
  reference preparation, Newick parsing, feature-space math, classification, the
  top-level API, and the CLI, runnable without PyTorch/MAFFT/IQ-TREE.
- GitHub Actions CI running `ruff` and `pytest` across Python 3.9–3.12.
- `--log-file` global CLI flag to mirror logs to a file (without ANSI colours).
- `py.typed` marker so downstream type checkers see Ariadne as typed.
- Packaging metadata in `pyproject.toml` (authors, URLs, classifiers, keywords)
  and `test`/`dev` optional-dependency groups.
- `CONTRIBUTING.md`, `CITATION.cff`, issue/PR templates, and `examples/`.

### Changed
- `ariadne run`/`classify`/`phylogeny` help is now organised into argument groups
  (classification, CeeSs ESM scoring, CeeSs advanced tuning, phylogeny), and the
  shared flag definitions were de-duplicated into helpers. No flags changed.
- `write_fasta` and `write_tsv` now create missing parent directories.
- `install.sh` defaults to the standard PyPI index (the previous mirror is still
  available via `PIP_INDEX_URL`).

### Fixed
- Corrected the package docstring to reference the real module names
  (`search.py`, `embed.py`).
- Replaced the broken `code.sh` (which called a non-existent `prepare-demo`
  command) with a working `examples/run_example.sh`.
- Removed a duplicated mid-file `import math` in `model.py`.
- `model.py`: removed a duplicate, non-returning `_save_mlp_classifier_checkpoint`
  definition that shadowed the real one, causing the saved MLP checkpoint path to
  be dropped from the returned `ceess_classifier_checkpoint` output.
- `data.py`: the insect-workbook parser used `-1` as a "missing column" sentinel,
  which is a valid index and silently read the last column when `Accession`,
  `Species`, or `Clade` were absent; it now falls back to the correct placeholders.
- `tree.py`: IQ-TREE no longer receives `--fast` together with `-B` (ultrafast
  bootstrap), a combination IQ-TREE 2 rejects; `--fast` is dropped when bootstrap
  replicates are requested.

### Removed
- Generated artifacts and OS cruft (`test_run/`, `output/`, `.DS_Store`) are no
  longer tracked in version control.

## [1.0.0] - 2026-04-07

Initial public release of the four-stage Ariadne platform: HMM-guided discovery,
coverage/length/identity filtering, TPS feature-space classification with an
optional ESM2 CeeSs scoring layer, and MAFFT + IQ-TREE phylogeny.
