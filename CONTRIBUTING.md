# Contributing to Ariadne

Thanks for your interest in improving Ariadne! This document explains how to set
up a development environment, run the checks, and submit changes.

## Development setup

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne

# Create and activate a virtual environment (Python 3.9+)
python -m venv .venv
source .venv/bin/activate

# Install the package together with the development tooling
pip install -e ".[dev]"
```

The optional ESM2 CeeSs scoring stage and the phylogeny stage require extra
software:

- ESM2 scoring: `pip install -e ".[esm]"` (PyTorch + Transformers).
- Phylogeny: `mafft` and `iqtree`/`iqtree2` on your `PATH`.

These are **not** required to run the test suite.

## Running the checks

```bash
# Lint
ruff check ariadne tests

# Tests (with coverage)
pytest --cov=ariadne --cov-report=term-missing
```

Tests that depend on optional components (PyTorch, MAFFT, IQ-TREE) are skipped
automatically when those components are unavailable, so the suite runs anywhere.

## Pull request guidelines

1. Create a topic branch off `main`.
2. Keep changes focused; one logical change per pull request.
3. Add or update tests for any behaviour you change.
4. Make sure `ruff check` and `pytest` pass locally.
5. Update `CHANGELOG.md` under the `[Unreleased]` heading.
6. Follow the existing code style: type hints, descriptive docstrings, and the
   project's `typing.Optional`/`typing.Union` conventions (the codebase targets
   Python 3.9).

## Reporting bugs and requesting features

Please use the GitHub issue templates. Include the Ariadne version
(`ariadne --version`), your platform, and a minimal reproduction when possible.

## Code of conduct

Be respectful and constructive. Harassment or discrimination of any kind is not
tolerated.
