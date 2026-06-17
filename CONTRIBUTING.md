# Contributing to Ariadne

Thank you for your interest in improving Ariadne. This guide covers how to set up
a development environment, run the static checks, and submit changes.

## Development setup

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne

# Create and activate a virtual environment (Python 3.9+)
python -m venv .venv
source .venv/bin/activate

# Install the package together with the development tooling
pip install -e ".[dev]"
```

Two stages depend on additional software that is **not** required for everyday
development:

- **CeeSs scoring** — `pip install -e ".[esm]"` (PyTorch + Transformers).
- **Phylogeny** — `mafft` and `iqtree`/`iqtree2` on your `PATH`.

## Static checks

```bash
ruff check ariadne
```

The maintained unit-test suite runs against the optional, lazily-imported
components and is exercised by the maintainers before each release. If your
change touches behaviour, please describe how you verified it in the pull
request (a minimal reproduction or an example command is ideal).

## Pull request guidelines

1. Create a topic branch off `main`.
2. Keep changes focused — one logical change per pull request.
3. Make sure `ruff check ariadne` passes locally.
4. Update `CHANGELOG.md` under the `[Unreleased]` heading.
5. Follow the existing style: type hints, descriptive docstrings, and the
   project's `typing.Optional` / `typing.Union` conventions (the codebase
   targets Python 3.9).

## Reporting bugs and requesting features

Please open a [GitHub issue](https://github.com/zhaorui-bi/Ariadne/issues) and
include the Ariadne version (`ariadne --version`), your platform, and a minimal
reproduction whenever possible.

## Code of conduct

Be respectful and constructive. Harassment or discrimination of any kind is not
tolerated.
