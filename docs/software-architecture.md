# Software Architecture

Ariadne is structured as a reusable Python package first and a command-line tool second. The CLI is a thin layer over the same stage functions that are available from the Python API.

## Package Layout

| Module | Responsibility |
| --- | --- |
| `ariadne.search` | Stage I discovery: HMM construction, transcriptome ORF prediction, and protein-mode HMM search |
| `ariadne.filter` | Stage II filtering: coverage/length QC, near-duplicate clustering, and review tables |
| `ariadne.embed` | Stages III-IV: profile-space classification, nearest-reference evidence, PCA/LDA visualization |
| `ariadne.model` | Optional ESM2 CeeSs scoring and supervised TPS type diagnostics |
| `ariadne.data` | Reference FASTA preparation and metadata loading |
| `ariadne.utils` | Logging, FASTA/TSV I/O, and shared sequence helpers |
| `ariadne.cli` | CLI argument parsing and stage orchestration |

## Public API Contract

The package top level re-exports the stable functions most users need:

```python
import ariadne as ad

protein_files = ad.collect_protein_files("my_proteins")
discovery = ad.discover_candidates_from_proteins(
    protein_files,
    "query.hmm",
    "results/01_discovery",
)
filtering = ad.filter_candidates(
    discovery["candidate_proteins"],
    "results/02_filtering",
    reference_dir="reference_fastas",
)
classification = ad.classify_candidates(
    filtering["filtered_fasta"],
    "reference_fastas",
    "results/03_classification",
    hmm_dir="tps_hmms",
    ceess_xlsx="TPS.xlsx",
)
```

Public functions return dictionaries of concrete output paths. This is intentional: every stage writes inspectable artifacts that can be archived, compared, or rerun independently.

## Dependency Strategy

Ariadne keeps the core installation small:

- core bioinformatics and table dependencies are declared in `pyproject.toml`;
- ESM2 dependencies live behind the `[esm]` optional extra;
- package-level imports are lazy where possible, so `import ariadne` remains lightweight.

Install modes:

```bash
pip install -e .
pip install -e '.[esm]'
pip install -e '.[dev]'
```

## Output Contract

The release workflow writes durable TSV/SVG artifacts rather than transient objects:

| Stage | Contract files |
| --- | --- |
| Discovery | `candidates.protein.faa`, `candidates.hits.tsv` |
| Filtering | `candidates.filtered.faa`, `filter_report.tsv`, `dedupe_clusters.tsv` |
| Classification | `classification.tsv`, `nearest_neighbors.tsv`, `assignment_summary.tsv` |
| Visualization | `embedding.tsv`, `embedding.svg`, `embedding_3d_sections.svg` |
| CeeSs scoring | `ceess_predictions.tsv`, `ceess_candidates.tsv`, `ceess_candidates.fasta` |

Downstream notebooks and manuscript scripts should consume these files rather than importing private helper functions.

## Release Checklist

Before tagging a release:

```bash
python -m pip install -e '.[dev]'
python -m pytest
python -m ruff check ariadne
mkdocs build --strict
python -m build
```

Also verify:

- `README.md` and `README_ZH.md` describe the same workflow boundary;
- `mkdocs.yml` navigation includes new user-facing pages;
- every committed image is under `docs/images/` or another tracked docs path (`fig/` is local-only);
- optional ESM2 examples clearly mention `pip install -e '.[esm]'`;
- release notes document any output-schema changes.
