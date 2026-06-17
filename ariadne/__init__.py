"""Ariadne — coral terpene synthase (TPS) discovery, CeeSs prioritization & phylogeny.

Ariadne ships both a command-line interface (``ariadne ...`` / ``python -m
ariadne``) and an importable Python API. The most useful pipeline functions are
re-exported at the package top level so the platform can be driven directly from
a script or a notebook::

    import ariadne as ad

    # 1. discovery
    candidates = ad.discover_candidates_from_proteins(
        protein_files=ad.collect_protein_files("input"),
        query_hmm="tps.hmm",
        output_dir="results/01_discovery",
    )

    # 2. filtering
    ad.filter_candidates("results/01_discovery/candidates.protein.faa",
                         reference_dir="tree", output_dir="results/02_filtering")

    # 3. classification (HMM profile feature space)
    ad.classify_candidates("results/02_filtering/candidates.filtered.faa",
                           reference_dir="tree", output_dir="results/03_classification")

Heavy, optional dependencies (PyTorch, pyhmmer, pyrodigal) are imported lazily,
so a bare ``import ariadne`` always succeeds even in a minimal environment. The
underlying submodule is imported only the first time you access a function that
needs it (PEP 562). The names available on the package are listed in ``__all__``.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

__version__     = "1.1.0"
__author__      = "Zhaorui Jiang"
__description__ = "Coral TPS Discovery, Candidate CeeSs Prioritization & Phylogeny Platform"
__url__         = "https://github.com/zhaoruijiang26/Ariadne"

# Public API: exported name -> submodule that defines it. Resolved lazily on
# first access (see ``__getattr__``) so importing ``ariadne`` never eagerly pulls
# in torch / pyhmmer / pyrodigal. Grouped by pipeline stage for readability.
_LAZY_EXPORTS = {
    # --- core data types & helpers (ariadne.utils) ---------------------------
    "FastaRecord": "utils",
    "read_fasta": "utils",
    "write_fasta": "utils",
    "setup_logging": "utils",
    # --- reference preparation (ariadne.data) --------------------------------
    "prepare_coral_reference": "data",
    "prepare_insect_reference": "data",
    "prepare_extra_reference": "data",
    "load_reference_records": "data",
    "write_reference_metadata": "data",
    # --- discovery (ariadne.search) ------------------------------------------
    "build_hmm": "search",
    "search_proteins_with_hmm": "search",
    "collect_protein_files": "search",
    "discover_candidates": "search",
    "discover_candidates_from_proteins": "search",
    # --- filtering (ariadne.filter) ------------------------------------------
    "filter_candidates": "filter",
    "filter_by_coverage": "filter",
    "filter_by_length": "filter",
    "deduplicate_exact": "filter",
    "near_duplicate": "filter",
    # --- classification in HMM profile space (ariadne.embed) -----------------
    "classify_candidates": "embed",
    # --- ESM2 / CeeSs modelling (ariadne.model) ------------------------------
    "resolve_esm_model_name": "model",
    "load_tps_xlsx": "model",
    "compute_esm_embeddings": "model",
    "analyze_tps_types_with_esm": "model",
    "classify_ceess_candidates_with_esm": "model",
    "classify_ceess_candidates_with_supcon": "model",
    # --- phylogeny (ariadne.tree) --------------------------------------------
    "prepare_phylogeny_input": "tree",
    "run_mafft": "tree",
    "run_iqtree": "tree",
    "build_phylogeny": "tree",
    "render_phylogeny_preview": "tree",
}

__all__ = [
    "__version__",
    "__author__",
    "__description__",
    "__url__",
    *sorted(_LAZY_EXPORTS),
]


def __getattr__(name: str):
    """Lazily import and cache a re-exported pipeline function (PEP 562)."""
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = importlib.import_module(f"{__name__}.{module_name}")
    value = getattr(module, name)
    globals()[name] = value  # cache so subsequent access skips __getattr__
    return value


def __dir__():
    return sorted(__all__)


if TYPE_CHECKING:  # pragma: no cover - aids IDEs / type checkers only
    from ariadne.data import (
        load_reference_records,
        prepare_coral_reference,
        prepare_extra_reference,
        prepare_insect_reference,
        write_reference_metadata,
    )
    from ariadne.embed import classify_candidates
    from ariadne.filter import (
        deduplicate_exact,
        filter_by_coverage,
        filter_by_length,
        filter_candidates,
        near_duplicate,
    )
    from ariadne.model import (
        analyze_tps_types_with_esm,
        classify_ceess_candidates_with_esm,
        classify_ceess_candidates_with_supcon,
        compute_esm_embeddings,
        load_tps_xlsx,
        resolve_esm_model_name,
    )
    from ariadne.search import (
        build_hmm,
        collect_protein_files,
        discover_candidates,
        discover_candidates_from_proteins,
        search_proteins_with_hmm,
    )
    from ariadne.tree import (
        build_phylogeny,
        prepare_phylogeny_input,
        render_phylogeny_preview,
        run_iqtree,
        run_mafft,
    )
    from ariadne.utils import FastaRecord, read_fasta, setup_logging, write_fasta
