"""Command-line interface for the Ariadne TPS discovery pipeline.

The CLI is organised around the current four-stage workflow:
discovery -> filtering -> classification -> phylogeny.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Tuple, Union

from ariadne import __version__
from ariadne import utils as _log

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


def _repo_root() -> Path:
    """Return the repository root, used for bundled resource lookup."""
    return Path(__file__).resolve().parent.parent


def _first_existing_path(*candidates: PathLike) -> Optional[Path]:
    """Return the first existing path from a prioritized list of candidates."""
    for candidate in candidates:
        if candidate is None:
            continue
        path = Path(candidate).expanduser().resolve()
        if path.exists():
            return path
    return None


def _default_coral() -> Optional[Path]:
    """Locate the default coral reference FASTA when present."""
    return _first_existing_path(
        Path.cwd() / "coralTPS (modified)-cembrene.fasta",
        _repo_root() / "coralTPS (modified)-cembrene.fasta",
        Path(__file__).resolve().parent / "coralTPS (modified)-cembrene.fasta",
    )


def _default_insect() -> Optional[Path]:
    """Locate the default insect workbook when present."""
    return _first_existing_path(
        Path.cwd() / "Insecta TPS.xlsx",
        _repo_root() / "Insecta TPS.xlsx",
        Path(__file__).resolve().parent / "Insecta TPS.xlsx",
    )


def _default_bacteria() -> Optional[Path]:
    """Locate an optional bacterial reference FASTA."""
    return _first_existing_path(
        Path.cwd() / "bacteria.fasta",
        Path.cwd() / "bacteria.fa",
        Path.cwd() / "bacteria.faa",
        _repo_root() / "bacteria.fasta",
        _repo_root() / "bacteria.fa",
        _repo_root() / "bacteria.faa",
    )


def _default_fungal() -> Optional[Path]:
    """Locate an optional fungal reference FASTA."""
    return _first_existing_path(
        Path.cwd() / "fungal.fasta",
        Path.cwd() / "fungal.fa",
        Path.cwd() / "fungal.faa",
        Path.cwd() / "fungi.fasta",
        Path.cwd() / "fungi.fa",
        Path.cwd() / "fungi.faa",
        _repo_root() / "fungal.fasta",
        _repo_root() / "fungal.fa",
        _repo_root() / "fungal.faa",
        _repo_root() / "fungi.fasta",
        _repo_root() / "fungi.fa",
        _repo_root() / "fungi.faa",
    )


def _default_plant() -> Optional[Path]:
    """Locate an optional plant reference FASTA."""
    return _first_existing_path(
        Path.cwd() / "plant.fasta",
        Path.cwd() / "plant.fa",
        Path.cwd() / "plant.faa",
        _repo_root() / "plant.fasta",
        _repo_root() / "plant.fa",
        _repo_root() / "plant.faa",
    )


def _default_reference_dir() -> Optional[Path]:
    """Locate the default tree/reference directory when present."""
    return _first_existing_path(
        Path.cwd() / "tree",
        _repo_root() / "tree",
    )


def _default_tps_xlsx() -> Optional[Path]:
    """Locate the coral TPS spreadsheet used by the integrated CeeSs workflow."""
    return _first_existing_path(
        Path.cwd() / "TPS" / "TPS.xlsx",
        _repo_root() / "TPS" / "TPS.xlsx",
    )


def _default_reference_alignment(reference_dir: Optional[PathLike] = None) -> Optional[Path]:
    """Locate the default coral reference FASTA under ``tree/`` when available."""
    candidates: list[Path] = []
    if reference_dir is not None:
        directory = Path(reference_dir)
        candidates.extend(
            [
                directory / "coral.fasta",
                directory / "coral.fa",
                directory / "coral.faa",
            ]
        )
    default_reference_dir = _default_reference_dir()
    if default_reference_dir is not None:
        candidates.extend(
            [
                default_reference_dir / "coral.fasta",
                default_reference_dir / "coral.fa",
                default_reference_dir / "coral.faa",
            ]
        )
    return _first_existing_path(*candidates)


def _find_reference_alignment(reference_dir: PathLike) -> Path:
    """Infer a usable coral reference FASTA from the reference directory."""
    default_alignment = _default_reference_alignment(reference_dir)
    if default_alignment is not None and default_alignment.exists():
        return default_alignment
    directory = Path(reference_dir)
    for filename in ("coral.fasta", "coral.fa", "coral.faa"):
        direct = directory / filename
        if direct.exists():
            return direct
    for fasta_path in sorted(directory.glob("*.fa*")):
        lowered = fasta_path.name.lower()
        if "coral" in lowered:
            return fasta_path
    raise FileNotFoundError(
        f"Could not locate a coral reference FASTA under {directory}. "
        "Expected something like coral.fasta."
    )


def _reference_fasta_paths(reference_dir: PathLike) -> list[Path]:
    """Collect FASTA files from a prepared reference directory."""
    directory = Path(reference_dir)
    fasta_paths = sorted(directory.glob("*.fa*"))
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files were found in reference directory: {directory}")
    return fasta_paths


def _bundled_hmm_dir() -> Path:
    """Return the bundled HMM directory shipped with the Ariadne package."""
    package_hmm = Path(__file__).resolve().parent / "hmm"
    if package_hmm.exists():
        return package_hmm
    # Backward-compatible fallback for older worktrees that still used repo-root hmm/.
    return Path(__file__).resolve().parent.parent / "hmm"


def _bundled_query_hmm() -> Path | None:
    """Return the bundled default discovery HMM when it is available."""
    path = _bundled_hmm_dir() / "query.hmm"
    return path if path.exists() else None


def _bundled_tps_hmm_dir() -> Path | None:
    """Return the bundled TPS HMM library directory when it is available."""
    directory = _bundled_hmm_dir()
    if not directory.exists():
        return None
    hmm_files = sorted(path for path in directory.glob("*.hmm") if path.name != "query.hmm")
    return directory if hmm_files else None


def _auto_build_query_hmm(reference_dir: PathLike, output_path: PathLike, *, name: str) -> Path:
    """Build the discovery query HMM from the default reference FASTA in ``reference_dir``."""
    from ariadne.search import build_hmm

    source_fasta = _find_reference_alignment(reference_dir)
    logger.info("No --query-hmm provided; building query HMM from %s", source_fasta)
    return build_hmm(source_fasta, output_path, name=name)


def _auto_build_tps_hmm_library(reference_dir: PathLike, output_dir: PathLike) -> Path:
    """Build a TPS HMM library from all FASTA files in ``reference_dir``."""
    from ariadne.search import build_hmm
    from ariadne.utils import ensure_directory

    destination = ensure_directory(output_dir)
    built_any = False
    for fasta_path in _reference_fasta_paths(reference_dir):
        output_path = destination / f"{fasta_path.stem}.hmm"
        build_hmm(fasta_path, output_path, name=fasta_path.stem)
        built_any = True
    if not built_any:
        raise FileNotFoundError(f"Could not build any HMMs from reference directory: {reference_dir}")
    logger.info("Built TPS HMM library from %s -> %s", reference_dir, destination)
    return destination


def _parse_extra_reference(spec: str) -> Tuple[str, str]:
    """Parse ``SOURCE=PATH`` syntax used for extra references."""
    if "=" not in spec:
        raise ValueError(f"Expected SOURCE=PATH for extra references, got: {spec}")
    source, path = spec.split("=", 1)
    return source.strip(), path.strip()


def _existing_path_or_none(value: Optional[PathLike], *, label: str) -> Optional[Path]:
    """Return an existing path or log a warning when the file is missing."""
    if value is None:
        return None
    path = Path(value).expanduser()
    if path.exists():
        return path
    logger.warning("%s file not found, skipped: %s", label, path)
    return None


def cmd_prepare_references(args: argparse.Namespace) -> int:
    """Prepare reference FASTA files and metadata tables."""
    from ariadne.utils import ensure_directory
    from ariadne.data import (
        prepare_coral_reference,
        prepare_extra_reference,
        prepare_insect_reference,
        write_reference_metadata,
    )

    output_dir = ensure_directory(args.output_dir)
    all_records = []
    coral_path = _existing_path_or_none(args.coral, label="coral reference")
    if coral_path is not None:
        _, coral_records = prepare_coral_reference(coral_path, output_dir, limit=args.coral_limit)
        all_records.extend(coral_records)

    insect_path = _existing_path_or_none(args.insect_xlsx, label="insect workbook")
    if insect_path is not None:
        _, insect_records = prepare_insect_reference(insect_path, output_dir, limit=args.insect_limit)
        all_records.extend(insect_records)

    fungal_path = args.fungi_fasta if args.fungi_fasta is not None else args.fungal_fasta
    for source_name, source_path in (
        ("bacteria", args.bacteria_fasta),
        ("fungal", fungal_path),
        ("plant", args.plant_fasta),
    ):
        prepared_path = _existing_path_or_none(source_path, label=f"{source_name} FASTA")
        if prepared_path is None:
            continue
        _, extra_records = prepare_extra_reference(prepared_path, output_dir, source=source_name)
        all_records.extend(extra_records)

    for spec in args.extra_fasta or []:
        source, path = _parse_extra_reference(spec)
        _, extra_records = prepare_extra_reference(path, output_dir, source=source)
        all_records.extend(extra_records)

    if not all_records:
        raise ValueError(
            "No reference records were prepared. Provide at least one valid input via "
            "--coral / --insect-xlsx / --bacteria-fasta / --fungal-fasta / --plant-fasta / --extra-fasta."
        )

    metadata_path = write_reference_metadata(all_records, output_dir)
    logger.info("Prepared %d reference sequences in %s", len(all_records), output_dir)
    logger.info("Metadata: %s", metadata_path)
    return 0


def cmd_build_hmm(args: argparse.Namespace) -> int:
    """Build one query HMM from an aligned FASTA/MSA."""
    from ariadne.search import build_hmm

    hmm_path = build_hmm(args.alignment, args.output, name=args.name)
    logger.info("HMM written to %s", hmm_path)
    return 0


def cmd_build_tps_hmm_library(args: argparse.Namespace) -> int:
    """Build a directory of TPS HMM profiles from one or more alignments."""
    from ariadne.search import build_hmm
    from ariadne.utils import ensure_directory

    output_dir = ensure_directory(args.output_dir)
    built_paths: list[Path] = []
    for spec in args.alignment:
        if "=" in spec:
            name, alignment_path_text = spec.split("=", 1)
            alignment_path = Path(alignment_path_text).expanduser()
            hmm_name = name.strip()
        else:
            alignment_path = Path(spec).expanduser()
            hmm_name = alignment_path.stem
        output_path = output_dir / f"{hmm_name}.hmm"
        built_paths.append(build_hmm(alignment_path, output_path, name=hmm_name))
    logger.info("Built %d TPS HMM profiles in %s", len(built_paths), output_dir)
    for path in built_paths:
        logger.info("  %s", path)
    return 0


def cmd_discover(args: argparse.Namespace) -> int:
    """Run stage 1 candidate discovery from proteins or transcriptomes."""
    from ariadne.search import collect_protein_files, discover_candidates, discover_candidates_from_proteins

    if args.protein_folder:
        protein_paths = collect_protein_files(args.protein_folder, protein_glob=args.protein_glob)
        if not protein_paths:
            raise FileNotFoundError(
                f"No protein FASTA files were found in {args.protein_folder} "
                f"with patterns: {', '.join(args.protein_glob)}"
            )
        outputs = discover_candidates_from_proteins(
            protein_paths,
            args.hmm,
            args.output_dir,
            min_score=args.min_score,
            max_evalue=args.max_evalue,
        )
    elif args.transcriptomes:
        outputs = discover_candidates(
            args.transcriptomes,
            args.hmm,
            args.output_dir,
            min_score=args.min_score,
            max_evalue=args.max_evalue,
        )
    else:
        raise ValueError("Please provide either --protein-folder (preferred) or --transcriptomes for discover.")
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_filter(args: argparse.Namespace) -> int:
    """Run stage 2 candidate filtering and deduplication."""
    from ariadne.filter import filter_candidates

    outputs = filter_candidates(
        args.input_fasta,
        args.output_dir,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        identity_threshold=args.identity_threshold,
        reference_dir=args.reference_dir,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_classify(args: argparse.Namespace) -> int:
    """Run stage 3 feature-space classification."""
    from ariadne.embed import classify_candidates

    if args.tps_hmm_dir is not None:
        hmm_dir = Path(args.tps_hmm_dir)
    elif _bundled_tps_hmm_dir() is not None:
        hmm_dir = _bundled_tps_hmm_dir()
        logger.info("No --tps-hmm-dir provided; using bundled TPS HMM library from %s", hmm_dir)
    else:
        hmm_dir = _auto_build_tps_hmm_library(
            args.reference_dir,
            Path(args.output_dir) / "_auto_tps_hmms",
        )
    outputs = classify_candidates(
        args.candidates,
        args.reference_dir,
        args.output_dir,
        hmm_dir=hmm_dir,
        top_k=args.top_k,
        tree_neighbors=args.tree_neighbors,
        ceess_xlsx=None if args.skip_ceess_model else args.ceess_xlsx,
        ceess_model_name=args.ceess_model_name,
        ceess_batch_size=args.ceess_batch_size,
        ceess_max_length=args.ceess_max_length,
        ceess_device=args.ceess_device,
        ceess_cv_folds=args.ceess_cv_folds,
        ceess_random_state=args.ceess_random_state,
        ceess_threshold=args.ceess_threshold,
        ceess_classifier=args.ceess_classifier,
        ceess_epochs=args.ceess_epochs,
        ceess_hidden_dim=args.ceess_hidden_dim,
        ceess_representation_dim=args.ceess_barlow_representation_dim,
        ceess_projection_dim=args.ceess_barlow_projection_dim,
        ceess_dropout=args.ceess_dropout,
        ceess_learning_rate=args.ceess_learning_rate,
        ceess_weight_decay=args.ceess_weight_decay,
        ceess_train_batch_size=args.ceess_train_batch_size,
        ceess_barlow_redundancy_weight=args.ceess_barlow_redundancy_weight,
        ceess_mlp_checkpoint=args.ceess_mlp_checkpoint,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    """Execute the full Ariadne workflow end to end."""
    from ariadne.embed import classify_candidates
    from ariadne.search import collect_protein_files, discover_candidates, discover_candidates_from_proteins
    from ariadne.utils import ensure_directory, write_tsv
    from ariadne.filter import filter_candidates
    from ariadne.tree import build_phylogeny

    root = ensure_directory(args.output_dir)
    discovery_dir = ensure_directory(root / "01_discovery")
    filtering_dir = ensure_directory(root / "02_filtering")
    classification_dir = ensure_directory(root / "03_classification")
    phylogeny_dir = ensure_directory(root / "04_phylogeny") if not args.skip_phylogeny else None

    if args.query_hmm:
        hmm_path = Path(args.query_hmm)
    elif _bundled_query_hmm() is not None:
        hmm_path = _bundled_query_hmm()
        logger.info("No --query-hmm provided; using bundled discovery HMM from %s", hmm_path)
    else:
        hmm_path = _auto_build_query_hmm(args.reference_dir, discovery_dir / "query.hmm", name=args.hmm_name)

    if args.protein_folder:
        protein_paths = collect_protein_files(args.protein_folder, protein_glob=args.protein_glob)
        if not protein_paths:
            raise FileNotFoundError(
                f"No protein FASTA files were found in {args.protein_folder} "
                f"with patterns: {', '.join(args.protein_glob)}"
            )
        discovery_outputs = discover_candidates_from_proteins(
            protein_paths,
            hmm_path,
            discovery_dir,
            min_score=args.discovery_min_score,
            max_evalue=args.discovery_max_evalue,
        )
    elif args.transcriptomes:
        discovery_outputs = discover_candidates(
            args.transcriptomes,
            hmm_path,
            discovery_dir,
            min_score=args.discovery_min_score,
            max_evalue=args.discovery_max_evalue,
        )
    else:
        raise ValueError("Please provide either --protein-folder (preferred) or --transcriptomes for run.")
    filtering_outputs = filter_candidates(
        discovery_outputs["candidate_proteins"],
        filtering_dir,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        identity_threshold=args.identity_threshold,
        reference_dir=args.reference_dir,
    )
    if args.tps_hmm_dir is not None:
        tps_hmm_dir = Path(args.tps_hmm_dir)
    elif _bundled_tps_hmm_dir() is not None:
        tps_hmm_dir = _bundled_tps_hmm_dir()
        logger.info("No --tps-hmm-dir provided; using bundled TPS HMM library from %s", tps_hmm_dir)
    else:
        tps_hmm_dir = _auto_build_tps_hmm_library(
            args.reference_dir,
            classification_dir / "_auto_tps_hmms",
        )
    classification_outputs = classify_candidates(
        filtering_outputs["filtered_fasta"],
        args.reference_dir,
        classification_dir,
        hmm_dir=tps_hmm_dir,
        top_k=args.top_k,
        tree_neighbors=args.tree_neighbors,
        ceess_xlsx=None if args.skip_ceess_model else args.ceess_xlsx,
        ceess_model_name=args.ceess_model_name,
        ceess_batch_size=args.ceess_batch_size,
        ceess_max_length=args.ceess_max_length,
        ceess_device=args.ceess_device,
        ceess_cv_folds=args.ceess_cv_folds,
        ceess_random_state=args.ceess_random_state,
        ceess_threshold=args.ceess_threshold,
        ceess_classifier=args.ceess_classifier,
        ceess_epochs=args.ceess_epochs,
        ceess_hidden_dim=args.ceess_hidden_dim,
        ceess_representation_dim=args.ceess_barlow_representation_dim,
        ceess_projection_dim=args.ceess_barlow_projection_dim,
        ceess_dropout=args.ceess_dropout,
        ceess_learning_rate=args.ceess_learning_rate,
        ceess_weight_decay=args.ceess_weight_decay,
        ceess_train_batch_size=args.ceess_train_batch_size,
        ceess_barlow_redundancy_weight=args.ceess_barlow_redundancy_weight,
        ceess_mlp_checkpoint=args.ceess_mlp_checkpoint,
    )
    phylogeny_outputs = {}
    if not args.skip_phylogeny:
        phylogeny_outputs = build_phylogeny(
            filtering_outputs["filtered_fasta"],
            args.reference_dir,
            phylogeny_dir,
            mafft_bin=args.mafft_bin,
            mafft_mode=args.mafft_mode,
            iqtree_bin=args.iqtree_bin,
            iqtree_model=args.iqtree_model,
            iqtree_threads=args.iqtree_threads,
            iqtree_bootstrap=args.iqtree_bootstrap,
            iqtree_fast=not args.no_iqtree_fast,
        )
    summary_rows = []
    for label, outputs in (
        ("discovery", discovery_outputs),
        ("filtering", filtering_outputs),
        ("classification", classification_outputs),
        ("phylogeny", phylogeny_outputs),
    ):
        for key, value in outputs.items():
            summary_rows.append({"stage": label, "artifact": key, "path": value})
    summary_path = write_tsv(summary_rows, root / "pipeline_summary.tsv")
    logger.info("Pipeline completed. Summary: %s", summary_path)
    return 0


def cmd_phylogeny(args: argparse.Namespace) -> int:
    """Build a MAFFT alignment and IQ-TREE phylogeny from candidates plus references."""
    from ariadne.tree import build_phylogeny

    outputs = build_phylogeny(
        args.candidates,
        args.reference_dir,
        args.output_dir,
        mafft_bin=args.mafft_bin,
        mafft_mode=args.mafft_mode,
        iqtree_bin=args.iqtree_bin,
        iqtree_model=args.iqtree_model,
        iqtree_threads=args.iqtree_threads,
        iqtree_bootstrap=args.iqtree_bootstrap,
        iqtree_fast=not args.no_iqtree_fast,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the top-level parser and all Ariadne subcommands."""
    from ariadne.model import (
        DEFAULT_ESM_MODEL_NAME,
        esm_model_help_text,
    )

    parser = argparse.ArgumentParser(
        prog="ariadne",
        description="Ariadne — Terpene Synthase Discovery & Annotation Pipeline.",
    )
    parser.add_argument(
        "--version", "-V",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        default=False,
        help="Enable verbose (DEBUG level) logging.",
    )
    subparsers = parser.add_subparsers(dest="command", required=False)

    prepare_refs = subparsers.add_parser("prepare-references", help="Prepare clean reference FASTA files from the bundled coral and insect resources.")
    prepare_refs.add_argument("--coral", default=_default_coral(), type=Path)
    prepare_refs.add_argument("--coral-limit", type=int, default=None)
    prepare_refs.add_argument("--insect-xlsx", default=_default_insect(), type=Path)
    prepare_refs.add_argument("--insect-limit", type=int, default=None)
    prepare_refs.add_argument("--bacteria-fasta", type=Path, default=_default_bacteria(), help="Optional bacteria TPS reference FASTA.")
    prepare_refs.add_argument("--fungal-fasta", type=Path, default=_default_fungal(), help="Optional fungal TPS reference FASTA.")
    prepare_refs.add_argument("--fungi-fasta", type=Path, default=None, help="Alias of --fungal-fasta.")
    prepare_refs.add_argument("--plant-fasta", type=Path, default=_default_plant(), help="Optional plant TPS reference FASTA.")
    prepare_refs.add_argument("--extra-fasta", action="append", default=[], help="Additional reference FASTA in SOURCE=PATH form.")
    prepare_refs.add_argument("--output-dir", required=True, type=Path)
    prepare_refs.set_defaults(func=cmd_prepare_references)


    build_hmm_parser = subparsers.add_parser("build-hmm", help="Build a HMM from a reference FASTA/MSA source.")
    build_hmm_parser.add_argument("--alignment", required=True, type=Path)
    build_hmm_parser.add_argument("--output", required=True, type=Path)
    build_hmm_parser.add_argument("--name", default=None)
    build_hmm_parser.set_defaults(func=cmd_build_hmm)

    build_tps_lib = subparsers.add_parser(
        "build-tps-hmm-library",
        help="Build a TPS HMM library directory from reference FASTA/MSA files.",
    )
    build_tps_lib.add_argument(
        "--alignment",
        nargs="+",
        required=True,
        help="Reference FASTA/MSA path or NAME=PATH. Example: --alignment coral=tree/coral.fasta insect=tree/insect.fasta",
    )
    build_tps_lib.add_argument("--output-dir", required=True, type=Path)
    build_tps_lib.set_defaults(func=cmd_build_tps_hmm_library)

    discover = subparsers.add_parser("discover", help="Predict ORFs from transcriptomes and search them with a HMM.")
    discover.add_argument("--transcriptomes", nargs="+", default=None, type=Path)
    discover.add_argument("--protein-folder", type=Path, default=None, help="Folder containing Prodigal-predicted protein FASTA files.")
    discover.add_argument(
        "--protein-glob",
        nargs="+",
        default=["*.faa", "*.fa", "*.fasta", "*.pep", "*.prot"],
        help="Glob pattern(s) used to find protein FASTA files under --protein-folder.",
    )
    discover.add_argument("--hmm", required=True, type=Path)
    discover.add_argument("--output-dir", required=True, type=Path)
    discover.add_argument("--min-score", type=float, default=None)
    discover.add_argument("--max-evalue", type=float, default=None)
    discover.set_defaults(func=cmd_discover)

    filter_parser = subparsers.add_parser("filter", help="Coverage/length filtering and 95%% similarity deduplication.")
    filter_parser.add_argument("--input-fasta", required=True, type=Path)
    filter_parser.add_argument("--output-dir", required=True, type=Path)
    filter_parser.add_argument("--min-coverage", type=float, default=10.0)
    filter_parser.add_argument("--min-length", type=int, default=300)
    filter_parser.add_argument("--identity-threshold", type=float, default=0.95)
    filter_parser.add_argument("--reference-dir", type=Path, default=None, help="Optional reference FASTA directory used to remove candidates that already match known references.")
    filter_parser.set_defaults(func=cmd_filter)

    classify = subparsers.add_parser("classify", help="Classify TPS candidates in profile feature space.")
    classify.add_argument("--candidates", required=True, type=Path)
    classify.add_argument("--reference-dir", required=True, type=Path)
    classify.add_argument("--output-dir", required=True, type=Path)
    classify.add_argument("--tps-hmm-dir", default=None, type=Path, help="Optional directory containing TPS HMM profiles (*.hmm). When omitted, Ariadne builds them from --reference-dir.")
    classify.add_argument("--top-k", type=int, default=5)
    classify.add_argument("--tree-neighbors", type=int, default=12)
    classify.add_argument("--ceess-xlsx", type=Path, default=_default_tps_xlsx(), help="Optional coral TPS workbook used to train the ESM CeeSs model. Defaults to TPS/TPS.xlsx when present.")
    classify.add_argument("--skip-ceess-model", action="store_true", help="Skip the optional ESM-based CeeSs scoring stage.")
    classify.add_argument("--ceess-threshold", type=float, default=0.9, help="Probability threshold used to keep predicted CeeSs candidates.")
    classify.add_argument("--ceess-classifier", choices=["mlp", "logreg", "contrastive"], default="mlp", help="Classifier pipeline used on top of frozen ESM embeddings.")
    classify.add_argument("--ceess-model-name", default=DEFAULT_ESM_MODEL_NAME, help=esm_model_help_text())
    classify.add_argument("--ceess-batch-size", type=int, default=4)
    classify.add_argument("--ceess-max-length", type=int, default=2048)
    classify.add_argument("--ceess-device", default=None, help="Optional torch device for the ESM CeeSs model, for example cpu or cuda.")
    classify.add_argument("--ceess-cv-folds", type=int, default=5)
    classify.add_argument("--ceess-random-state", type=int, default=0)
    classify.add_argument("--ceess-epochs", type=int, default=200, help="Training epochs for the MLP CeeSs head.")
    classify.add_argument("--ceess-hidden-dim", type=int, default=128, help="Hidden layer width for the MLP CeeSs head.")
    classify.add_argument("--ceess-barlow-representation-dim", type=int, default=None, help="Optional encoded representation width for the Barlow Twins CeeSs pipeline.")
    classify.add_argument("--ceess-barlow-projection-dim", type=int, default=None, help="Optional projection-head output width for the Barlow Twins CeeSs pipeline.")
    classify.add_argument("--ceess-barlow-redundancy-weight", type=float, default=0.005, help="Off-diagonal redundancy penalty used when --ceess-classifier=contrastive.")
    classify.add_argument("--ceess-dropout", type=float, default=0.1, help="Dropout rate for the MLP CeeSs head.")
    classify.add_argument("--ceess-learning-rate", type=float, default=1e-3, help="Learning rate for the MLP CeeSs head.")
    classify.add_argument("--ceess-weight-decay", type=float, default=1e-4, help="Weight decay for the MLP CeeSs head.")
    classify.add_argument("--ceess-train-batch-size", type=int, default=8, help="Training batch size for the MLP CeeSs head.")
    classify.add_argument("--ceess-mlp-checkpoint", type=Path, default=None, help="Optional pretrained Torch MLP checkpoint (.pt) for --ceess-classifier mlp. When provided, Ariadne skips final MLP training and loads this classifier directly.")
    classify.set_defaults(func=cmd_classify)

    run = subparsers.add_parser("run", help="Execute Ariadne end-to-end: discovery, filtering, classification, and optional phylogeny.")
    run.add_argument("--transcriptomes", nargs="+", default=None, type=Path)
    run.add_argument("--protein-folder", type=Path, default=None, help="Preferred input mode: folder of predicted protein FASTA files.")
    run.add_argument(
        "--protein-glob",
        nargs="+",
        default=["*.faa", "*.fa", "*.fasta", "*.pep", "*.prot"],
        help="Glob pattern(s) used to find protein FASTA files under --protein-folder.",
    )
    run.add_argument("--query-hmm", type=Path, default=None, help="Optional prebuilt query HMM. When omitted, Ariadne builds one from the reference-dir FASTA set.")
    run.add_argument("--reference-dir", required=True, type=Path)
    run.add_argument("--output-dir", required=True, type=Path)
    run.add_argument("--hmm-name", default="ariadne_query")
    run.add_argument("--discovery-min-score", type=float, default=None)
    run.add_argument("--discovery-max-evalue", type=float, default=None)
    run.add_argument("--min-coverage", type=float, default=10.0)
    run.add_argument("--min-length", type=int, default=300)
    run.add_argument("--identity-threshold", type=float, default=0.95)
    run.add_argument("--tps-hmm-dir", default=None, type=Path, help="Optional directory containing TPS HMM profiles (*.hmm). When omitted, Ariadne builds them from --reference-dir.")
    run.add_argument("--top-k", type=int, default=5)
    run.add_argument("--tree-neighbors", type=int, default=12)
    run.add_argument("--ceess-xlsx", type=Path, default=_default_tps_xlsx(), help="Optional coral TPS workbook used to train the ESM CeeSs model. Defaults to TPS/TPS.xlsx when present.")
    run.add_argument("--skip-ceess-model", action="store_true", help="Skip the optional ESM-based CeeSs scoring stage.")
    run.add_argument("--ceess-threshold", type=float, default=0.9, help="Probability threshold used to keep predicted CeeSs candidates.")
    run.add_argument("--ceess-classifier", choices=["mlp", "logreg", "contrastive"], default="mlp", help="Classifier pipeline used on top of frozen ESM embeddings.")
    run.add_argument("--ceess-model-name", default=DEFAULT_ESM_MODEL_NAME, help=esm_model_help_text())
    run.add_argument("--ceess-batch-size", type=int, default=4)
    run.add_argument("--ceess-max-length", type=int, default=2048)
    run.add_argument("--ceess-device", default=None, help="Optional torch device for the ESM CeeSs model, for example cpu or cuda.")
    run.add_argument("--ceess-cv-folds", type=int, default=5)
    run.add_argument("--ceess-random-state", type=int, default=0)
    run.add_argument("--ceess-epochs", type=int, default=200, help="Training epochs for the MLP CeeSs head.")
    run.add_argument("--ceess-hidden-dim", type=int, default=128, help="Hidden layer width for the MLP CeeSs head.")
    run.add_argument("--ceess-barlow-representation-dim", type=int, default=None, help="Optional encoded representation width for the Barlow Twins CeeSs pipeline.")
    run.add_argument("--ceess-barlow-projection-dim", type=int, default=None, help="Optional projection-head output width for the Barlow Twins CeeSs pipeline.")
    run.add_argument("--ceess-barlow-redundancy-weight", type=float, default=0.005, help="Off-diagonal redundancy penalty used when --ceess-classifier=contrastive.")
    run.add_argument("--ceess-dropout", type=float, default=0.1, help="Dropout rate for the MLP CeeSs head.")
    run.add_argument("--ceess-learning-rate", type=float, default=1e-3, help="Learning rate for the MLP CeeSs head.")
    run.add_argument("--ceess-weight-decay", type=float, default=1e-4, help="Weight decay for the MLP CeeSs head.")
    run.add_argument("--ceess-train-batch-size", type=int, default=8, help="Training batch size for the MLP CeeSs head.")
    run.add_argument("--ceess-mlp-checkpoint", type=Path, default=None, help="Optional pretrained Torch MLP checkpoint (.pt) for --ceess-classifier mlp. When provided, Ariadne skips final MLP training and loads this classifier directly.")
    run.add_argument("--skip-phylogeny", action="store_true", help="Skip the MAFFT + IQ-TREE phylogeny step.")
    run.add_argument("--mafft-bin", default=None, help="Path or executable name for MAFFT.")
    run.add_argument("--mafft-mode", default="--auto", help="MAFFT mode flag, for example --auto or --localpair.")
    run.add_argument("--iqtree-bin", default=None, help="Path or executable name for IQ-TREE.")
    run.add_argument("--iqtree-model", default="LG", help="IQ-TREE substitution model setting.")
    run.add_argument("--iqtree-threads", default="AUTO", help="IQ-TREE thread setting, for example AUTO or 8.")
    run.add_argument("--iqtree-bootstrap", type=int, default=None, help="Optional IQ-TREE ultrafast bootstrap replicates.")
    run.add_argument(
        "--no-iqtree-fast",
        action="store_true",
        help="Disable IQ-TREE fast mode. By default Ariadne uses --fast for practical end-to-end runs.",
    )
    run.set_defaults(func=cmd_run)

    phylogeny = subparsers.add_parser("phylogeny", help="Build a MAFFT alignment and IQ-TREE phylogeny from candidates plus references.")
    phylogeny.add_argument("--candidates", required=True, type=Path)
    phylogeny.add_argument("--reference-dir", required=True, type=Path)
    phylogeny.add_argument("--output-dir", required=True, type=Path)
    phylogeny.add_argument("--mafft-bin", default=None, help="Path or executable name for MAFFT.")
    phylogeny.add_argument("--mafft-mode", default="--auto", help="MAFFT mode flag, for example --auto or --localpair.")
    phylogeny.add_argument("--iqtree-bin", default=None, help="Path or executable name for IQ-TREE.")
    phylogeny.add_argument("--iqtree-model", default="LG", help="IQ-TREE substitution model setting.")
    phylogeny.add_argument("--iqtree-threads", default="AUTO", help="IQ-TREE thread setting, for example AUTO or 8.")
    phylogeny.add_argument("--iqtree-bootstrap", type=int, default=None, help="Optional IQ-TREE ultrafast bootstrap replicates.")
    phylogeny.add_argument(
        "--no-iqtree-fast",
        action="store_true",
        help="Disable IQ-TREE fast mode. By default Ariadne uses --fast for practical end-to-end runs.",
    )
    phylogeny.set_defaults(func=cmd_phylogeny)

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    """Entry point shared by ``python -m ariadne`` and the console script."""
    parser = build_parser()
    args = parser.parse_args(argv)

    _log.setup_logging(verbose=getattr(args, "verbose", False))

    if not hasattr(args, "func"):
        # No subcommand supplied: print banner + help then exit cleanly.
        _log.print_banner(__version__)
        parser.print_help(sys.stderr)
        return 0

    _log.print_banner(__version__)
    return args.func(args)
