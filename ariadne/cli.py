"""Command-line interface for the Ariadne TPS discovery pipeline.

The CLI is organised around pipeline stages so users can either run the full
workflow with ``ariadne run`` or execute individual modules for debugging and
method development.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Tuple, Union

from ariadne import __version__
from ariadne import log as _log

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


def _default_tps_hmms() -> Path:
    """Return the bundled TPS HMM directory."""
    bundled = Path(__file__).resolve().parent / "tps_hmm"
    return bundled


def _default_reference_dir() -> Optional[Path]:
    """Locate the default tree/reference directory when present."""
    return _first_existing_path(
        Path.cwd() / "tree",
        _repo_root() / "tree",
    )


def _is_bundled_tps_hmms(path: PathLike) -> bool:
    """Check whether a path points at Ariadne's bundled placeholder HMMs."""
    return Path(path).resolve() == _default_tps_hmms().resolve()


def _warn_legacy_tps_hmms(path: PathLike) -> None:
    """Warn users that bundled HMMs are only placeholders."""
    if not _is_bundled_tps_hmms(path):
        return
    logger.warning(
        "Using bundled TPS HMMs under ariadne/tps_hmm. "
        "These profiles are legacy AFLP-derived placeholders. "
        "For publication-grade TPS analysis, build and pass your own --tps-hmm-dir."
    )


def _default_reference_alignment(reference_dir: Optional[PathLike] = None) -> Optional[Path]:
    """Locate the default motif/reference FASTA, preferring coral data under ``tree/``."""
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
    """Infer a usable motif/reference FASTA from the reference directory."""
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
        f"Could not locate a reference FASTA for motif analysis under {directory}. "
        "Expected something like coral.fasta."
    )


def _reference_fasta_paths(reference_dir: PathLike) -> list[Path]:
    """Collect FASTA files from a prepared reference directory."""
    directory = Path(reference_dir)
    fasta_paths = sorted(directory.glob("*.fa*"))
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files were found in reference directory: {directory}")
    return fasta_paths


def _auto_build_query_hmm(reference_dir: PathLike, output_path: PathLike, *, name: str) -> Path:
    """Build the discovery query HMM from the default reference FASTA in ``reference_dir``."""
    from ariadne.discovery import build_hmm

    source_fasta = _find_reference_alignment(reference_dir)
    logger.info("No --query-hmm provided; building query HMM from %s", source_fasta)
    return build_hmm(source_fasta, output_path, name=name)


def _auto_build_tps_hmm_library(reference_dir: PathLike, output_dir: PathLike) -> Path:
    """Build a TPS HMM library from all FASTA files in ``reference_dir``."""
    from ariadne.discovery import build_hmm
    from ariadne.fasta_utils import ensure_directory

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
    from ariadne.fasta_utils import ensure_directory
    from ariadne.references import (
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


def cmd_prepare_demo(args: argparse.Namespace) -> int:
    """Create a compact demo workspace with synthetic inputs and prepared refs."""
    from ariadne.demo import prepare_demo_workspace

    if args.coral is None or not Path(args.coral).exists():
        raise FileNotFoundError(
            "prepare-demo requires a coral alignment FASTA. "
            "Please pass --coral explicitly (e.g., --coral './coralTPS (modified)-cembrene.fasta')."
        )
    if args.insect_xlsx is None or not Path(args.insect_xlsx).exists():
        raise FileNotFoundError(
            "prepare-demo requires an insect workbook. "
            "Please pass --insect-xlsx explicitly (e.g., --insect-xlsx './Insecta TPS.xlsx')."
        )

    outputs = prepare_demo_workspace(
        args.output_dir,
        coral_alignment_fasta=args.coral,
        insect_xlsx=args.insect_xlsx,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_build_hmm(args: argparse.Namespace) -> int:
    """Build one query HMM from an aligned FASTA/MSA."""
    from ariadne.discovery import build_hmm

    hmm_path = build_hmm(args.alignment, args.output, name=args.name)
    logger.info("HMM written to %s", hmm_path)
    return 0


def cmd_build_tps_hmm_library(args: argparse.Namespace) -> int:
    """Build a directory of TPS HMM profiles from one or more alignments."""
    from ariadne.discovery import build_hmm
    from ariadne.fasta_utils import ensure_directory

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
    from ariadne.discovery import collect_protein_files, discover_candidates, discover_candidates_from_proteins

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
    from ariadne.filtering import filter_candidates

    outputs = filter_candidates(
        args.input_fasta,
        args.output_dir,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        identity_threshold=args.identity_threshold,
        motif_anchor=args.motif_anchor,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_classify(args: argparse.Namespace) -> int:
    """Run stage 3 feature-space classification."""
    from ariadne.classification import classify_candidates

    hmm_dir = Path(args.tps_hmm_dir) if args.tps_hmm_dir is not None else _auto_build_tps_hmm_library(
        args.reference_dir,
        Path(args.output_dir) / "_auto_tps_hmms",
    )
    _warn_legacy_tps_hmms(hmm_dir)
    outputs = classify_candidates(
        args.candidates,
        args.reference_dir,
        args.output_dir,
        hmm_dir=hmm_dir,
        candidate_annotations=args.candidate_annotations,
        top_k=args.top_k,
        tree_neighbors=args.tree_neighbors,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_motif(args: argparse.Namespace) -> int:
    """Run stage 4 motif-based TPS and cembrene analysis."""
    from ariadne.motif import analyze_motifs

    if args.reference_alignment is None:
        raise ValueError("Please provide --reference-alignment for motif analysis.")
    outputs = analyze_motifs(
        args.candidates,
        args.reference_alignment,
        args.output_dir,
        validated_cess_fasta=args.validated_cess_fasta,
        tps_anchor_pattern=args.tps_pattern,
        tps_center_position=args.tps_center_position,
        tps_search_radius=args.tps_search_radius,
        anchor_pattern=args.anchor_pattern,
        flank=args.flank,
        center_position=args.center_position,
        allow_center_fallback=args.allow_center_fallback,
        validated_cess_identity_threshold=args.validated_cess_identity_threshold,
        cembrene_identity_threshold=args.cembrene_identity_threshold,
        cembrene_margin_threshold=args.cembrene_margin_threshold,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    """Execute the full Ariadne workflow end to end."""
    from ariadne.benchmark import compare_fasta_sets
    from ariadne.classification import classify_candidates
    from ariadne.discovery import collect_protein_files, discover_candidates, discover_candidates_from_proteins
    from ariadne.fasta_utils import ensure_directory, write_tsv
    from ariadne.filtering import filter_candidates
    from ariadne.motif import analyze_motifs
    from ariadne.phylogeny import build_phylogeny

    root = ensure_directory(args.output_dir)
    discovery_dir = ensure_directory(root / "01_discovery")
    filtering_dir = ensure_directory(root / "02_filtering")
    classification_dir = ensure_directory(root / "03_classification")
    motif_dir = ensure_directory(root / "04_motif")
    if args.enable_benchmark and args.expected_fasta is None:
        raise ValueError("Benchmark mode requires --expected-fasta. Pass --enable-benchmark together with a benchmark FASTA.")

    if args.query_hmm:
        hmm_path = Path(args.query_hmm)
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
        motif_anchor=args.motif_anchor,
    )
    reference_alignment = Path(args.reference_alignment) if args.reference_alignment else _find_reference_alignment(args.reference_dir)
    motif_outputs = analyze_motifs(
        filtering_outputs["filtered_fasta"],
        reference_alignment,
        motif_dir,
        validated_cess_fasta=args.validated_cess_fasta if args.validated_cess_fasta is not None else args.expected_fasta,
        tps_anchor_pattern=args.tps_pattern,
        tps_center_position=args.tps_center_position,
        tps_search_radius=args.tps_search_radius,
        anchor_pattern=args.anchor_pattern,
        flank=args.flank,
        center_position=args.center_position,
        allow_center_fallback=args.allow_center_fallback,
        validated_cess_identity_threshold=args.validated_cess_identity_threshold,
        cembrene_identity_threshold=args.cembrene_identity_threshold,
        cembrene_margin_threshold=args.cembrene_margin_threshold,
    )
    tps_hmm_dir = Path(args.tps_hmm_dir) if args.tps_hmm_dir is not None else _auto_build_tps_hmm_library(
        args.reference_dir,
        classification_dir / "_auto_tps_hmms",
    )
    _warn_legacy_tps_hmms(tps_hmm_dir)
    classification_outputs = classify_candidates(
        filtering_outputs["filtered_fasta"],
        args.reference_dir,
        classification_dir,
        hmm_dir=tps_hmm_dir,
        candidate_annotations=motif_outputs["motif_summary"],
        top_k=args.top_k,
        tree_neighbors=args.tree_neighbors,
    )
    phylogeny_outputs = {}
    if not args.skip_phylogeny:
        phylogeny_dir = ensure_directory(root / "06_phylogeny")
        phylogeny_candidates = (
            filtering_outputs["filtered_fasta"]
            if args.phylogeny_candidates == "filtered"
            else motif_outputs["cembrene_fasta"]
        )
        phylogeny_outputs = build_phylogeny(
            phylogeny_candidates,
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
    benchmark_outputs = {}
    if args.enable_benchmark and args.expected_fasta:
        benchmark_dir = ensure_directory(root / "05_benchmark")
        benchmark_outputs = compare_fasta_sets(
            motif_outputs["cembrene_fasta"],
            args.expected_fasta,
            benchmark_dir,
        )
    summary_rows = []
    for label, outputs in (
        ("discovery", discovery_outputs),
        ("filtering", filtering_outputs),
        ("motif", motif_outputs),
        ("classification", classification_outputs),
        ("benchmark", benchmark_outputs),
        ("phylogeny", phylogeny_outputs),
    ):
        for key, value in outputs.items():
            summary_rows.append({"stage": label, "artifact": key, "path": value})
    summary_path = write_tsv(summary_rows, root / "pipeline_summary.tsv")
    logger.info("Pipeline completed. Summary: %s", summary_path)
    return 0


def cmd_compare_fasta(args: argparse.Namespace) -> int:
    """Benchmark one predicted FASTA against an expected FASTA set."""
    from ariadne.benchmark import compare_fasta_sets

    outputs = compare_fasta_sets(args.predicted_fasta, args.expected_fasta, args.output_dir)
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_visualize(args: argparse.Namespace) -> int:
    """Generate legacy t-SNE style visualisations from feature tables."""
    from ariadne.visualization import visualize_profiles

    outputs = visualize_profiles(
        args.input_table,
        args.output_dir,
        perplexities=args.perplexities,
        min_points=args.min_points,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_phylogeny(args: argparse.Namespace) -> int:
    """Build a MAFFT alignment and IQ-TREE phylogeny from candidates plus references."""
    from ariadne.phylogeny import build_phylogeny

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


def cmd_filter_coverage_only(args: argparse.Namespace) -> int:
    """Compatibility wrapper for coverage-only filtering."""
    from ariadne.fasta_utils import read_fasta, write_fasta
    from ariadne.filtering import filter_by_coverage

    records = filter_by_coverage(read_fasta(args.input_fasta), args.min_coverage)
    output_path = write_fasta(records, args.output)
    logger.info("Output: %s", output_path)
    return 0


def cmd_filter_length_only(args: argparse.Namespace) -> int:
    """Compatibility wrapper for length-only filtering."""
    from ariadne.fasta_utils import read_fasta, write_fasta
    from ariadne.filtering import filter_by_length

    records = filter_by_length(read_fasta(args.input_fasta), args.min_length)
    output_path = write_fasta(records, args.output)
    logger.info("Output: %s", output_path)
    return 0


def cmd_dedupe_exact(args: argparse.Namespace) -> int:
    """Compatibility wrapper for exact-sequence deduplication."""
    from ariadne.fasta_utils import read_fasta, write_fasta
    from ariadne.filtering import deduplicate_exact

    records = deduplicate_exact(read_fasta(args.input_fasta))
    output_path = write_fasta(records, args.output)
    logger.info("Output: %s", output_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    """Create the top-level parser and all Ariadne subcommands."""
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

    prepare_demo = subparsers.add_parser("prepare-demo", help="Create a small fully-runnable demo workspace.")
    prepare_demo.add_argument("--output-dir", required=True, type=Path)
    prepare_demo.add_argument("--coral", default=_default_coral(), type=Path)
    prepare_demo.add_argument("--insect-xlsx", default=_default_insect(), type=Path)
    prepare_demo.set_defaults(func=cmd_prepare_demo)

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
    filter_parser.add_argument("--motif-anchor", default="CFDVL")
    filter_parser.set_defaults(func=cmd_filter)

    classify = subparsers.add_parser("classify", help="Classify TPS candidates in profile feature space.")
    classify.add_argument("--candidates", required=True, type=Path)
    classify.add_argument("--reference-dir", required=True, type=Path)
    classify.add_argument("--output-dir", required=True, type=Path)
    classify.add_argument("--tps-hmm-dir", default=None, type=Path, help="Optional directory containing TPS HMM profiles (*.hmm). When omitted, Ariadne builds them from --reference-dir.")
    classify.add_argument("--candidate-annotations", type=Path, default=None, help="Optional motif/benchmark annotation TSV keyed by sequence_id.")
    classify.add_argument("--top-k", type=int, default=5)
    classify.add_argument("--tree-neighbors", type=int, default=12)
    classify.set_defaults(func=cmd_classify)

    motif = subparsers.add_parser(
        "motif",
        help="Two-stage motif gate: first confirm TPS by DDXXD/E near 125 aa, then CeSS-oriented cembrene judgement (center fallback enabled by default).",
    )
    motif.add_argument("--candidates", required=True, type=Path)
    motif.add_argument(
        "--reference-alignment",
        default=_default_reference_alignment(),
        type=Path,
        help="Reference FASTA for motif analysis (defaults to tree/coral.fasta when available).",
    )
    motif.add_argument("--validated-cess-fasta", type=Path, default=None, help="Optional validated CeSS FASTA used to calibrate fallback calling and candidate ranking.")
    motif.add_argument("--output-dir", required=True, type=Path)
    motif.add_argument("--tps-pattern", default=r"DD..[DE]")
    motif.add_argument("--tps-center-position", type=int, default=125)
    motif.add_argument("--tps-search-radius", type=int, default=40)
    motif.add_argument("--anchor-pattern", default=r"CFDVL.")
    motif.add_argument("--flank", type=int, default=10)
    motif.add_argument("--center-position", type=int, default=210)
    motif.add_argument("--validated-cess-identity-threshold", type=float, default=0.95, help="Minimum motif-window identity to a validated CeSS reference required for high-confidence CeSS calls.")
    motif.add_argument("--cembrene-identity-threshold", type=float, default=0.75, help="Minimum motif-window identity to a cembrene-family reference for family-level support.")
    motif.add_argument("--cembrene-margin-threshold", type=float, default=0.1, help="Minimum identity margin over non-cembrene references for family-level support.")
    motif.set_defaults(allow_center_fallback=True)
    motif.add_argument(
        "--allow-center-fallback",
        dest="allow_center_fallback",
        action="store_true",
        help="Enable center-position fallback for cembrene comparison. This is the default.",
    )
    motif.add_argument(
        "--disable-center-fallback",
        dest="allow_center_fallback",
        action="store_false",
        help="Disable center-position fallback and require a direct anchor match.",
    )
    motif.set_defaults(func=cmd_motif)

    run = subparsers.add_parser("run", help="Execute Ariadne end-to-end, including optional phylogeny.")
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
    run.add_argument("--expected-fasta", type=Path, default=None, help="Optional benchmark FASTA. Ariadne only runs benchmark when --enable-benchmark is set.")
    run.add_argument(
        "--reference-alignment",
        type=Path,
        default=None,
        help="Reference FASTA for module 4. Defaults to coral.fasta under --reference-dir when available.",
    )
    run.add_argument("--validated-cess-fasta", type=Path, default=None, help="Optional validated CeSS FASTA for calibrated fallback calling. When omitted, --expected-fasta is reused if present.")
    run.add_argument("--hmm-name", default="ariadne_query")
    run.add_argument("--discovery-min-score", type=float, default=None)
    run.add_argument("--discovery-max-evalue", type=float, default=None)
    run.add_argument("--min-coverage", type=float, default=10.0)
    run.add_argument("--min-length", type=int, default=300)
    run.add_argument("--identity-threshold", type=float, default=0.95)
    run.add_argument("--motif-anchor", default="CFDVL")
    run.add_argument("--tps-pattern", default=r"DD..[DE]")
    run.add_argument("--tps-center-position", type=int, default=125)
    run.add_argument("--tps-search-radius", type=int, default=40)
    run.add_argument("--anchor-pattern", default=r"CFDVL.")
    run.add_argument("--flank", type=int, default=10)
    run.add_argument("--center-position", type=int, default=210)
    run.add_argument("--validated-cess-identity-threshold", type=float, default=0.95)
    run.add_argument("--cembrene-identity-threshold", type=float, default=0.75)
    run.add_argument("--cembrene-margin-threshold", type=float, default=0.1)
    run.add_argument("--enable-benchmark", action="store_true", help="Opt in to benchmark output generation under 05_benchmark/. Disabled by default for routine mining runs.")
    run.set_defaults(allow_center_fallback=True)
    run.add_argument(
        "--allow-center-fallback",
        dest="allow_center_fallback",
        action="store_true",
        help="Enable center-position fallback when anchor motif is absent. This is the default.",
    )
    run.add_argument(
        "--disable-center-fallback",
        dest="allow_center_fallback",
        action="store_false",
        help="Disable center-position fallback and require a direct anchor match.",
    )
    run.add_argument("--tps-hmm-dir", default=None, type=Path, help="Optional directory containing TPS HMM profiles (*.hmm). When omitted, Ariadne builds them from --reference-dir.")
    run.add_argument("--top-k", type=int, default=5)
    run.add_argument("--tree-neighbors", type=int, default=12)
    run.add_argument("--skip-phylogeny", action="store_true", help="Skip the MAFFT + IQ-TREE phylogeny step.")
    run.add_argument("--phylogeny-candidates", choices=["filtered", "cembrene"], default="filtered", help="Which candidate FASTA to use for MAFFT/IQ-TREE phylogeny.")
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

    compare_fasta = subparsers.add_parser("compare-fasta", help="Compare a predicted FASTA against an expected FASTA benchmark.")
    compare_fasta.add_argument("--predicted-fasta", required=True, type=Path)
    compare_fasta.add_argument("--expected-fasta", required=True, type=Path)
    compare_fasta.add_argument("--output-dir", required=True, type=Path)
    compare_fasta.set_defaults(func=cmd_compare_fasta)

    visualize = subparsers.add_parser(
        "visualize",
        help="Generate t-SNE style visualization from TPS HMM score tables.",
    )
    visualize.add_argument(
        "--input-table",
        required=True,
        type=Path,
        help="Input TPS score table: supports run/classify output tps_features.tsv and legacy hmm_tab.txt style.",
    )
    visualize.add_argument("--output-dir", required=True, type=Path)
    visualize.add_argument("--perplexities", nargs="+", type=int, default=None, help="Explicit t-SNE perplexities. Default follows legacy AFLP heuristic.")
    visualize.add_argument("--min-points", type=int, default=20, help="Clustering min points, equivalent to legacy minPts.")
    visualize.set_defaults(func=cmd_visualize)

    coverage_only = subparsers.add_parser("filter-coverage-only", help="Compatibility helper for the legacy coverage filter.")
    coverage_only.add_argument("min_coverage", type=float)
    coverage_only.add_argument("input_fasta", type=Path)
    coverage_only.add_argument("--output", required=True, type=Path)
    coverage_only.set_defaults(func=cmd_filter_coverage_only)

    length_only = subparsers.add_parser("filter-length-only", help="Compatibility helper for the legacy length filter.")
    length_only.add_argument("min_length", type=int)
    length_only.add_argument("input_fasta", type=Path)
    length_only.add_argument("--output", required=True, type=Path)
    length_only.set_defaults(func=cmd_filter_length_only)

    dedupe_only = subparsers.add_parser("dedupe-exact", help="Compatibility helper for exact duplicate removal.")
    dedupe_only.add_argument("input_fasta", type=Path)
    dedupe_only.add_argument("--output", required=True, type=Path)
    dedupe_only.set_defaults(func=cmd_dedupe_exact)
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
