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
    return Path(__file__).resolve().parent.parent


def _first_existing_path(*candidates: PathLike) -> Optional[Path]:
    for candidate in candidates:
        if candidate is None:
            continue
        path = Path(candidate).expanduser().resolve()
        if path.exists():
            return path
    return None


def _default_coral() -> Optional[Path]:
    return _first_existing_path(
        Path.cwd() / "coralTPS (modified)-cembrene.fasta",
        _repo_root() / "coralTPS (modified)-cembrene.fasta",
        Path(__file__).resolve().parent / "coralTPS (modified)-cembrene.fasta",
    )


def _default_insect() -> Optional[Path]:
    return _first_existing_path(
        Path.cwd() / "Insecta TPS.xlsx",
        _repo_root() / "Insecta TPS.xlsx",
        Path(__file__).resolve().parent / "Insecta TPS.xlsx",
    )


def _default_bacteria() -> Optional[Path]:
    return _first_existing_path(
        Path.cwd() / "bacteria.fasta",
        Path.cwd() / "bacteria.fa",
        Path.cwd() / "bacteria.faa",
        _repo_root() / "bacteria.fasta",
        _repo_root() / "bacteria.fa",
        _repo_root() / "bacteria.faa",
    )


def _default_fungal() -> Optional[Path]:
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
    return _first_existing_path(
        Path.cwd() / "plant.fasta",
        Path.cwd() / "plant.fa",
        Path.cwd() / "plant.faa",
        _repo_root() / "plant.fasta",
        _repo_root() / "plant.fa",
        _repo_root() / "plant.faa",
    )


def _default_tps_hmms() -> Path:
    bundled = Path(__file__).resolve().parent / "tps_hmm"
    return bundled


def _is_bundled_tps_hmms(path: PathLike) -> bool:
    return Path(path).resolve() == _default_tps_hmms().resolve()


def _warn_legacy_tps_hmms(path: PathLike) -> None:
    if not _is_bundled_tps_hmms(path):
        return
    logger.warning(
        "Using bundled TPS HMMs under ariadne/tps_hmm. "
        "These profiles are legacy AFLP-derived placeholders. "
        "For publication-grade TPS analysis, build and pass your own --tps-hmm-dir."
    )


def _default_alignment_reference() -> Optional[Path]:
    candidates = [
        Path.cwd() / "Alignment.fasta",
        _repo_root() / "Alignment.fasta",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _find_reference_alignment(reference_dir: PathLike) -> Path:
    default_alignment = _default_alignment_reference()
    if default_alignment is not None and default_alignment.exists():
        return default_alignment
    directory = Path(reference_dir)
    for filename in ("Alignment.fasta", "alignment.fasta", "coral.fasta"):
        direct = directory / filename
        if direct.exists():
            return direct
    for fasta_path in sorted(directory.glob("*.fa*")):
        lowered = fasta_path.name.lower()
        if "alignment" in lowered or "coral" in lowered:
            return fasta_path
    raise FileNotFoundError(
        f"Could not locate reference alignment FASTA. "
        f"Checked current folder Alignment.fasta and {directory}."
    )


def _parse_extra_reference(spec: str) -> Tuple[str, str]:
    if "=" not in spec:
        raise ValueError(f"Expected SOURCE=PATH for extra references, got: {spec}")
    source, path = spec.split("=", 1)
    return source.strip(), path.strip()


def _existing_path_or_none(value: Optional[PathLike], *, label: str) -> Optional[Path]:
    if value is None:
        return None
    path = Path(value).expanduser()
    if path.exists():
        return path
    logger.warning("%s file not found, skipped: %s", label, path)
    return None


def cmd_prepare_references(args: argparse.Namespace) -> int:
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
    from ariadne.discovery import build_hmm

    hmm_path = build_hmm(args.alignment, args.output, name=args.name)
    logger.info("HMM written to %s", hmm_path)
    return 0


def cmd_build_tps_hmm_library(args: argparse.Namespace) -> int:
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
    from ariadne.classification import classify_candidates

    _warn_legacy_tps_hmms(args.tps_hmm_dir)
    outputs = classify_candidates(
        args.candidates,
        args.reference_dir,
        args.output_dir,
        hmm_dir=args.tps_hmm_dir,
        candidate_annotations=args.candidate_annotations,
        top_k=args.top_k,
        tree_neighbors=args.tree_neighbors,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_motif(args: argparse.Namespace) -> int:
    from ariadne.motif import analyze_motifs

    if args.reference_alignment is None:
        raise ValueError("Please provide --reference-alignment for motif analysis.")
    outputs = analyze_motifs(
        args.candidates,
        args.reference_alignment,
        args.output_dir,
        tps_anchor_pattern=args.tps_pattern,
        tps_center_position=args.tps_center_position,
        tps_search_radius=args.tps_search_radius,
        anchor_pattern=args.anchor_pattern,
        flank=args.flank,
        center_position=args.center_position,
        allow_center_fallback=args.allow_center_fallback,
    )
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    from ariadne.benchmark import compare_fasta_sets
    from ariadne.classification import classify_candidates
    from ariadne.discovery import build_hmm, collect_protein_files, discover_candidates, discover_candidates_from_proteins
    from ariadne.fasta_utils import ensure_directory, write_tsv
    from ariadne.filtering import filter_candidates
    from ariadne.motif import analyze_motifs

    root = ensure_directory(args.output_dir)
    discovery_dir = ensure_directory(root / "01_discovery")
    filtering_dir = ensure_directory(root / "02_filtering")
    classification_dir = ensure_directory(root / "03_classification")
    motif_dir = ensure_directory(root / "04_motif")

    if args.query_hmm:
        hmm_path = Path(args.query_hmm)
    elif args.seed_alignment:
        hmm_path = build_hmm(args.seed_alignment, discovery_dir / "query.hmm", name=args.hmm_name)
    else:
        raise ValueError("Please provide --query-hmm or --seed-alignment.")

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
        tps_anchor_pattern=args.tps_pattern,
        tps_center_position=args.tps_center_position,
        tps_search_radius=args.tps_search_radius,
        anchor_pattern=args.anchor_pattern,
        flank=args.flank,
        center_position=args.center_position,
        allow_center_fallback=args.allow_center_fallback,
    )
    _warn_legacy_tps_hmms(args.tps_hmm_dir)
    classification_outputs = classify_candidates(
        filtering_outputs["filtered_fasta"],
        args.reference_dir,
        classification_dir,
        hmm_dir=args.tps_hmm_dir,
        candidate_annotations=motif_outputs["motif_summary"],
        top_k=args.top_k,
        tree_neighbors=args.tree_neighbors,
    )
    benchmark_outputs = {}
    if args.expected_fasta:
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
    ):
        for key, value in outputs.items():
            summary_rows.append({"stage": label, "artifact": key, "path": value})
    summary_path = write_tsv(summary_rows, root / "pipeline_summary.tsv")
    logger.info("Pipeline completed. Summary: %s", summary_path)
    return 0


def cmd_compare_fasta(args: argparse.Namespace) -> int:
    from ariadne.benchmark import compare_fasta_sets

    outputs = compare_fasta_sets(args.predicted_fasta, args.expected_fasta, args.output_dir)
    for key, value in outputs.items():
        logger.info("  %-28s %s", key + ":", value)
    return 0


def cmd_visualize(args: argparse.Namespace) -> int:
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


def cmd_filter_coverage_only(args: argparse.Namespace) -> int:
    from ariadne.fasta_utils import read_fasta, write_fasta
    from ariadne.filtering import filter_by_coverage

    records = filter_by_coverage(read_fasta(args.input_fasta), args.min_coverage)
    output_path = write_fasta(records, args.output)
    logger.info("Output: %s", output_path)
    return 0


def cmd_filter_length_only(args: argparse.Namespace) -> int:
    from ariadne.fasta_utils import read_fasta, write_fasta
    from ariadne.filtering import filter_by_length

    records = filter_by_length(read_fasta(args.input_fasta), args.min_length)
    output_path = write_fasta(records, args.output)
    logger.info("Output: %s", output_path)
    return 0


def cmd_dedupe_exact(args: argparse.Namespace) -> int:
    from ariadne.fasta_utils import read_fasta, write_fasta
    from ariadne.filtering import deduplicate_exact

    records = deduplicate_exact(read_fasta(args.input_fasta))
    output_path = write_fasta(records, args.output)
    logger.info("Output: %s", output_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
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

    build_hmm_parser = subparsers.add_parser("build-hmm", help="Build a HMM from an aligned protein FASTA/MSA.")
    build_hmm_parser.add_argument("--alignment", required=True, type=Path)
    build_hmm_parser.add_argument("--output", required=True, type=Path)
    build_hmm_parser.add_argument("--name", default=None)
    build_hmm_parser.set_defaults(func=cmd_build_hmm)

    build_tps_lib = subparsers.add_parser(
        "build-tps-hmm-library",
        help="Build a TPS HMM library directory from one or more aligned FASTA/MSA files.",
    )
    build_tps_lib.add_argument(
        "--alignment",
        nargs="+",
        required=True,
        help="Alignment file path or NAME=PATH. Example: --alignment clade1=align1.fasta clade2=align2.fasta",
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
    classify.add_argument("--tps-hmm-dir", default=_default_tps_hmms(), type=Path, help="Directory containing TPS HMM profiles (*.hmm).")
    classify.add_argument("--candidate-annotations", type=Path, default=None, help="Optional motif/benchmark annotation TSV keyed by sequence_id.")
    classify.add_argument("--top-k", type=int, default=5)
    classify.add_argument("--tree-neighbors", type=int, default=12)
    classify.set_defaults(func=cmd_classify)

    motif = subparsers.add_parser(
        "motif",
        help="Two-stage motif gate: first confirm TPS by DDXXD/E near 125 aa, then cembrene motif judgment (strict by default, optional center fallback).",
    )
    motif.add_argument("--candidates", required=True, type=Path)
    motif.add_argument(
        "--reference-alignment",
        default=_default_alignment_reference(),
        type=Path,
        help="Reference alignment FASTA (defaults to ./Alignment.fasta when available).",
    )
    motif.add_argument("--output-dir", required=True, type=Path)
    motif.add_argument("--tps-pattern", default=r"DD..[DE]")
    motif.add_argument("--tps-center-position", type=int, default=125)
    motif.add_argument("--tps-search-radius", type=int, default=40)
    motif.add_argument("--anchor-pattern", default=r"CFDVL.")
    motif.add_argument("--flank", type=int, default=10)
    motif.add_argument("--center-position", type=int, default=210)
    motif.add_argument(
        "--allow-center-fallback",
        action="store_true",
        help="Enable lenient mode: if anchor motif is absent, fall back to center-position window for cembrene comparison.",
    )
    motif.set_defaults(func=cmd_motif)

    run = subparsers.add_parser("run", help="Execute the four Ariadne modules end-to-end.")
    run.add_argument("--transcriptomes", nargs="+", default=None, type=Path)
    run.add_argument("--protein-folder", type=Path, default=None, help="Preferred input mode: folder of predicted protein FASTA files.")
    run.add_argument(
        "--protein-glob",
        nargs="+",
        default=["*.faa", "*.fa", "*.fasta", "*.pep", "*.prot"],
        help="Glob pattern(s) used to find protein FASTA files under --protein-folder.",
    )
    run.add_argument("--seed-alignment", type=Path, default=None)
    run.add_argument("--query-hmm", type=Path, default=None, help="Use an existing query HMM instead of building from --seed-alignment.")
    run.add_argument("--reference-dir", required=True, type=Path)
    run.add_argument("--output-dir", required=True, type=Path)
    run.add_argument("--expected-fasta", type=Path, default=None, help="Optional benchmark FASTA used to compare against the final cembrene candidate FASTA.")
    run.add_argument(
        "--reference-alignment",
        type=Path,
        default=_default_alignment_reference(),
        help="Reference alignment FASTA for module 4 (defaults to ./Alignment.fasta when available).",
    )
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
    run.add_argument(
        "--allow-center-fallback",
        action="store_true",
        help="Enable lenient motif mode for module 4 when anchor motif is absent.",
    )
    run.add_argument("--tps-hmm-dir", default=_default_tps_hmms(), type=Path, help="Directory containing TPS HMM profiles (*.hmm).")
    run.add_argument("--top-k", type=int, default=5)
    run.add_argument("--tree-neighbors", type=int, default=12)
    run.set_defaults(func=cmd_run)

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
