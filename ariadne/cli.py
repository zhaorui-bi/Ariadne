from __future__ import annotations

import argparse
from pathlib import Path

from ariadne.classification import classify_candidates
from ariadne.demo import prepare_demo_workspace
from ariadne.discovery import build_hmm, discover_candidates
from ariadne.filtering import deduplicate_exact, filter_by_coverage, filter_by_length, filter_candidates
from ariadne.fasta_utils import ensure_directory, read_fasta, write_fasta, write_tsv
from ariadne.motif import analyze_motifs
from ariadne.references import (
    load_reference_records,
    prepare_coral_reference,
    prepare_extra_reference,
    prepare_insect_reference,
    write_reference_metadata,
)


def _repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _default_coral() -> Path:
    return Path(__file__).resolve().parent / "coralTPS (modified)-cembrene.fasta"


def _default_insect() -> Path:
    return Path(__file__).resolve().parent / "Insecta TPS.xlsx"


def _default_aflp_hmms() -> Path:
    return _repo_root() / "AFLP_finder-main" / "hmm"


def _find_coral_reference(reference_dir: str | Path) -> Path:
    directory = Path(reference_dir)
    direct = directory / "coral.fasta"
    if direct.exists():
        return direct
    for fasta_path in sorted(directory.glob("*.fa*")):
        if "coral" in fasta_path.name.lower():
            return fasta_path
    raise FileNotFoundError(f"Could not locate a coral reference FASTA in {directory}.")


def _parse_extra_reference(spec: str) -> tuple[str, str]:
    if "=" not in spec:
        raise ValueError(f"Expected SOURCE=PATH for extra references, got: {spec}")
    source, path = spec.split("=", 1)
    return source.strip(), path.strip()


def cmd_prepare_references(args: argparse.Namespace) -> int:
    output_dir = ensure_directory(args.output_dir)
    all_records = []
    if args.coral:
        _, coral_records = prepare_coral_reference(args.coral, output_dir, limit=args.coral_limit)
        all_records.extend(coral_records)
    if args.insect_xlsx:
        _, insect_records = prepare_insect_reference(args.insect_xlsx, output_dir, limit=args.insect_limit)
        all_records.extend(insect_records)
    for spec in args.extra_fasta or []:
        source, path = _parse_extra_reference(spec)
        _, extra_records = prepare_extra_reference(path, output_dir, source=source)
        all_records.extend(extra_records)
    metadata_path = write_reference_metadata(all_records, output_dir)
    print(f"Prepared {len(all_records)} reference sequences in {output_dir}")
    print(f"Metadata: {metadata_path}")
    return 0


def cmd_prepare_demo(args: argparse.Namespace) -> int:
    outputs = prepare_demo_workspace(
        args.output_dir,
        coral_alignment_fasta=args.coral,
        insect_xlsx=args.insect_xlsx,
    )
    for key, value in outputs.items():
        print(f"{key}: {value}")
    return 0


def cmd_build_hmm(args: argparse.Namespace) -> int:
    hmm_path = build_hmm(args.alignment, args.output, name=args.name)
    print(f"HMM written to {hmm_path}")
    return 0


def cmd_discover(args: argparse.Namespace) -> int:
    outputs = discover_candidates(
        args.transcriptomes,
        args.hmm,
        args.output_dir,
        min_score=args.min_score,
        max_evalue=args.max_evalue,
    )
    for key, value in outputs.items():
        print(f"{key}: {value}")
    return 0


def cmd_filter(args: argparse.Namespace) -> int:
    outputs = filter_candidates(
        args.input_fasta,
        args.output_dir,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        identity_threshold=args.identity_threshold,
        motif_anchor=args.motif_anchor,
    )
    for key, value in outputs.items():
        print(f"{key}: {value}")
    return 0


def cmd_classify(args: argparse.Namespace) -> int:
    outputs = classify_candidates(
        args.candidates,
        args.reference_dir,
        args.output_dir,
        hmm_dir=args.hmm_dir,
        top_k=args.top_k,
        tree_neighbors=args.tree_neighbors,
    )
    for key, value in outputs.items():
        print(f"{key}: {value}")
    return 0


def cmd_motif(args: argparse.Namespace) -> int:
    outputs = analyze_motifs(
        args.candidates,
        args.coral_reference,
        args.output_dir,
        tps_anchor_pattern=args.tps_pattern,
        tps_center_position=args.tps_center_position,
        tps_search_radius=args.tps_search_radius,
        anchor_pattern=args.anchor_pattern,
        flank=args.flank,
        center_position=args.center_position,
    )
    for key, value in outputs.items():
        print(f"{key}: {value}")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    root = ensure_directory(args.output_dir)
    discovery_dir = ensure_directory(root / "01_discovery")
    filtering_dir = ensure_directory(root / "02_filtering")
    classification_dir = ensure_directory(root / "03_classification")
    motif_dir = ensure_directory(root / "04_motif")

    hmm_path = build_hmm(args.seed_alignment, discovery_dir / "query.hmm", name=args.hmm_name)
    discovery_outputs = discover_candidates(
        args.transcriptomes,
        hmm_path,
        discovery_dir,
        min_score=args.discovery_min_score,
        max_evalue=args.discovery_max_evalue,
    )
    filtering_outputs = filter_candidates(
        discovery_outputs["candidate_proteins"],
        filtering_dir,
        min_coverage=args.min_coverage,
        min_length=args.min_length,
        identity_threshold=args.identity_threshold,
        motif_anchor=args.motif_anchor,
    )
    classification_outputs = classify_candidates(
        filtering_outputs["filtered_fasta"],
        args.reference_dir,
        classification_dir,
        hmm_dir=args.aflp_hmm_dir,
        top_k=args.top_k,
        tree_neighbors=args.tree_neighbors,
    )
    coral_reference = Path(args.coral_reference) if args.coral_reference else _find_coral_reference(args.reference_dir)
    motif_outputs = analyze_motifs(
        filtering_outputs["filtered_fasta"],
        coral_reference,
        motif_dir,
        tps_anchor_pattern=args.tps_pattern,
        tps_center_position=args.tps_center_position,
        tps_search_radius=args.tps_search_radius,
        anchor_pattern=args.anchor_pattern,
        flank=args.flank,
        center_position=args.center_position,
    )
    summary_rows = []
    for label, outputs in (
        ("discovery", discovery_outputs),
        ("filtering", filtering_outputs),
        ("classification", classification_outputs),
        ("motif", motif_outputs),
    ):
        for key, value in outputs.items():
            summary_rows.append({"stage": label, "artifact": key, "path": value})
    summary_path = write_tsv(summary_rows, root / "pipeline_summary.tsv")
    print(f"Pipeline completed. Summary: {summary_path}")
    return 0


def cmd_filter_coverage_only(args: argparse.Namespace) -> int:
    records = filter_by_coverage(read_fasta(args.input_fasta), args.min_coverage)
    output_path = write_fasta(records, args.output)
    print(output_path)
    return 0


def cmd_filter_length_only(args: argparse.Namespace) -> int:
    records = filter_by_length(read_fasta(args.input_fasta), args.min_length)
    output_path = write_fasta(records, args.output)
    print(output_path)
    return 0


def cmd_dedupe_exact(args: argparse.Namespace) -> int:
    records = deduplicate_exact(read_fasta(args.input_fasta))
    output_path = write_fasta(records, args.output)
    print(output_path)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="ariadne", description="Ariadne TPS discovery and annotation pipeline.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_refs = subparsers.add_parser("prepare-references", help="Prepare clean reference FASTA files from the bundled coral and insect resources.")
    prepare_refs.add_argument("--coral", default=_default_coral(), type=Path)
    prepare_refs.add_argument("--coral-limit", type=int, default=None)
    prepare_refs.add_argument("--insect-xlsx", default=_default_insect(), type=Path)
    prepare_refs.add_argument("--insect-limit", type=int, default=None)
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

    discover = subparsers.add_parser("discover", help="Predict ORFs from transcriptomes and search them with a HMM.")
    discover.add_argument("--transcriptomes", nargs="+", required=True, type=Path)
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
    classify.add_argument("--hmm-dir", default=_default_aflp_hmms(), type=Path)
    classify.add_argument("--top-k", type=int, default=5)
    classify.add_argument("--tree-neighbors", type=int, default=12)
    classify.set_defaults(func=cmd_classify)

    motif = subparsers.add_parser("motif", help="First confirm TPS by DDXXD/E near 125 aa, then compare 210 aa motif windows against coral cembrene references.")
    motif.add_argument("--candidates", required=True, type=Path)
    motif.add_argument("--coral-reference", required=True, type=Path)
    motif.add_argument("--output-dir", required=True, type=Path)
    motif.add_argument("--tps-pattern", default=r"DD..[DE]")
    motif.add_argument("--tps-center-position", type=int, default=125)
    motif.add_argument("--tps-search-radius", type=int, default=40)
    motif.add_argument("--anchor-pattern", default=r"CFDVL.")
    motif.add_argument("--flank", type=int, default=10)
    motif.add_argument("--center-position", type=int, default=210)
    motif.set_defaults(func=cmd_motif)

    run = subparsers.add_parser("run", help="Execute the four Ariadne modules end-to-end.")
    run.add_argument("--transcriptomes", nargs="+", required=True, type=Path)
    run.add_argument("--seed-alignment", required=True, type=Path)
    run.add_argument("--reference-dir", required=True, type=Path)
    run.add_argument("--output-dir", required=True, type=Path)
    run.add_argument("--coral-reference", type=Path, default=None)
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
    run.add_argument("--aflp-hmm-dir", default=_default_aflp_hmms(), type=Path)
    run.add_argument("--top-k", type=int, default=5)
    run.add_argument("--tree-neighbors", type=int, default=12)
    run.set_defaults(func=cmd_run)

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


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)
