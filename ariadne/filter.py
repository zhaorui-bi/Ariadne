"""Stage 2: candidate quality filtering and near-duplicate collapsing."""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Optional, Union

from ariadne.utils import FastaRecord, parse_coverage, read_fasta, write_fasta, write_tsv
from ariadne.data import load_reference_records

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


def _edit_distance_with_limit(sequence_a: str, sequence_b: str, max_edits: int) -> int:
    """Compute bounded edit distance and stop early once the limit is exceeded."""
    if max_edits < 0:
        return max_edits + 1
    if abs(len(sequence_a) - len(sequence_b)) > max_edits:
        return max_edits + 1
    if len(sequence_a) < len(sequence_b):
        sequence_a, sequence_b = sequence_b, sequence_a
    previous = {index: index for index in range(min(len(sequence_b), max_edits) + 1)}
    for index_a, char_a in enumerate(sequence_a, start=1):
        current: dict[int, int] = {}
        start = max(0, index_a - max_edits)
        end = min(len(sequence_b), index_a + max_edits)
        if start == 0:
            current[0] = index_a
        for index_b in range(max(1, start), end + 1):
            delete_cost = previous.get(index_b, max_edits + 1) + 1
            insert_cost = current.get(index_b - 1, max_edits + 1) + 1
            replace_cost = previous.get(index_b - 1, max_edits + 1) + (char_a != sequence_b[index_b - 1])
            cost = min(delete_cost, insert_cost, replace_cost)
            if cost <= max_edits:
                current[index_b] = cost
        if not current:
            return max_edits + 1
        previous = current
    return previous.get(len(sequence_b), max_edits + 1)


def near_duplicate(sequence_a: str, sequence_b: str, identity_threshold: float) -> bool:
    """Return ``True`` when two sequences meet the configured identity cutoff."""
    if not sequence_a and not sequence_b:
        return True
    max_length = max(len(sequence_a), len(sequence_b))
    allowed_edits = math.floor((1.0 - identity_threshold) * max_length)
    distance = _edit_distance_with_limit(sequence_a, sequence_b, allowed_edits)
    return distance <= allowed_edits


def record_priority(record: FastaRecord) -> tuple[float, int, str]:
    """Rank sequences so higher-coverage and longer representatives are kept first."""
    coverage = parse_coverage(record.header)
    return (coverage if coverage is not None else -1.0, len(record.sequence), record.id)


def filter_by_coverage(records: list[FastaRecord], min_coverage: float) -> list[FastaRecord]:
    """Keep sequences whose parsed coverage is missing or above the threshold."""
    filtered: list[FastaRecord] = []
    for record in records:
        coverage = parse_coverage(record.header)
        if coverage is None or coverage >= min_coverage:
            filtered.append(record)
    return filtered


def filter_by_length(records: list[FastaRecord], min_length: int) -> list[FastaRecord]:
    """Keep sequences longer than or equal to ``min_length`` amino acids."""
    return [record for record in records if len(record.sequence) >= min_length]


def deduplicate_exact(records: list[FastaRecord]) -> list[FastaRecord]:
    """Remove exact duplicate sequences while keeping the first occurrence."""
    seen: dict[str, FastaRecord] = {}
    for record in records:
        seen.setdefault(record.sequence, record)
    return list(seen.values())


def _matching_reference_ids(
    record: FastaRecord,
    reference_records: list[FastaRecord],
    identity_threshold: float,
) -> list[str]:
    """Return reference ids whose sequences meet the configured identity cutoff."""
    matches: list[str] = []
    for reference in reference_records:
        if near_duplicate(record.sequence, reference.sequence, identity_threshold):
            matches.append(reference.id)
    return matches


def filter_candidates(
    input_fasta: PathLike,
    output_dir: PathLike,
    *,
    min_coverage: float = 10.0,
    min_length: int = 300,
    identity_threshold: float = 0.95,
    reference_dir: Optional[PathLike] = None,
) -> dict[str, Path]:
    """Apply basic QC, near-duplicate clustering, and manual-review summaries."""
    records = read_fasta(input_fasta)
    reference_records = load_reference_records(reference_dir) if reference_dir is not None else []
    reference_sequences = {reference.id: reference.sequence for reference in reference_records}
    kept_after_basic_filters: list[FastaRecord] = []
    removed_rows: list[dict[str, object]] = []
    reference_match_rows: list[dict[str, object]] = []

    for record in records:
        reasons: list[str] = []
        coverage = parse_coverage(record.header)
        if coverage is not None and coverage < min_coverage:
            reasons.append("low_coverage")
        if len(record.sequence) < min_length:
            reasons.append("too_short")
        if reasons:
            removed_rows.append(
                {
                    "sequence_id": record.id,
                    "status": "removed",
                    "reason": ",".join(reasons),
                    "coverage": coverage if coverage is not None else "",
                    "length": len(record.sequence),
                }
            )
            continue
        kept_after_basic_filters.append(record)

    representatives: list[FastaRecord] = []
    clusters: list[list[FastaRecord]] = []
    for record in sorted(kept_after_basic_filters, key=record_priority, reverse=True):
        reference_matches = _matching_reference_ids(record, reference_records, identity_threshold)
        for reference_id in reference_matches:
            reference_match_rows.append(
                {
                    "sequence_id": record.id,
                    "reference_id": reference_id,
                    "match_type": "exact" if record.sequence == reference_sequences[reference_id] else "near_duplicate",
                    "coverage": parse_coverage(record.header) or "",
                    "length": len(record.sequence),
                }
            )
        matched_cluster: Optional[list[FastaRecord]] = None
        for cluster in clusters:
            if near_duplicate(record.sequence, cluster[0].sequence, identity_threshold):
                matched_cluster = cluster
                break
        if matched_cluster is None:
            clusters.append([record])
            representatives.append(record)
            removed_rows.append(
                {
                    "sequence_id": record.id,
                    "status": "kept",
                    "reason": "representative",
                    "coverage": parse_coverage(record.header) or "",
                    "length": len(record.sequence),
                }
            )
            continue
        matched_cluster.append(record)
        removed_rows.append(
            {
                "sequence_id": record.id,
                "status": "removed",
                "reason": f"deduplicated_against:{matched_cluster[0].id}",
                "coverage": parse_coverage(record.header) or "",
                "length": len(record.sequence),
            }
        )

    cluster_rows: list[dict[str, object]] = []
    for cluster_index, cluster in enumerate(clusters, start=1):
        representative = cluster[0]
        for member in cluster:
            cluster_rows.append(
                {
                    "cluster": cluster_index,
                    "representative": representative.id,
                    "sequence_id": member.id,
                    "length": len(member.sequence),
                    "coverage": parse_coverage(member.header) or "",
                }
            )

    manual_review_rows = []
    for record in representatives:
        manual_review_rows.append(
            {
                "sequence_id": record.id,
                "length": len(record.sequence),
                "coverage": parse_coverage(record.header) or "",
                "starts_with_m": "yes" if record.sequence.startswith("M") else "no",
                "notes": "manual_visual_check_required",
            }
        )

    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    filtered_fasta_path = write_fasta(representatives, destination / "candidates.filtered.faa")
    report_path = write_tsv(removed_rows, destination / "filter_report.tsv")
    cluster_path = write_tsv(cluster_rows, destination / "dedupe_clusters.tsv")
    review_path = write_tsv(manual_review_rows, destination / "manual_review.tsv")
    reference_match_path = write_tsv(reference_match_rows, destination / "reference_matches.tsv")
    return {
        "filtered_fasta": filtered_fasta_path,
        "filter_report": report_path,
        "dedupe_clusters": cluster_path,
        "manual_review": review_path,
        "reference_matches": reference_match_path,
    }
