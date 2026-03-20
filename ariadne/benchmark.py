from __future__ import annotations

from pathlib import Path
from typing import Union

from ariadne.fasta_utils import FastaRecord, ensure_directory, pairwise_identity, read_fasta, write_fasta, write_tsv

PathLike = Union[str, Path]


def _normalize_sequence(sequence: str) -> str:
    return sequence.replace("-", "").replace(".", "").replace("*", "").upper()


def _best_match(record: FastaRecord, targets: list[FastaRecord]) -> tuple[str, float]:
    if not targets:
        return "", 0.0
    best_target = max(targets, key=lambda target: pairwise_identity(record.sequence, target.sequence))
    return best_target.id, pairwise_identity(record.sequence, best_target.sequence)


def compare_fasta_sets(
    predicted_fasta: PathLike,
    expected_fasta: PathLike,
    output_dir: PathLike,
) -> dict[str, Path]:
    predicted = read_fasta(predicted_fasta)
    expected = read_fasta(expected_fasta, keep_gaps=True)
    normalized_predicted_records = [record.clone(sequence=_normalize_sequence(record.sequence)) for record in predicted]
    normalized_expected_records = [record.clone(sequence=_normalize_sequence(record.sequence)) for record in expected]
    normalized_predicted = {
        _normalize_sequence(record.sequence): record
        for record in normalized_predicted_records
        if _normalize_sequence(record.sequence)
    }
    normalized_expected = {
        _normalize_sequence(record.sequence): record
        for record in normalized_expected_records
        if _normalize_sequence(record.sequence)
    }

    shared_sequences = sorted(set(normalized_predicted) & set(normalized_expected))
    predicted_only_sequences = sorted(set(normalized_predicted) - set(normalized_expected))
    expected_only_sequences = sorted(set(normalized_expected) - set(normalized_predicted))
    predicted_ids = {record.id for record in predicted}
    expected_ids = {record.id for record in expected}
    shared_ids = sorted(predicted_ids & expected_ids)

    predicted_rows: list[dict[str, object]] = []
    for original_record, record in zip(predicted, normalized_predicted_records):
        normalized = record.sequence
        exact = normalized in normalized_expected
        best_expected_id, best_identity = _best_match(record, list(normalized_expected.values()))
        predicted_rows.append(
            {
                "direction": "predicted_to_expected",
                "sequence_id": original_record.id,
                "length": len(record.sequence),
                "exact_sequence_match": "yes" if exact else "no",
                "exact_match_id": normalized_expected[normalized].id if exact else "",
                "best_match_id": best_expected_id,
                "best_identity": round(best_identity, 4),
            }
        )

    expected_rows: list[dict[str, object]] = []
    for record in normalized_expected_records:
        normalized = record.sequence
        exact = normalized in normalized_predicted
        best_predicted_id, best_identity = _best_match(record, list(normalized_predicted.values()))
        expected_rows.append(
            {
                "direction": "expected_to_predicted",
                "sequence_id": record.id,
                "length": len(normalized),
                "exact_sequence_match": "yes" if exact else "no",
                "exact_match_id": normalized_predicted[normalized].id if exact else "",
                "best_match_id": best_predicted_id,
                "best_identity": round(best_identity, 4),
            }
        )

    output_root = ensure_directory(output_dir)
    summary_path = write_tsv(
        [
            {"metric": "predicted_count", "value": len(predicted)},
            {"metric": "expected_count", "value": len(expected)},
            {"metric": "shared_exact_sequence_count", "value": len(shared_sequences)},
            {"metric": "shared_sequence_id_count", "value": len(shared_ids)},
            {"metric": "predicted_only_count", "value": len(predicted_only_sequences)},
            {"metric": "expected_only_count", "value": len(expected_only_sequences)},
        ],
        output_root / "benchmark_summary.tsv",
    )
    comparison_path = write_tsv(predicted_rows + expected_rows, output_root / "benchmark_comparison.tsv")
    exact_matches_path = write_tsv(
        [
            {
                "normalized_sequence": sequence,
                "predicted_id": normalized_predicted[sequence].id,
                "expected_id": normalized_expected[sequence].id,
            }
            for sequence in shared_sequences
        ],
        output_root / "exact_matches.tsv",
    )
    predicted_only_path = write_fasta(
        [normalized_predicted[sequence] for sequence in predicted_only_sequences],
        output_root / "only_in_predicted.fasta",
    )
    expected_only_path = write_fasta(
        [normalized_expected[sequence].clone(sequence=sequence) for sequence in expected_only_sequences],
        output_root / "only_in_expected.fasta",
    )
    return {
        "benchmark_summary": summary_path,
        "benchmark_comparison": comparison_path,
        "exact_matches": exact_matches_path,
        "only_in_predicted": predicted_only_path,
        "only_in_expected": expected_only_path,
    }
