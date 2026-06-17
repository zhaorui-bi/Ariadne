"""Tests for ariadne.filter — QC filtering and near-duplicate collapsing."""

from __future__ import annotations

from pathlib import Path

from ariadne.filter import (
    _edit_distance_with_limit,
    deduplicate_exact,
    filter_by_coverage,
    filter_by_length,
    filter_candidates,
    near_duplicate,
    record_priority,
)
from ariadne.utils import FastaRecord, write_fasta


class TestEditDistance:
    def test_identical(self):
        assert _edit_distance_with_limit("HELLO", "HELLO", 2) == 0

    def test_single_substitution(self):
        assert _edit_distance_with_limit("HELLO", "HELLP", 2) == 1

    def test_exceeds_limit_returns_limit_plus_one(self):
        assert _edit_distance_with_limit("HELLO", "WORLD", 2) == 3

    def test_length_gap_short_circuits(self):
        assert _edit_distance_with_limit("AAAA", "AAAAAAAA", 1) == 2


class TestNearDuplicate:
    def test_identical_sequences(self):
        assert near_duplicate("A" * 100, "A" * 100, 0.95) is True

    def test_within_threshold(self):
        a = "A" * 100
        b = "A" * 99 + "C"  # one mismatch, allowed_edits = 5
        assert near_duplicate(a, b, 0.95) is True

    def test_outside_threshold(self):
        a = "A" * 100
        b = "A" * 90 + "C" * 10  # ten mismatches, allowed_edits = 5
        assert near_duplicate(a, b, 0.95) is False

    def test_both_empty(self):
        assert near_duplicate("", "", 0.95) is True


class TestBasicFilters:
    def test_record_priority_prefers_higher_coverage(self):
        high = FastaRecord(header="a cov_50", sequence="ACDEF")
        low = FastaRecord(header="b cov_5", sequence="ACDEF")
        assert record_priority(high) > record_priority(low)

    def test_filter_by_coverage_keeps_missing_and_high(self):
        records = [
            FastaRecord(header="a cov_50", sequence="A"),
            FastaRecord(header="b cov_2", sequence="A"),
            FastaRecord(header="c no-cov", sequence="A"),
        ]
        kept = {r.id for r in filter_by_coverage(records, min_coverage=10.0)}
        assert kept == {"a", "c"}

    def test_filter_by_length(self):
        records = [
            FastaRecord(header="short", sequence="A" * 5),
            FastaRecord(header="long", sequence="A" * 50),
        ]
        kept = {r.id for r in filter_by_length(records, min_length=10)}
        assert kept == {"long"}

    def test_deduplicate_exact(self):
        records = [
            FastaRecord(header="a", sequence="ACDEF"),
            FastaRecord(header="b", sequence="ACDEF"),
            FastaRecord(header="c", sequence="MKLMN"),
        ]
        kept = deduplicate_exact(records)
        assert len(kept) == 2


class TestFilterCandidatesEndToEnd:
    def test_pipeline_outputs(self, tmp_path: Path):
        records = [
            FastaRecord(header="keep1 cov_40", sequence="M" + "ACDEFGHIKL" * 4),
            # near-duplicate of keep1 (collapsed away)
            FastaRecord(header="dup1 cov_30", sequence="M" + "ACDEFGHIKL" * 4),
            # distinct representative
            FastaRecord(header="keep2 cov_40", sequence="M" + "MNPQRSTVWY" * 4),
            # removed: too short
            FastaRecord(header="tooshort cov_40", sequence="ACDE"),
            # removed: low coverage
            FastaRecord(header="lowcov cov_2", sequence="M" + "WYACDEFGHI" * 4),
        ]
        input_fasta = tmp_path / "in.faa"
        write_fasta(records, input_fasta)
        out_dir = tmp_path / "out"

        outputs = filter_candidates(
            input_fasta,
            out_dir,
            min_coverage=10.0,
            min_length=20,
            identity_threshold=0.95,
        )

        assert outputs["filtered_fasta"].exists()
        from ariadne.utils import read_fasta

        kept_ids = {r.id for r in read_fasta(outputs["filtered_fasta"])}
        assert kept_ids == {"keep1", "keep2"}

        report = outputs["filter_report"].read_text()
        assert "tooshort" in report and "too_short" in report
        assert "lowcov" in report and "low_coverage" in report
        assert outputs["dedupe_clusters"].exists()
        assert outputs["manual_review"].exists()
