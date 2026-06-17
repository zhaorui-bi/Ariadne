"""Tests for ariadne.utils — FASTA I/O, sequence helpers, logging, TSV."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from ariadne import utils
from ariadne.utils import (
    FastaRecord,
    clean_sequence,
    ensure_directory,
    first_existing,
    pad_sequence,
    pairwise_identity,
    parse_coverage,
    read_fasta,
    sanitize_newick_name,
    slugify,
    ungap,
    write_fasta,
    write_tsv,
)


class TestCleanSequence:
    def test_uppercase_and_whitespace_removed(self):
        assert clean_sequence("ac def\n") == "ACDEF"

    def test_gaps_removed_by_default(self):
        assert clean_sequence("AC-DE.F") == "ACDEF"

    def test_gaps_kept_and_normalised(self):
        assert clean_sequence("AC-DE.F", keep_gaps=True) == "AC-DE-F"

    def test_non_amino_acid_characters_dropped(self):
        assert clean_sequence("AC1DE2") == "ACDE"


class TestFastaRecord:
    def test_id_and_description(self):
        record = FastaRecord(header="seq1 some description here", sequence="ACDEF")
        assert record.id == "seq1"
        assert record.description == "some description here"

    def test_description_empty_when_no_space(self):
        assert FastaRecord(header="seq1", sequence="ACDEF").description == ""

    def test_clone_is_independent(self):
        record = FastaRecord(header="a", sequence="ACD", metadata={"k": "v"})
        clone = record.clone(sequence="EFG")
        clone.metadata["k"] = "changed"
        assert record.sequence == "ACD"
        assert record.metadata["k"] == "v"
        assert clone.sequence == "EFG"


class TestFastaIO:
    def test_roundtrip(self, tmp_path: Path):
        records = [
            FastaRecord(header="s1 desc", sequence="ACDEFGHIKL"),
            FastaRecord(header="s2", sequence="MNPQRSTVWY"),
        ]
        path = tmp_path / "x.faa"
        write_fasta(records, path)
        loaded = read_fasta(path)
        assert [r.header for r in loaded] == ["s1 desc", "s2"]
        assert [r.sequence for r in loaded] == ["ACDEFGHIKL", "MNPQRSTVWY"]

    def test_multiline_sequence_is_joined(self, tmp_path: Path):
        path = tmp_path / "ml.faa"
        path.write_text(">m1\nACDE\nFGHI\nKLMN\n")
        loaded = read_fasta(path)
        assert len(loaded) == 1
        assert loaded[0].sequence == "ACDEFGHIKLMN"

    def test_write_wraps_long_lines(self, tmp_path: Path):
        record = FastaRecord(header="long", sequence="A" * 200)
        path = tmp_path / "long.faa"
        write_fasta([record], path, width=80)
        body_lines = [ln for ln in path.read_text().splitlines() if not ln.startswith(">")]
        assert all(len(ln) <= 80 for ln in body_lines)
        assert "".join(body_lines) == "A" * 200


class TestSequenceHelpers:
    @pytest.mark.parametrize(
        "text,expected",
        [("NODE_1_cov_12.5", 12.5), ("cov_3", 3.0), ("no coverage here", None)],
    )
    def test_parse_coverage(self, text, expected):
        assert parse_coverage(text) == expected

    def test_ungap(self):
        assert ungap("AC-DE.F") == "ACDEF"

    def test_slugify(self):
        assert slugify("hello world!") == "hello_world"
        assert slugify("   ") == "unknown"

    def test_sanitize_newick_name(self):
        assert sanitize_newick_name("a (b):c,d") == "a_b_c_d"
        # non-empty but unsafe characters collapse to a single underscore
        assert sanitize_newick_name("###") == "_"
        # only a fully empty result falls back to the placeholder
        assert sanitize_newick_name("") == "node"

    def test_pad_sequence(self):
        assert pad_sequence("ABC", 5) == "ABC--"
        assert pad_sequence("ABCDEF", 3) == "ABC"

    def test_pairwise_identity(self):
        assert pairwise_identity("ABCDE", "ABCDE") == 1.0
        assert pairwise_identity("", "") == 1.0
        assert pairwise_identity("ABCD", "ABXD") == 0.75


class TestFilesystemHelpers:
    def test_ensure_directory_creates(self, tmp_path: Path):
        target = tmp_path / "a" / "b"
        out = ensure_directory(target)
        assert out.is_dir()

    def test_first_existing(self, tmp_path: Path):
        present = tmp_path / "here.txt"
        present.write_text("x")
        assert first_existing(tmp_path / "missing", present) == present
        assert first_existing(tmp_path / "nope") is None


class TestWriteTsv:
    def test_empty_rows_writes_empty_file(self, tmp_path: Path):
        path = write_tsv([], tmp_path / "empty.tsv")
        assert path.read_text() == ""

    def test_header_and_values(self, tmp_path: Path):
        rows = [{"a": 1, "b": 2}, {"a": 3, "b": 4}]
        path = write_tsv(rows, tmp_path / "t.tsv")
        lines = path.read_text().splitlines()
        assert lines[0] == "a\tb"
        assert lines[1] == "1\t2"

    def test_ragged_rows_union_of_keys(self, tmp_path: Path):
        rows = [{"a": 1}, {"a": 2, "b": 3}]
        path = write_tsv(rows, tmp_path / "r.tsv")
        header = path.read_text().splitlines()[0]
        assert header == "a\tb"


class TestLogging:
    def test_log_file_is_written(self, tmp_path: Path):
        log_path = tmp_path / "nested" / "run.log"
        utils.setup_logging(verbose=True, log_file=log_path)
        logging.getLogger("ariadne.test").info("hello-from-test")
        for handler in logging.getLogger("ariadne").handlers:
            handler.flush()
        assert log_path.exists()
        assert "hello-from-test" in log_path.read_text()
        # Reset handlers so the file handle does not leak into other tests.
        utils.setup_logging(verbose=False)
