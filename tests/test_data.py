"""Tests for ariadne.data — reference preparation and loading."""

from __future__ import annotations

from pathlib import Path

import pytest

from ariadne.data import (
    _prepare_record,
    load_reference_records,
    prepare_extra_reference,
    prepare_insect_reference,
    write_reference_metadata,
)
from ariadne.utils import FastaRecord, write_fasta


def _make_insect_workbook(path: Path, header: list[str], rows: list[list[object]]) -> None:
    """Write a minimal insect-reference workbook for parser tests."""
    openpyxl = pytest.importorskip("openpyxl")
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "Protein Science"
    sheet.append(header)
    for row in rows:
        sheet.append(row)
    workbook.save(path)


class TestPrepareRecord:
    def test_cembrene_header_marked_as_ceess(self):
        record = FastaRecord(header="x1 cembrene synthase", sequence="AC-DE.F")
        prepared = _prepare_record(record, source="coral")
        assert prepared.sequence == "ACDEF"  # gaps removed
        assert prepared.metadata["source"] == "coral"
        assert prepared.metadata["label"] == "CeeSs"
        assert prepared.metadata["is_ceess"] == "yes"

    def test_non_cembrene_header(self):
        record = FastaRecord(header="x2 some other synthase", sequence="ACDEF")
        prepared = _prepare_record(record, source="plant")
        assert prepared.metadata["label"] == "other"
        assert prepared.metadata["is_ceess"] == "no"


class TestPrepareExtraReference:
    def test_writes_named_fasta_with_source_metadata(self, tmp_path: Path):
        src = tmp_path / "src.fasta"
        write_fasta(
            [
                FastaRecord(header="p1", sequence="ACDEFGHIKL"),
                FastaRecord(header="p2", sequence="MNPQRSTVWY"),
            ],
            src,
        )
        out_dir = tmp_path / "refs"
        path, records = prepare_extra_reference(src, out_dir, source="plant")
        assert path.name == "plant.fasta"
        assert path.exists()
        assert len(records) == 2
        assert all(r.metadata["source"] == "plant" for r in records)


class TestMetadataRoundTrip:
    def test_write_and_load(self, tmp_path: Path):
        refs = tmp_path / "refs"
        write_fasta([FastaRecord(header="c1 cembrene", sequence="ACDEFGHIKL")], refs / "coral.fasta")
        coral_records = [
            _prepare_record(r, source="coral")
            for r in [FastaRecord(header="c1 cembrene", sequence="ACDEFGHIKL")]
        ]
        meta_path = write_reference_metadata(coral_records, refs)
        assert meta_path.exists()

        loaded = load_reference_records(refs)
        assert len(loaded) == 1
        record = loaded[0]
        assert record.metadata["source"] == "coral"
        # metadata.tsv merge should carry through the CeeSs label
        assert record.metadata.get("is_ceess") == "yes"

    def test_source_inferred_from_filename(self, tmp_path: Path):
        refs = tmp_path / "refs"
        write_fasta([FastaRecord(header="i1", sequence="ACDEFGHIKL")], refs / "insect.fasta")
        loaded = load_reference_records(refs)
        assert loaded[0].metadata["source"] == "insect"


class TestPrepareInsectReference:
    def test_missing_optional_columns_fall_back(self, tmp_path: Path):
        # Workbook lacks Accession / Species / Clade columns. The parser must
        # fall back to placeholders instead of silently grabbing the last column.
        xlsx = tmp_path / "insect.xlsx"
        _make_insect_workbook(
            xlsx,
            header=["Sequence ID", "Sequence"],
            rows=[["seqA", "ACDEFGHIKLMNPQRSTVWY"]],
        )
        _, records = prepare_insect_reference(xlsx, tmp_path / "out")
        assert len(records) == 1
        record = records[0]
        assert record.metadata["species"] == "unknown_species"
        assert record.metadata["clade"] == "unknown_clade"
        # accession falls back to the Sequence ID column, NOT the last column
        assert record.id == "seqA"

    def test_present_columns_are_read(self, tmp_path: Path):
        xlsx = tmp_path / "insect.xlsx"
        _make_insect_workbook(
            xlsx,
            header=["Sequence ID", "Accession", "Species", "Clade", "Sequence"],
            rows=[["id1", "ACC1", "CoralSp", "CladeX", "ACDEFGHIKLMNPQRSTVWY"]],
        )
        _, records = prepare_insect_reference(xlsx, tmp_path / "out")
        record = records[0]
        assert record.id == "ACC1"
        assert record.metadata["species"] == "CoralSp"
        assert record.metadata["clade"] == "CladeX"
