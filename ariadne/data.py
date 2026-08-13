"""Reference-data preparation utilities for Ariadne.

Reference FASTA files are the stable boundary between user data and the
classification workflow. This module normalizes heterogeneous sources into a
small directory of FASTA files plus an optional ``metadata.tsv`` table. Later
stages consume that directory through :func:`load_reference_records`.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union

from ariadne.utils import (
    FastaRecord,
    ensure_directory,
    is_fasta_path,
    read_fasta,
    slugify,
    ungap,
    write_fasta,
    write_tsv,
)

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


def _prepare_record(record: FastaRecord, *, source: str, extra: Optional[dict[str, str]] = None) -> FastaRecord:
    """Normalise a reference record and attach source-level metadata."""
    prepared = record.clone(sequence=ungap(record.sequence).replace("*", ""))
    prepared.metadata["source"] = source
    prepared.metadata["header"] = prepared.header
    prepared.metadata["label"] = "CeeSs" if "cembrene" in prepared.header.lower() else "other"
    prepared.metadata["is_cembrene"] = "yes" if "cembrene" in prepared.header.lower() else "no"
    prepared.metadata["is_ceess"] = "yes" if "cembrene" in prepared.header.lower() else "no"
    if extra:
        prepared.metadata.update(extra)
    return prepared


def prepare_coral_reference(
    input_fasta: PathLike,
    output_dir: PathLike,
    *,
    filename: str = "coral.fasta",
    limit: Optional[int] = None,
) -> tuple[Path, list[FastaRecord]]:
    """Prepare coral TPS references from a FASTA alignment or source file.

    Gaps and terminal stop symbols are removed from the output FASTA, while
    source metadata is retained on the in-memory records returned to callers.
    """
    destination = ensure_directory(output_dir)
    records = [_prepare_record(record, source="coral") for record in read_fasta(input_fasta, keep_gaps=True)]
    if limit is not None:
        records = records[:limit]
    output_path = destination / filename
    write_fasta(records, output_path)
    return output_path, records


def prepare_insect_reference(
    input_xlsx: PathLike,
    output_dir: PathLike,
    *,
    filename: str = "insect.fasta",
    sheet_name: str = "Protein Science",
    limit: Optional[int] = None,
) -> tuple[Path, list[FastaRecord]]:
    """Extract insect TPS references from the curated Excel workbook.

    The workbook parser searches for the row containing ``Sequence ID`` and
    ``Sequence`` so modest header offsets in curated spreadsheets do not break
    reference preparation.
    """
    try:
        import openpyxl
    except ImportError as exc:
        raise RuntimeError("openpyxl is required to parse the insect reference workbook.") from exc

    workbook = openpyxl.load_workbook(Path(input_xlsx), read_only=True, data_only=True)
    sheet = workbook[sheet_name]

    header_map: dict[str, int] = {}
    header_row_index = None
    for row_index, row in enumerate(sheet.iter_rows(values_only=True), start=1):
        values = [value.strip() if isinstance(value, str) else value for value in row]
        if "Sequence ID" in values and "Sequence" in values:
            header_map = {str(value): index for index, value in enumerate(values) if value}
            header_row_index = row_index
            break
    if header_row_index is None:
        raise ValueError(f"Could not find the sequence table header in sheet '{sheet_name}'.")

    def optional_cell(row: tuple, column_name: str) -> Optional[object]:
        """Return a cell value by column name, or None when the column is absent/empty."""
        index = header_map.get(column_name)
        if index is None:
            return None
        value = row[index]
        return value if value not in (None, "") else None

    records: list[FastaRecord] = []
    for row in sheet.iter_rows(min_row=header_row_index + 1, values_only=True):
        sequence = row[header_map["Sequence"]]
        if not sequence:
            continue
        sequence_id = optional_cell(row, "Accession") or row[header_map["Sequence ID"]]
        species = optional_cell(row, "Species") or "unknown_species"
        clade = optional_cell(row, "Clade") or "unknown_clade"
        original_id = row[header_map["Sequence ID"]]
        accession = str(sequence_id).strip()
        header = f"{accession} {species} [{clade}]"
        record = FastaRecord(header=header, sequence=str(sequence).strip())
        record.metadata["source"] = "insect"
        record.metadata["header"] = header
        record.metadata["original_id"] = str(original_id).strip()
        record.metadata["species"] = str(species).strip()
        record.metadata["clade"] = str(clade).strip()
        record.metadata["label"] = str(clade).strip()
        record.metadata["is_cembrene"] = "no"
        record.metadata["is_ceess"] = "no"
        records.append(record)
        if limit is not None and len(records) >= limit:
            break

    destination = ensure_directory(output_dir)
    output_path = destination / filename
    write_fasta(records, output_path)
    return output_path, records


def prepare_extra_reference(input_fasta: PathLike, output_dir: PathLike, *, source: str) -> tuple[Path, list[FastaRecord]]:
    """Prepare an extra reference FASTA such as plant, fungal, or bacterial TPSs."""
    destination = ensure_directory(output_dir)
    filename = f"{slugify(source)}.fasta"
    records = [_prepare_record(record, source=source) for record in read_fasta(input_fasta, keep_gaps=True)]
    output_path = destination / filename
    write_fasta(records, output_path)
    return output_path, records


def write_reference_metadata(records: list[FastaRecord], output_dir: PathLike, *, filename: str = "metadata.tsv") -> Path:
    """Write a flattened metadata table for all prepared reference records."""
    rows: list[dict[str, object]] = []
    for record in records:
        row = {"sequence_id": record.id, "header": record.header}
        row.update(record.metadata)
        rows.append(row)
    return write_tsv(rows, ensure_directory(output_dir) / filename)


def load_reference_records(reference_dir: PathLike) -> list[FastaRecord]:
    """Load prepared reference FASTA files and merge metadata when present.

    Each ``*.fa*`` file contributes records whose default ``source`` is derived
    from the filename stem. When ``metadata.tsv`` exists, matching rows override
    or extend the record metadata, preserving labels created by
    ``prepare-references``.
    """
    directory = Path(reference_dir)
    if not directory.exists():
        raise FileNotFoundError(f"Reference directory does not exist: {directory}")
    metadata_map: dict[str, dict[str, str]] = {}
    metadata_path = directory / "metadata.tsv"
    if metadata_path.exists():
        import csv

        with metadata_path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                metadata_map[row["sequence_id"]] = {key: value for key, value in row.items() if value}

    records: list[FastaRecord] = []
    for fasta_path in sorted(path for path in directory.iterdir() if path.is_file() and is_fasta_path(path)):
        source_name = fasta_path.stem.split(".", 1)[0]
        if source_name in {"fungal"}:
            source_name = "fungi"
        for record in read_fasta(fasta_path):
            if not record.sequence:
                continue
            record.metadata.setdefault("source", source_name)
            record.metadata.setdefault("header", record.header)
            if record.id in metadata_map:
                record.metadata.update(metadata_map[record.id])
            records.append(record)
    return records
