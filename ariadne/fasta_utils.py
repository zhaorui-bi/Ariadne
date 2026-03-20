"""Shared FASTA parsing and lightweight sequence utility helpers.

These helpers are intentionally dependency-light because almost every stage in
the pipeline needs them. The functions here standardise how Ariadne reads,
normalises, writes, and compares protein FASTA records.
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional, Union

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
AA_PATTERN = re.compile(r"[^A-Z*\-]")
COVERAGE_PATTERN = re.compile(r"cov_([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)


@dataclass
class FastaRecord:
    """Simple in-memory FASTA record used across the pipeline."""

    header: str
    sequence: str
    metadata: dict[str, str] = field(default_factory=dict)

    @property
    def id(self) -> str:
        return self.header.split()[0]

    @property
    def description(self) -> str:
        if " " not in self.header:
            return ""
        return self.header.split(" ", 1)[1]

    def clone(self, *, header: Optional[str] = None, sequence: Optional[str] = None) -> "FastaRecord":
        """Return a shallow copy while preserving metadata."""
        return FastaRecord(
            header=header if header is not None else self.header,
            sequence=sequence if sequence is not None else self.sequence,
            metadata=dict(self.metadata),
        )


def ensure_directory(path: PathLike) -> Path:
    """Create *path* when needed and return it as a :class:`Path`."""
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def clean_sequence(sequence: str, *, keep_gaps: bool = False) -> str:
    """Normalise FASTA sequence text to Ariadne's canonical representation."""
    sequence = sequence.replace(" ", "").replace("\r", "").replace("\n", "").upper()
    if keep_gaps:
        sequence = sequence.replace(".", "-")
        return AA_PATTERN.sub("", sequence)
    sequence = sequence.replace("-", "").replace(".", "")
    return AA_PATTERN.sub("", sequence)


def _record_from_header_and_chunks(header: str, chunks: list[str], *, keep_gaps: bool) -> FastaRecord:
    """Build a FASTA record, rescuing malformed header-embedded sequences when possible."""
    if chunks:
        sequence = clean_sequence("".join(chunks), keep_gaps=keep_gaps)
        return FastaRecord(header=header, sequence=sequence)

    parts = header.split(maxsplit=1)
    if len(parts) == 2:
        candidate_header, embedded = parts
        cleaned_embedded = clean_sequence(embedded, keep_gaps=keep_gaps)
        if len(cleaned_embedded) >= 30 and cleaned_embedded == clean_sequence(cleaned_embedded, keep_gaps=keep_gaps):
            logger.warning("Recovered malformed FASTA record with header-embedded sequence: %s", candidate_header)
            return FastaRecord(header=candidate_header, sequence=cleaned_embedded)
    return FastaRecord(header=header, sequence="")


def read_fasta(path: PathLike, *, keep_gaps: bool = False) -> list[FastaRecord]:
    """Parse a FASTA file into :class:`FastaRecord` objects."""
    records: list[FastaRecord] = []
    header: Optional[str] = None
    chunks: list[str] = []
    with Path(path).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(_record_from_header_and_chunks(header, chunks, keep_gaps=keep_gaps))
                header = line[1:].strip()
                chunks = []
                continue
            chunks.append(line)
    if header is not None:
        records.append(_record_from_header_and_chunks(header, chunks, keep_gaps=keep_gaps))
    return records


def write_fasta(records: Iterable[FastaRecord], path: PathLike, *, width: int = 80) -> Path:
    """Write FASTA records using wrapped sequence lines."""
    target = Path(path)
    with target.open("w") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            for offset in range(0, len(record.sequence), width):
                handle.write(record.sequence[offset : offset + width] + "\n")
    return target


def parse_coverage(text: str) -> Optional[float]:
    """Extract ``cov_<value>`` style coverage metadata from a FASTA header."""
    match = COVERAGE_PATTERN.search(text)
    if match is None:
        return None
    return float(match.group(1))


def ungap(sequence: str) -> str:
    """Remove common alignment gap symbols from a sequence string."""
    return sequence.replace("-", "").replace(".", "")


def slugify(text: str) -> str:
    """Convert free text into a filename-safe identifier."""
    compact = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    return compact.strip("_") or "unknown"


def first_existing(*paths: PathLike) -> Optional[Path]:
    """Return the first existing path from a list of candidates."""
    for path in paths:
        candidate = Path(path)
        if candidate.exists():
            return candidate
    return None


def write_tsv(rows: Iterable[dict[str, object]], path: PathLike) -> Path:
    """Write a list of dictionaries to a tab-separated table."""
    rows = list(rows)
    target = Path(path)
    if not rows:
        target.write_text("")
        return target
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with target.open("w", newline="") as handle:
        import csv

        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return target


def sanitize_newick_name(name: str) -> str:
    """Replace unsupported characters so a label is safe in a Newick tree."""
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return safe or "node"


def pad_sequence(sequence: str, length: int, fill: str = "-") -> str:
    """Right-pad a sequence to a fixed length for simple comparisons."""
    if len(sequence) >= length:
        return sequence[:length]
    return sequence + (fill * (length - len(sequence)))


def pairwise_identity(sequence_a: str, sequence_b: str) -> float:
    """Compute a naive position-wise identity after right-padding sequences."""
    length = max(len(sequence_a), len(sequence_b))
    if length == 0:
        return 1.0
    matches = 0
    for char_a, char_b in zip(pad_sequence(sequence_a, length), pad_sequence(sequence_b, length)):
        if char_a == char_b:
            matches += 1
    return matches / length
