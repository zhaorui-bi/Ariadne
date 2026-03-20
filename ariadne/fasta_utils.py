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
        return FastaRecord(
            header=header if header is not None else self.header,
            sequence=sequence if sequence is not None else self.sequence,
            metadata=dict(self.metadata),
        )


def ensure_directory(path: PathLike) -> Path:
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def clean_sequence(sequence: str, *, keep_gaps: bool = False) -> str:
    sequence = sequence.replace(" ", "").replace("\r", "").replace("\n", "").upper()
    if keep_gaps:
        sequence = sequence.replace(".", "-")
        return AA_PATTERN.sub("", sequence)
    sequence = sequence.replace("-", "").replace(".", "")
    return AA_PATTERN.sub("", sequence)


def read_fasta(path: PathLike, *, keep_gaps: bool = False) -> list[FastaRecord]:
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
                    records.append(FastaRecord(header=header, sequence=clean_sequence("".join(chunks), keep_gaps=keep_gaps)))
                header = line[1:].strip()
                chunks = []
                continue
            chunks.append(line)
    if header is not None:
        records.append(FastaRecord(header=header, sequence=clean_sequence("".join(chunks), keep_gaps=keep_gaps)))
    return records


def write_fasta(records: Iterable[FastaRecord], path: PathLike, *, width: int = 80) -> Path:
    target = Path(path)
    with target.open("w") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            for offset in range(0, len(record.sequence), width):
                handle.write(record.sequence[offset : offset + width] + "\n")
    return target


def parse_coverage(text: str) -> Optional[float]:
    match = COVERAGE_PATTERN.search(text)
    if match is None:
        return None
    return float(match.group(1))


def ungap(sequence: str) -> str:
    return sequence.replace("-", "").replace(".", "")


def slugify(text: str) -> str:
    compact = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    return compact.strip("_") or "unknown"


def first_existing(*paths: PathLike) -> Optional[Path]:
    for path in paths:
        candidate = Path(path)
        if candidate.exists():
            return candidate
    return None


def write_tsv(rows: Iterable[dict[str, object]], path: PathLike) -> Path:
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
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return safe or "node"


def pad_sequence(sequence: str, length: int, fill: str = "-") -> str:
    if len(sequence) >= length:
        return sequence[:length]
    return sequence + (fill * (length - len(sequence)))


def pairwise_identity(sequence_a: str, sequence_b: str) -> float:
    length = max(len(sequence_a), len(sequence_b))
    if length == 0:
        return 1.0
    matches = 0
    for char_a, char_b in zip(pad_sequence(sequence_a, length), pad_sequence(sequence_b, length)):
        if char_a == char_b:
            matches += 1
    return matches / length
