"""Shared logging, terminal output, and sequence utility helpers for Ariadne.

This module consolidates two formerly separate concerns:

* Terminal / log formatting (ANSI colours, ASCII banner, ColourFormatter).
* FASTA I/O and lightweight sequence helpers used across every pipeline stage.

Usage
-----
In every module::

    import logging
    logger = logging.getLogger(__name__)

At program entry point::

    from ariadne import utils
    utils.setup_logging(verbose=args.verbose)
    utils.print_banner(__version__)
"""

from __future__ import annotations

import csv
import logging
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional, Union

PathLike = Union[str, Path]

# ---------------------------------------------------------------------------
# ANSI colour helpers
# ---------------------------------------------------------------------------


def _supports_colour() -> bool:
    """Return True when *stderr* is a colour-capable TTY."""
    return hasattr(sys.stderr, "isatty") and sys.stderr.isatty()


class _Colours:
    """Container for ANSI escape sequences; empty strings when unsupported."""

    def __init__(self) -> None:
        on = _supports_colour()
        self.RESET   = "\033[0m"   if on else ""
        self.BOLD    = "\033[1m"   if on else ""
        self.DIM     = "\033[2m"   if on else ""
        self.RED     = "\033[91m"  if on else ""
        self.YELLOW  = "\033[93m"  if on else ""
        self.GREEN   = "\033[92m"  if on else ""
        self.CYAN    = "\033[96m"  if on else ""
        self.BLUE    = "\033[94m"  if on else ""
        self.MAGENTA = "\033[95m"  if on else ""


C = _Colours()

# ---------------------------------------------------------------------------
# ASCII banner
# ---------------------------------------------------------------------------

_BANNER_ART = (
    r"    _    ____  ___    _    ____  _   _ _____  ",
    r"   / \  |  _ \|_ _|  / \  |  _ \| \ | | ____|",
    r"  / _ \ | |_) || |  / _ \ | | | |  \| |  _|  ",
    r" / ___ \|  _ < | | / ___ \| |_| | |\  | |___ ",
    r"/_/   \_\_| \_\___/_/   \_\____/|_| \_|_____|",
)

_SUBTITLE = "Terpene Synthase Discovery & Annotation Pipeline"
_BOX_TOTAL = 62


def _box_rule() -> str:
    return "─" * (_BOX_TOTAL - 2)


def _box_line(content: str = "", visible_len: int | None = None) -> str:
    if visible_len is None:
        visible_len = len(content)
    padding = " " * max(0, _BOX_TOTAL - 4 - visible_len)
    return f"│ {content}{padding} │"


def _render_banner(version: str) -> str:
    rule = _box_rule()
    rows: list[str] = [f"╭{rule}╮", _box_line()]
    for art_line in _BANNER_ART:
        coloured = f"  {C.CYAN}{art_line}{C.RESET}"
        rows.append(_box_line(coloured, visible_len=2 + len(art_line)))
    rows.append(_box_line())
    subtitle_coloured = f"  {C.DIM}{_SUBTITLE}{C.RESET}"
    rows.append(_box_line(subtitle_coloured, visible_len=2 + len(_SUBTITLE)))
    ver_str = f"version {version}"
    ver_coloured = f"  {C.BOLD}{ver_str}{C.RESET}"
    rows.append(_box_line(ver_coloured, visible_len=2 + len(ver_str)))
    rows.append(_box_line())
    rows.append(f"╰{rule}╯")
    return "\n".join(rows)


def print_banner(version: str) -> None:
    """Print the Ariadne banner to *stderr*."""
    print(_render_banner(version), file=sys.stderr)
    print(file=sys.stderr)


# ---------------------------------------------------------------------------
# Coloured log formatter
# ---------------------------------------------------------------------------


class _ColourFormatter(logging.Formatter):
    """Log formatter that adds ANSI colours to level names and module names."""

    _LEVEL_LABEL: dict[int, str] = {}

    def __init__(self) -> None:
        super().__init__()
        self._LEVEL_LABEL = {
            logging.DEBUG:    C.DIM    + "DEBUG" + C.RESET,
            logging.INFO:     C.CYAN   + "INFO " + C.RESET,
            logging.WARNING:  C.YELLOW + "WARN " + C.RESET,
            logging.ERROR:    C.RED    + "ERROR" + C.RESET,
            logging.CRITICAL: C.BOLD + C.RED + "CRIT " + C.RESET,
        }

    def format(self, record: logging.LogRecord) -> str:
        level  = self._LEVEL_LABEL.get(record.levelno, record.levelname)
        module = C.MAGENTA + record.name.split(".")[-1] + C.RESET
        return f"[{level}] [{module}] {record.getMessage()}"


def setup_logging(verbose: bool = False) -> None:
    """Configure the *ariadne* logger hierarchy."""
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_ColourFormatter())
    root = logging.getLogger("ariadne")
    root.setLevel(level)
    root.handlers.clear()
    root.addHandler(handler)
    root.propagate = False


# ---------------------------------------------------------------------------
# FASTA record
# ---------------------------------------------------------------------------

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
        return FastaRecord(
            header=header if header is not None else self.header,
            sequence=sequence if sequence is not None else self.sequence,
            metadata=dict(self.metadata),
        )


# ---------------------------------------------------------------------------
# FASTA I/O
# ---------------------------------------------------------------------------


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


_logger = logging.getLogger(__name__)


def _record_from_header_and_chunks(header: str, chunks: list[str], *, keep_gaps: bool) -> FastaRecord:
    if chunks:
        sequence = clean_sequence("".join(chunks), keep_gaps=keep_gaps)
        return FastaRecord(header=header, sequence=sequence)
    parts = header.split(maxsplit=1)
    if len(parts) == 2:
        candidate_header, embedded = parts
        cleaned_embedded = clean_sequence(embedded, keep_gaps=keep_gaps)
        if len(cleaned_embedded) >= 30 and cleaned_embedded == clean_sequence(cleaned_embedded, keep_gaps=keep_gaps):
            _logger.warning("Recovered malformed FASTA record with header-embedded sequence: %s", candidate_header)
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


logger = logging.getLogger(__name__)
