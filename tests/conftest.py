"""Shared pytest fixtures for the Ariadne test suite."""

from __future__ import annotations

from pathlib import Path

import pytest

from ariadne.utils import FastaRecord, write_fasta

REPO_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture
def repo_root() -> Path:
    """Path to the repository root (for access to bundled data/HMMs)."""
    return REPO_ROOT


@pytest.fixture
def sample_records() -> list[FastaRecord]:
    """A small set of in-memory FASTA records used across tests."""
    return [
        FastaRecord(header="seq1 cov_12.5 sample", sequence="MABCDEFGHIKLMNPQRSTVWY" * 3),
        FastaRecord(header="seq2 cov_3.0", sequence="MABCDEFGHIKLMNPQRSTVWY" * 3),
        FastaRecord(header="seq3", sequence="MKLMNPQRSTVWYACDEFGHIK" * 3),
    ]


@pytest.fixture
def fasta_file(tmp_path: Path, sample_records: list[FastaRecord]) -> Path:
    """Write ``sample_records`` to a temporary FASTA file and return its path."""
    path = tmp_path / "sample.faa"
    write_fasta(sample_records, path)
    return path
