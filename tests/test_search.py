"""Tests for ariadne.search — protein file collection and HMM build/search."""

from __future__ import annotations

from pathlib import Path

import pytest

from ariadne.search import build_hmm, collect_protein_files, search_proteins_with_hmm
from ariadne.utils import FastaRecord, write_fasta

# A small, equal-length (already aligned) protein set so build_hmm does not need MAFFT.
_BASE = "MKLVACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTV"  # length 60


def _aligned_records() -> list[FastaRecord]:
    records = [FastaRecord(header="ref0", sequence=_BASE)]
    # introduce a couple of point substitutions while keeping equal length
    for i, (pos, aa) in enumerate([(5, "W"), (10, "Y"), (20, "F"), (40, "M")], start=1):
        seq = list(_BASE)
        seq[pos] = aa
        records.append(FastaRecord(header=f"ref{i}", sequence="".join(seq)))
    return records


class TestCollectProteinFiles:
    def test_collects_recursively_and_sorted(self, tmp_path: Path):
        (tmp_path / "sub").mkdir()
        write_fasta([FastaRecord(header="a", sequence="ACDE")], tmp_path / "b.faa")
        write_fasta([FastaRecord(header="c", sequence="MNPQ")], tmp_path / "sub" / "a.fasta")
        found = collect_protein_files(tmp_path)
        assert len(found) == 2
        assert found == sorted(found)

    def test_missing_directory_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError):
            collect_protein_files(tmp_path / "does-not-exist")


class TestBuildAndSearchHmm:
    def test_build_then_search(self, tmp_path: Path):
        pytest.importorskip("pyhmmer")
        msa = tmp_path / "aln.afa"
        write_fasta(_aligned_records(), msa)

        hmm_path = build_hmm(msa, tmp_path / "profile.hmm", name="test_profile")
        assert hmm_path.exists()
        assert hmm_path.read_text().startswith("HMMER")

        query = tmp_path / "query.faa"
        write_fasta([FastaRecord(header="ref0", sequence=_BASE)], query)
        rows = search_proteins_with_hmm(query, hmm_path)
        assert len(rows) >= 1
        assert "score" in rows[0]
        assert float(rows[0]["score"]) > 0
