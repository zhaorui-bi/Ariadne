"""Tests for ariadne.embed — feature math, UPGMA, and classification."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from ariadne.embed import (
    _candidate_group,
    _cluster_label_count,
    _distance_matrix,
    _fmt_tick,
    _nice_ticks,
    _normalize_matrix,
    _pca_coordinates,
    _upgma_newick,
    _zscore_matrix,
    classify_candidates,
    sorted_hmm_paths,
)
from ariadne.tree import _NewickPreviewParser
from ariadne.utils import FastaRecord, read_fasta, write_fasta


class TestFeatureMath:
    def test_normalize_matrix_by_column_mean(self):
        matrix = np.array([[2.0, 0.0], [4.0, 0.0]])
        out = _normalize_matrix(matrix)
        # column 0 mean is 3 -> [2/3, 4/3]; zero-mean column stays zero
        assert np.allclose(out[:, 0], [2 / 3, 4 / 3])
        assert np.allclose(out[:, 1], [0.0, 0.0])

    def test_zscore_constant_column_is_zero(self):
        matrix = np.array([[5.0, 1.0], [5.0, 3.0]])
        out = _zscore_matrix(matrix)
        assert np.allclose(out[:, 0], [0.0, 0.0])
        assert abs(out[:, 1].mean()) < 1e-9

    def test_pca_shape(self):
        matrix = np.random.RandomState(0).rand(6, 4)
        coords, explained = _pca_coordinates(matrix, n_components=3)
        assert coords.shape == (6, 3)
        assert explained.shape == (3,)

    def test_pca_single_row(self):
        coords, explained = _pca_coordinates(np.array([[1.0, 2.0, 3.0]]), n_components=3)
        assert coords.shape == (1, 3)

    def test_distance_matrix(self):
        features = np.array([[0.0, 0.0], [3.0, 4.0]])
        dist = _distance_matrix(features)
        assert np.allclose(dist, [[0.0, 5.0], [5.0, 0.0]])


class TestUpgma:
    def test_single_leaf(self):
        assert _upgma_newick(["A"], np.zeros((1, 1))) == "A;"

    def test_two_leaves_parseable(self):
        names = ["alpha", "beta"]
        distances = np.array([[0.0, 2.0], [2.0, 0.0]])
        newick = _upgma_newick(names, distances)
        root = _NewickPreviewParser(newick).parse()
        assert {child.name for child in root.children} == {"alpha", "beta"}


class TestClusterLabelCount:
    @pytest.mark.parametrize(
        "size,expected",
        [(5, 1), (24, 2), (70, 3), (140, 4), (220, 5)],
    )
    def test_thresholds(self, size, expected):
        assert _cluster_label_count(size) == expected


class TestTicks:
    def test_nice_ticks_span(self):
        ticks = _nice_ticks(0.0, 10.0)
        assert ticks[0] >= 0.0
        assert ticks[-1] <= 10.0 + 1e-6
        assert len(ticks) >= 2

    def test_fmt_tick(self):
        assert _fmt_tick(0.0) == "0"
        assert _fmt_tick(3.0) == "3"


class TestSortedHmmPaths:
    def test_numeric_then_alpha_order(self, tmp_path: Path):
        for name in ("10.hmm", "2.hmm", "coral.hmm"):
            (tmp_path / name).write_text("x")
        stems = [p.stem for p in sorted_hmm_paths(tmp_path)]
        assert stems == ["2", "10", "coral"]

    def test_missing_dir_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError):
            sorted_hmm_paths(tmp_path / "no_hmms_here")


class TestCandidateGroup:
    def test_reference(self):
        rec = FastaRecord(header="r", sequence="A", metadata={"source": "coral"})
        assert _candidate_group(rec) == "ref:coral"

    def test_candidate_ceess(self):
        rec = FastaRecord(
            header="c", sequence="A", metadata={"source": "candidate", "is_ceess_candidate": "yes"}
        )
        assert _candidate_group(rec) == "candidate_ceess"

    def test_candidate_coral_like(self):
        rec = FastaRecord(
            header="c", sequence="A", metadata={"source": "candidate", "is_coral_like": "yes"}
        )
        assert _candidate_group(rec) == "candidate_non_ceess"


class TestClassifyCandidatesIntegration:
    """End-to-end classification using the bundled HMM library (needs pyhmmer)."""

    def test_classify(self, tmp_path: Path, repo_root: Path):
        pytest.importorskip("pyhmmer")
        hmm_dir = repo_root / "ariadne" / "hmm"
        if not any(hmm_dir.glob("*.hmm")):
            pytest.skip("bundled HMM library not available")

        coral = read_fasta(repo_root / "tree" / "coral.fasta")
        insect = read_fasta(repo_root / "tree" / "insect.fasta")
        if len(coral) < 8 or len(insect) < 6:
            pytest.skip("bundled reference FASTAs are smaller than expected")

        ref_dir = tmp_path / "refs"
        write_fasta(coral[:6], ref_dir / "coral.fasta")
        write_fasta(insect[:6], ref_dir / "insect.fasta")

        candidates = tmp_path / "cand.faa"
        write_fasta(coral[6:8], candidates)

        out_dir = tmp_path / "out"
        outputs = classify_candidates(
            candidates,
            ref_dir,
            out_dir,
            hmm_dir=hmm_dir,
            top_k=3,
            tree_neighbors=4,
            ceess_xlsx=None,
        )

        assert outputs["classification"].exists()
        assert outputs["embedding_svg"].read_text().startswith("<svg")
        classification = outputs["classification"].read_text().splitlines()
        # header + 2 candidate rows
        assert len(classification) == 3
        assert "predicted_source" in classification[0]
