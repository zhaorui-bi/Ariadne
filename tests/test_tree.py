"""Tests for ariadne.tree — Newick preview parsing and phylogeny input prep."""

from __future__ import annotations

from pathlib import Path

import pytest

from ariadne.tree import (
    _assign_preview_layout,
    _NewickPreviewParser,
    _preview_category_for_leaf,
    _unique_header,
    prepare_phylogeny_input,
    render_phylogeny_preview,
)
from ariadne.utils import FastaRecord, write_fasta


class TestNewickParser:
    def test_parse_simple_tree(self):
        root = _NewickPreviewParser("(A:1.0,B:2.0)root:0.0;").parse()
        assert len(root.children) == 2
        names = {child.name for child in root.children}
        assert names == {"A", "B"}
        assert root.children[0].length == 1.0

    def test_nested_tree(self):
        root = _NewickPreviewParser("((A:1,B:1):2,C:3);").parse()
        assert len(root.children) == 2
        # one child is an internal node with two leaves, the other is leaf C
        leaf_counts = sorted(len(child.children) for child in root.children)
        assert leaf_counts == [0, 2]

    def test_trailing_content_raises(self):
        with pytest.raises(ValueError):
            _NewickPreviewParser("(A:1,B:1); extra").parse()


class TestPreviewLayout:
    def test_leaf_y_coordinates_are_sequential(self):
        root = _NewickPreviewParser("(A:1,B:1,C:1);").parse()
        leaves = _assign_preview_layout(root)
        assert [leaf.y for leaf in leaves] == [0.0, 1.0, 2.0]


class TestPreviewCategory:
    @pytest.mark.parametrize(
        "name,label",
        [
            ("candidate_seq1", "Candidates"),
            ("ref_coral_x", "Coral refs"),
            ("ref_insect_y", "Insect refs"),
            ("ref_unknownsource_z", "Other refs"),
        ],
    )
    def test_categories(self, name, label):
        assert _preview_category_for_leaf(name)[0] == label


class TestUniqueHeader:
    def test_deduplicates(self):
        seen: set[str] = set()
        first = _unique_header("ref coral x", seen)
        second = _unique_header("ref coral x", seen)
        assert first != second
        assert first in seen and second in seen


class TestRenderPreview:
    def test_render_writes_svg(self, tmp_path: Path):
        tree_file = tmp_path / "t.treefile"
        tree_file.write_text("(candidate_a:1,(ref_coral_b:1,ref_insect_c:1):1);")
        out = render_phylogeny_preview(tree_file, tmp_path / "preview.svg")
        assert out.exists()
        content = out.read_text()
        assert content.startswith("<svg")
        assert "</svg>" in content


class TestPreparePhylogenyInput:
    def test_combines_references_and_candidates(self, tmp_path: Path):
        ref_dir = tmp_path / "refs"
        write_fasta([FastaRecord(header="r1", sequence="ACDEFGHIKL")], ref_dir / "coral.fasta")
        cand = tmp_path / "cand.faa"
        write_fasta([FastaRecord(header="c1", sequence="MNPQRSTVWY")], cand)

        outputs = prepare_phylogeny_input(cand, ref_dir, tmp_path / "out")
        assert outputs["phylogeny_input"].exists()
        assert outputs["sequence_map"].exists()

        from ariadne.utils import read_fasta

        combined = read_fasta(outputs["phylogeny_input"])
        headers = [r.header for r in combined]
        assert any(h.startswith("ref_coral_") for h in headers)
        assert any(h.startswith("candidate_") for h in headers)

    def test_empty_candidates_raises(self, tmp_path: Path):
        ref_dir = tmp_path / "refs"
        write_fasta([FastaRecord(header="r1", sequence="ACDEFGHIKL")], ref_dir / "coral.fasta")
        empty = tmp_path / "empty.faa"
        empty.write_text("")
        with pytest.raises(ValueError):
            prepare_phylogeny_input(empty, ref_dir, tmp_path / "out")
