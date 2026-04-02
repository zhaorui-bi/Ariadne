"""Sequence-based phylogeny utilities powered by MAFFT and IQ-TREE."""

from __future__ import annotations

import html
import logging
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

from ariadne.utils import FastaRecord, ensure_directory, slugify, write_fasta, write_tsv
from ariadne.data import load_reference_records

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


@dataclass
class _PreviewNode:
    """Minimal Newick node used to render a static SVG tree preview."""

    name: str | None = None
    length: float = 0.0
    children: list["_PreviewNode"] = field(default_factory=list)
    x: float = 0.0
    y: float = 0.0

    @property
    def is_leaf(self) -> bool:
        """Return True when the preview node is a leaf."""
        return not self.children


class _NewickPreviewParser:
    """Parse the subset of Newick emitted by IQ-TREE for static preview rendering."""

    def __init__(self, text: str) -> None:
        self.text = text.strip()
        self.index = 0

    def parse(self) -> _PreviewNode:
        """Parse the full Newick string into a preview tree."""
        node = self._parse_node()
        self._consume_whitespace()
        if self._peek() == ";":
            self.index += 1
        self._consume_whitespace()
        if self.index != len(self.text):
            raise ValueError(f"Unexpected trailing Newick content at position {self.index}")
        return node

    def _parse_node(self) -> _PreviewNode:
        self._consume_whitespace()
        if self._peek() == "(":
            self.index += 1
            children: list[_PreviewNode] = []
            while True:
                children.append(self._parse_node())
                self._consume_whitespace()
                token = self._peek()
                if token == ",":
                    self.index += 1
                    continue
                if token == ")":
                    self.index += 1
                    break
                raise ValueError(f"Unexpected token {token!r} at position {self.index}")
            name = self._parse_name()
            length = self._parse_length()
            return _PreviewNode(name=name or None, length=length, children=children)

        name = self._parse_name()
        length = self._parse_length()
        return _PreviewNode(name=name or None, length=length)

    def _parse_name(self) -> str:
        self._consume_whitespace()
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ",():;":
            self.index += 1
        return self.text[start:self.index].strip()

    def _parse_length(self) -> float:
        self._consume_whitespace()
        if self._peek() != ":":
            return 0.0
        self.index += 1
        start = self.index
        while self.index < len(self.text) and self.text[self.index] not in ",() ;":
            self.index += 1
        raw = self.text[start:self.index].strip()
        return float(raw) if raw else 0.0

    def _peek(self) -> str:
        if self.index >= len(self.text):
            return ""
        return self.text[self.index]

    def _consume_whitespace(self) -> None:
        while self.index < len(self.text) and self.text[self.index].isspace():
            self.index += 1


def _assign_preview_layout(root: _PreviewNode) -> list[_PreviewNode]:
    """Assign phylogram x/y coordinates and return leaves in plotting order."""
    leaves: list[_PreviewNode] = []

    def walk(node: _PreviewNode, parent_x: float) -> None:
        node.x = parent_x + max(node.length, 0.0)
        for child in node.children:
            walk(child, node.x)
        if node.is_leaf:
            node.y = float(len(leaves))
            leaves.append(node)
        else:
            node.y = sum(child.y for child in node.children) / len(node.children)

    walk(root, 0.0)
    return leaves


def _preview_category_for_leaf(name: str) -> tuple[str, str]:
    """Map tree tip ids to preview labels and colors."""
    if name.startswith("candidate_"):
        return "Candidates", "#F97316"
    if name.startswith("ref_coral_"):
        return "Coral refs", "#0F766E"
    if name.startswith("ref_bacteria_"):
        return "Bacteria refs", "#334155"
    if name.startswith("ref_insect_"):
        return "Insect refs", "#7C3AED"
    if name.startswith("ref_plant_"):
        return "Plant refs", "#16A34A"
    if name.startswith("ref_fungi_"):
        return "Fungi refs", "#DC2626"
    return "Other refs", "#64748B"


def _gather_preview_counts(leaves: list[_PreviewNode]) -> list[tuple[str, str, int]]:
    """Count leaf categories for the preview legend."""
    counts: dict[str, tuple[str, int]] = {}
    for leaf in leaves:
        label, color = _preview_category_for_leaf(leaf.name or "unknown")
        _, current = counts.get(label, (color, 0))
        counts[label] = (color, current + 1)
    order = [
        "Candidates",
        "Coral refs",
        "Bacteria refs",
        "Insect refs",
        "Plant refs",
        "Fungi refs",
        "Other refs",
    ]
    return [(label, counts[label][0], counts[label][1]) for label in order if label in counts]


def _collect_preview_edges(node: _PreviewNode) -> list[tuple[_PreviewNode, _PreviewNode]]:
    """Return parent-child edges for the preview renderer."""
    edges: list[tuple[_PreviewNode, _PreviewNode]] = []
    stack = [node]
    while stack:
        current = stack.pop()
        for child in current.children:
            edges.append((current, child))
            stack.append(child)
    return edges


def render_phylogeny_preview(treefile: PathLike, output_svg: PathLike) -> Path:
    """Render a compact SVG preview directly from an IQ-TREE Newick file."""
    tree_path = Path(treefile)
    output_path = Path(output_svg)
    root = _NewickPreviewParser(tree_path.read_text(encoding="utf-8")).parse()
    leaves = _assign_preview_layout(root)

    width = 1600
    height = min(1180, max(980, 220 + len(leaves)))
    left_margin = 120
    right_margin = 180
    bottom_margin = 70
    legend = _gather_preview_counts(leaves)
    total_leaves = len(leaves)
    candidate_count = next((count for label, _, count in legend if label == "Candidates"), 0)
    legend_columns = 3
    legend_rows = max(1, (len(legend) + legend_columns - 1) // legend_columns)
    legend_row_height = 24
    legend_top = 174
    header_y = 36
    header_height = max(148, legend_top - header_y + legend_rows * legend_row_height + 26)
    top_margin = header_y + header_height + 28
    plot_width = width - left_margin - right_margin
    plot_height = height - top_margin - bottom_margin
    max_x = max((leaf.x for leaf in leaves), default=1.0) or 1.0
    y_step = plot_height / max(len(leaves) - 1, 1)

    def px(node_x: float) -> float:
        return left_margin + (node_x / max_x) * plot_width

    def py(node_y: float) -> float:
        return top_margin + node_y * y_step

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '  <title id="title">Ariadne phylogeny preview</title>',
        '  <desc id="desc">Compact phylogram rendered from the IQ-TREE result for downstream preview and reporting.</desc>',
        '  <rect width="100%" height="100%" fill="#F8FAFC"/>',
        f'  <rect x="40" y="{header_y}" width="1520" height="{header_height}" rx="28" fill="#FFFFFF" stroke="#E2E8F0"/>',
        '  <text x="72" y="78" fill="#0F172A" font-family="Helvetica, Arial, sans-serif" font-size="34" font-weight="700">Ariadne Phylogeny Preview</text>',
        f'  <text x="72" y="112" fill="#475569" font-family="Helvetica, Arial, sans-serif" font-size="19">Rendered from {html.escape(tree_path.name)} with {total_leaves} leaves and {candidate_count} candidate sequences.</text>',
        '  <rect x="72" y="126" width="170" height="10" rx="5" fill="#0F766E" opacity="0.9"/>',
        '  <rect x="250" y="126" width="170" height="10" rx="5" fill="#F97316" opacity="0.9"/>',
        '  <text x="72" y="150" fill="#334155" font-family="Helvetica, Arial, sans-serif" font-size="15">Coral reference branches and candidate placements derived from the current phylogeny stage.</text>',
    ]

    lines.append(
        f'  <text x="72" y="{legend_top - 10}" fill="#0F172A" font-family="Helvetica, Arial, sans-serif" font-size="15" font-weight="700">Reference / candidate groups</text>'
    )
    legend_x = 72
    legend_y = legend_top + 6
    legend_col_width = 280
    for idx, (label, color, count) in enumerate(legend):
        row = idx // legend_columns
        col = idx % legend_columns
        x = legend_x + col * legend_col_width
        y = legend_y + row * legend_row_height
        lines.append(
            f'  <circle cx="{x}" cy="{y - 4}" r="6" fill="{color}" opacity="0.95"/>'
            f'<text x="{x + 14}" y="{y}" fill="#334155" font-family="Helvetica, Arial, sans-serif" font-size="13">{html.escape(label)} ({count})</text>'
        )

    lines.append(f'  <rect x="{left_margin}" y="{top_margin}" width="{plot_width}" height="{plot_height}" fill="#FFFFFF" stroke="#E2E8F0" rx="24"/>')

    for parent, child in _collect_preview_edges(root):
        x0 = px(parent.x)
        x1 = px(child.x)
        y0 = py(parent.y)
        y1 = py(child.y)
        lines.append(f'  <line x1="{x0:.2f}" y1="{y1:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" stroke="#CBD5E1" stroke-width="1.1"/>')
        lines.append(f'  <line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x0:.2f}" y2="{y1:.2f}" stroke="#CBD5E1" stroke-width="1.1"/>')

    for leaf in leaves:
        _, color = _preview_category_for_leaf(leaf.name or "")
        x = px(leaf.x)
        y = py(leaf.y)
        tip_x = width - right_margin + 20
        lines.append(f'  <line x1="{x:.2f}" y1="{y:.2f}" x2="{tip_x:.2f}" y2="{y:.2f}" stroke="{color}" stroke-width="1.4" opacity="0.55"/>')
        lines.append(f'  <circle cx="{tip_x:.2f}" cy="{y:.2f}" r="1.9" fill="{color}" opacity="0.9"/>')

    lines.extend(
        [
            f'  <text x="{left_margin}" y="{height - 28}" fill="#64748B" font-family="Helvetica, Arial, sans-serif" font-size="15">Preview view for Ariadne outputs. Use {html.escape(tree_path.name)} and iqtree.iqtree for full-resolution phylogenetic interpretation.</text>',
            '</svg>',
        ]
    )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def _resolve_binary(preferred: Optional[str], *fallbacks: str) -> str:
    """Resolve an external executable name from preferred and fallback options."""
    candidates = []
    if preferred:
        candidates.append(preferred)
    candidates.extend(fallbacks)
    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved:
            return resolved
    expected = ", ".join(candidates) if candidates else "<none>"
    raise FileNotFoundError(f"Could not find required executable. Tried: {expected}")


def _unique_header(base: str, seen: set[str]) -> str:
    """Generate a unique FASTA header token suitable for MAFFT/IQ-TREE."""
    header = slugify(base)
    if header not in seen:
        seen.add(header)
        return header
    index = 2
    while f"{header}_{index}" in seen:
        index += 1
    unique = f"{header}_{index}"
    seen.add(unique)
    return unique


def prepare_phylogeny_input(
    candidate_fasta: PathLike,
    reference_dir: PathLike,
    output_dir: PathLike,
) -> dict[str, Path]:
    """Combine references and candidates into a clean FASTA for tree building."""
    from ariadne.utils import read_fasta

    output_root = ensure_directory(output_dir)
    reference_records = load_reference_records(reference_dir)
    candidate_records = read_fasta(candidate_fasta)
    if not reference_records:
        raise ValueError(f"No reference sequences were found in {reference_dir}.")
    if not candidate_records:
        raise ValueError(f"No candidate sequences were found in {candidate_fasta}.")

    combined_records: list[FastaRecord] = []
    mapping_rows: list[dict[str, object]] = []
    seen_headers: set[str] = set()

    for record in reference_records:
        source = record.metadata.get("source", "reference")
        header = _unique_header(f"ref_{source}_{record.id}", seen_headers)
        combined_records.append(FastaRecord(header=header, sequence=record.sequence, metadata=dict(record.metadata)))
        mapping_rows.append(
            {
                "tree_id": header,
                "record_type": "reference",
                "source": source,
                "original_id": record.id,
                "original_header": record.header,
            }
        )

    for record in candidate_records:
        header = _unique_header(f"candidate_{record.id}", seen_headers)
        combined_records.append(FastaRecord(header=header, sequence=record.sequence, metadata=dict(record.metadata)))
        mapping_rows.append(
            {
                "tree_id": header,
                "record_type": "candidate",
                "source": "candidate",
                "original_id": record.id,
                "original_header": record.header,
            }
        )

    combined_fasta = write_fasta(combined_records, output_root / "phylogeny_input.fasta")
    mapping_tsv = write_tsv(mapping_rows, output_root / "phylogeny_sequence_map.tsv")
    return {
        "phylogeny_input": combined_fasta,
        "sequence_map": mapping_tsv,
    }


def run_mafft(
    input_fasta: PathLike,
    output_alignment: PathLike,
    *,
    mafft_bin: Optional[str] = None,
    mode: str = "--auto",
) -> Path:
    """Run MAFFT and write the aligned FASTA output."""
    mafft = _resolve_binary(mafft_bin, "mafft")
    output_path = Path(output_alignment)
    command = [mafft, mode, str(input_fasta)]
    logger.info("Running MAFFT: %s", " ".join(command))
    with output_path.open("w") as handle:
        subprocess.run(command, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)
    return output_path


def run_iqtree(
    alignment_fasta: PathLike,
    output_dir: PathLike,
    *,
    iqtree_bin: Optional[str] = None,
    model: str = "LG",
    threads: str = "AUTO",
    bootstrap: Optional[int] = None,
    fast: bool = True,
) -> dict[str, Path]:
    """Run IQ-TREE on an alignment and collect the main output files."""
    iqtree = _resolve_binary(iqtree_bin, "iqtree2", "iqtree")
    output_root = ensure_directory(output_dir)
    prefix = output_root / "iqtree"
    command = [
        iqtree,
        "-s",
        str(alignment_fasta),
        "--prefix",
        str(prefix),
        "-m",
        model,
        "-T",
        str(threads),
    ]
    if fast:
        command.append("--fast")
    if bootstrap is not None and bootstrap > 0:
        command.extend(["-B", str(bootstrap)])
    logger.info("Running IQ-TREE: %s", " ".join(command))
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    outputs = {
        "iqtree_prefix": prefix,
        "iqtree_treefile": prefix.with_suffix(".treefile"),
        "iqtree_report": prefix.with_suffix(".iqtree"),
        "iqtree_log": prefix.with_suffix(".log"),
    }
    optional = {
        "iqtree_consensus_tree": prefix.with_suffix(".contree"),
        "iqtree_checkpoint": prefix.with_suffix(".ckp.gz"),
    }
    for key, path in optional.items():
        if path.exists():
            outputs[key] = path
    return outputs


def build_phylogeny(
    candidate_fasta: PathLike,
    reference_dir: PathLike,
    output_dir: PathLike,
    *,
    mafft_bin: Optional[str] = None,
    mafft_mode: str = "--auto",
    iqtree_bin: Optional[str] = None,
    iqtree_model: str = "LG",
    iqtree_threads: str = "AUTO",
    iqtree_bootstrap: Optional[int] = None,
    iqtree_fast: bool = True,
) -> dict[str, Path]:
    """Build a sequence alignment and phylogenetic tree from references plus candidates."""
    output_root = ensure_directory(output_dir)
    prepared = prepare_phylogeny_input(candidate_fasta, reference_dir, output_root)
    alignment_path = run_mafft(
        prepared["phylogeny_input"],
        output_root / "phylogeny_alignment.fasta",
        mafft_bin=mafft_bin,
        mode=mafft_mode,
    )
    iqtree_outputs = run_iqtree(
        alignment_path,
        output_root,
        iqtree_bin=iqtree_bin,
        model=iqtree_model,
        threads=iqtree_threads,
        bootstrap=iqtree_bootstrap,
        fast=iqtree_fast,
    )
    outputs = dict(prepared)
    outputs["phylogeny_alignment"] = alignment_path
    outputs.update(iqtree_outputs)
    preview_path = render_phylogeny_preview(
        outputs["iqtree_treefile"],
        output_root / "phylogeny_preview.svg",
    )
    outputs["phylogeny_preview"] = preview_path
    return outputs
