from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional, Tuple, Union

from ariadne.fasta_utils import FastaRecord, pad_sequence, pairwise_identity, read_fasta, write_fasta, write_tsv

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


def _match_pattern_near_position(
    sequence: str,
    pattern: str,
    *,
    center_position: int,
    search_radius: int,
) -> Optional[Tuple[int, str]]:
    if not sequence:
        return None
    search_start = max(0, center_position - search_radius)
    search_end = min(len(sequence), center_position + search_radius + 1)
    search_region = sequence[search_start:search_end]
    matches = list(re.finditer(pattern, search_region))
    if not matches:
        return None
    best_match = min(
        matches,
        key=lambda match: abs((search_start + match.start()) - center_position),
    )
    start = search_start + best_match.start()
    return start, sequence[start : search_start + best_match.end()]


def _anchor_window(
    sequence: str,
    anchor_pattern: str,
    flank: int,
    *,
    center_position: int,
    fallback_span: int = 6,
    allow_center_fallback: bool = True,
) -> Optional[Tuple[str, int, str, str]]:
    match = re.search(anchor_pattern, sequence)
    if match is not None:
        start = match.start()
        end = match.end()
        method = "anchor_match"
    else:
        if not allow_center_fallback:
            return None
        if len(sequence) < center_position + 1:
            return None
        start = max(0, center_position)
        end = min(len(sequence), start + fallback_span)
        method = "center_fallback"
    window_start = max(0, start - flank)
    window_end = min(len(sequence), end + flank)
    window = sequence[window_start:window_end]
    window = pad_sequence(window, (flank * 2) + (end - start))
    return window, start, sequence[start:end], method


def _aa_color(amino_acid: str) -> str:
    amino_acid = amino_acid.upper()
    if amino_acid in {"D", "E"}:
        return "#d62728"
    if amino_acid in {"K", "R", "H"}:
        return "#9467bd"
    if amino_acid in {"S", "T", "N", "Q", "Y", "C"}:
        return "#1f77b4"
    if amino_acid in {"A", "V", "L", "I", "M", "F", "W"}:
        return "#2ca02c"
    if amino_acid == "G":
        return "#ff7f0e"
    if amino_acid == "P":
        return "#8c564b"
    if amino_acid == "-":
        return "#bbbbbb"
    return "#333333"


def _render_motif_svg(rows: list[tuple[str, str, str]], flank: int, output_path: PathLike) -> Path:
    cell_width = 16
    cell_height = 24
    label_width = 320
    width = label_width + (len(rows[0][2]) * cell_width) + 60 if rows else 640
    height = 70 + (len(rows) * cell_height) + 40
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="white" />',
        '<text x="20" y="30" font-size="22" font-weight="bold">Ariadne motif window view</text>',
    ]
    if rows:
        motif_width = len(rows[0][2])
        anchor_start = flank * cell_width
        anchor_width = max(cell_width, (motif_width - (2 * flank)) * cell_width)
        lines.append(
            f'<rect x="{label_width + anchor_start}" y="40" width="{anchor_width}" height="{height - 60}" '
            f'fill="#fff6d5" />'
        )
    for row_index, (group, label, motif_window) in enumerate(rows, start=0):
        y = 70 + (row_index * cell_height)
        lines.append(f'<text x="20" y="{y + 17}" font-size="13" font-family="monospace">{group}: {label}</text>')
        for column_index, amino_acid in enumerate(motif_window):
            x = label_width + (column_index * cell_width)
            lines.append(
                f'<rect x="{x}" y="{y}" width="{cell_width}" height="{cell_height - 4}" '
                f'fill="{_aa_color(amino_acid)}" opacity="0.22" />'
            )
            lines.append(
                f'<text x="{x + 4}" y="{y + 17}" font-size="13" font-family="monospace" '
                f'fill="{_aa_color(amino_acid)}">{amino_acid}</text>'
            )
    lines.append("</svg>")
    target = Path(output_path)
    target.write_text("".join(lines))
    return target


def analyze_motifs(
    candidate_fasta: PathLike,
    reference_alignment_fasta: PathLike,
    output_dir: PathLike,
    *,
    tps_anchor_pattern: str = r"DD..[DE]",
    tps_center_position: int = 125,
    tps_search_radius: int = 40,
    anchor_pattern: str = r"CFDVL.",
    flank: int = 10,
    center_position: int = 210,
    allow_center_fallback: bool = False,
) -> dict[str, Path]:
    candidates = read_fasta(candidate_fasta)
    reference_records = read_fasta(reference_alignment_fasta)
    cembrene_refs = [
        record
        for record in reference_records
        if "cembrene" in record.header.lower()
        and _match_pattern_near_position(
            record.sequence,
            tps_anchor_pattern,
            center_position=tps_center_position,
            search_radius=tps_search_radius,
        )
        is not None
    ]
    non_cembrene_refs = [
        record
        for record in reference_records
        if "cembrene" not in record.header.lower()
        and _match_pattern_near_position(
            record.sequence,
            tps_anchor_pattern,
            center_position=tps_center_position,
            search_radius=tps_search_radius,
        )
        is not None
    ]

    cembrene_windows = []
    non_cembrene_windows = []
    for record in cembrene_refs:
        anchor = _anchor_window(
            record.sequence,
            anchor_pattern,
            flank,
            center_position=center_position,
            allow_center_fallback=allow_center_fallback,
        )
        if anchor is not None:
            cembrene_windows.append((record, anchor[0], anchor[1], anchor[2]))
    for record in non_cembrene_refs:
        anchor = _anchor_window(
            record.sequence,
            anchor_pattern,
            flank,
            center_position=center_position,
            allow_center_fallback=allow_center_fallback,
        )
        if anchor is not None:
            non_cembrene_windows.append((record, anchor[0], anchor[1], anchor[2]))

    summary_rows = []
    motif_records: list[FastaRecord] = []
    cembrene_candidate_records: list[FastaRecord] = []
    display_rows: list[tuple[str, str, str]] = []

    for candidate in candidates:
        tps_match = _match_pattern_near_position(
            candidate.sequence,
            tps_anchor_pattern,
            center_position=tps_center_position,
            search_radius=tps_search_radius,
        )
        if tps_match is None:
            summary_rows.append(
                {
                    "sequence_id": candidate.id,
                    "is_tps": "no",
                    "tps_motif_found": "no",
                    "tps_motif_position": "",
                    "tps_motif_variant": "",
                    "anchor_found": "no",
                    "anchor_position": "",
                    "anchor_variant": "",
                    "cembrene_anchor_found": "no",
                    "cembrene_window_method": "not_assessed",
                    "cembrene_anchor_position": "",
                    "cembrene_anchor_variant": "",
                    "motif_window": "",
                    "best_cembrene_match": "",
                    "best_cembrene_identity": 0.0,
                    "best_non_cembrene_match": "",
                    "best_non_cembrene_identity": 0.0,
                    "predicted_cembrene_like": "no",
                }
            )
            continue

        tps_position, tps_variant = tps_match
        anchor = _anchor_window(
            candidate.sequence,
            anchor_pattern,
            flank,
            center_position=center_position,
            allow_center_fallback=allow_center_fallback,
        )
        if anchor is None:
            summary_rows.append(
                {
                    "sequence_id": candidate.id,
                    "is_tps": "yes",
                    "tps_motif_found": "yes",
                    "tps_motif_position": tps_position,
                    "tps_motif_variant": tps_variant,
                    "anchor_found": "no",
                    "anchor_position": "",
                    "anchor_variant": "",
                    "cembrene_anchor_found": "no",
                    "cembrene_window_method": "anchor_not_found",
                    "cembrene_anchor_position": "",
                    "cembrene_anchor_variant": "",
                    "motif_window": "",
                    "best_cembrene_match": "",
                    "best_cembrene_identity": 0.0,
                    "best_non_cembrene_match": "",
                    "best_non_cembrene_identity": 0.0,
                    "predicted_cembrene_like": "no",
                }
            )
            continue

        motif_window, motif_position, anchor_variant, anchor_method = anchor
        best_cembrene = max(
            ((record, pairwise_identity(motif_window, window)) for record, window, _, _ in cembrene_windows),
            key=lambda item: item[1],
            default=(None, 0.0),
        )
        best_non_cembrene = max(
            ((record, pairwise_identity(motif_window, window)) for record, window, _, _ in non_cembrene_windows),
            key=lambda item: item[1],
            default=(None, 0.0),
        )
        if best_cembrene[0] is None and best_non_cembrene[0] is None:
            predicted = "undetermined"
        elif best_cembrene[0] is None:
            predicted = "no"
        elif best_non_cembrene[0] is None:
            predicted = "yes"
        else:
            predicted = "yes" if best_cembrene[1] >= best_non_cembrene[1] else "no"
        summary_rows.append(
            {
                "sequence_id": candidate.id,
                "is_tps": "yes",
                "tps_motif_found": "yes",
                "tps_motif_position": tps_position,
                "tps_motif_variant": tps_variant,
                "anchor_found": "yes" if anchor_method == "anchor_match" else "fallback_window",
                "anchor_position": motif_position,
                "anchor_variant": anchor_variant,
                "cembrene_anchor_found": "yes" if anchor_method == "anchor_match" else "no",
                "cembrene_window_method": anchor_method,
                "cembrene_anchor_position": motif_position,
                "cembrene_anchor_variant": anchor_variant,
                "motif_window": motif_window,
                "best_cembrene_match": best_cembrene[0].id if best_cembrene[0] else "",
                "best_cembrene_identity": round(best_cembrene[1], 4),
                "best_non_cembrene_match": best_non_cembrene[0].id if best_non_cembrene[0] else "",
                "best_non_cembrene_identity": round(best_non_cembrene[1], 4),
                "predicted_cembrene_like": predicted,
            }
        )
        if predicted == "yes":
            cembrene_candidate_records.append(candidate.clone())
        motif_records.append(FastaRecord(header=f"{candidate.id}|motif_window", sequence=motif_window))
        display_rows.append(("candidate", candidate.id, motif_window))
        if best_cembrene[0] is not None:
            best_window = next(window for record, window, _, _ in cembrene_windows if record.id == best_cembrene[0].id)
            display_rows.append(("cembrene_ref", best_cembrene[0].id, best_window))
        if best_non_cembrene[0] is not None:
            best_window = next(window for record, window, _, _ in non_cembrene_windows if record.id == best_non_cembrene[0].id)
            display_rows.append(("other_ref", best_non_cembrene[0].id, best_window))

    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = write_tsv(summary_rows, output_root / "motif_summary.tsv")
    motif_fasta_path = write_fasta(motif_records, output_root / "motif_windows.fasta")
    motif_svg_path = _render_motif_svg(display_rows, flank, output_root / "motif_windows.svg")
    cembrene_candidates_path = write_tsv(
        [row for row in summary_rows if row.get("predicted_cembrene_like") == "yes"],
        output_root / "cembrene_candidates.tsv",
    )
    cembrene_fasta_path = write_fasta(cembrene_candidate_records, output_root / "cembrene_candidates.fasta")
    return {
        "motif_summary": summary_path,
        "motif_windows": motif_fasta_path,
        "motif_svg": motif_svg_path,
        "cembrene_candidates": cembrene_candidates_path,
        "cembrene_fasta": cembrene_fasta_path,
    }
