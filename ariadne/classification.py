from __future__ import annotations

import math
import tempfile
from collections import Counter
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pyhmmer

from ariadne.fasta_utils import FastaRecord, ensure_directory, read_fasta, sanitize_newick_name, write_fasta, write_tsv
from ariadne.references import load_reference_records

PathLike = Union[str, Path]


def sorted_hmm_paths(hmm_dir: PathLike) -> list[Path]:
    directory = Path(hmm_dir)
    hmm_paths = list(directory.glob("*.hmm"))
    if not hmm_paths:
        raise FileNotFoundError(
            f"No HMM profiles were found in {directory}. "
            "Please check whether the bundled ariadne/tps_hmm directory exists, or pass --tps-hmm-dir explicitly."
        )
    return sorted(hmm_paths, key=lambda path: (0, int(path.stem)) if path.stem.isdigit() else (1, path.stem))


def _score_records_against_hmms(records: list[FastaRecord], hmm_paths: list[Path], *, missing_score: float = 5.0) -> tuple[np.ndarray, list[str]]:
    if not records:
        return np.zeros((0, len(hmm_paths) + 1), dtype=float), [path.stem for path in hmm_paths] + ["length"]
    alphabet = pyhmmer.easel.Alphabet.amino()
    background = pyhmmer.plan7.Background(alphabet)
    pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)
    name_to_index = {record.id: index for index, record in enumerate(records)}
    matrix = np.full((len(records), len(hmm_paths) + 1), missing_score, dtype=float)
    matrix[:, -1] = np.array([len(record.sequence) for record in records], dtype=float)
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as handle:
        temp_path = Path(handle.name)
    try:
        write_fasta(records, temp_path)
        with pyhmmer.easel.SequenceFile(str(temp_path), digital=True, alphabet=alphabet) as sequence_file:
            sequences = sequence_file.read_block()
        for column_index, hmm_path in enumerate(hmm_paths):
            with pyhmmer.plan7.HMMFile(str(hmm_path)) as hmm_file:
                hmm = hmm_file.read()
            hits = pipeline.search_hmm(hmm, sequences)
            for hit in hits:
                row_index = name_to_index[hit.name]
                matrix[row_index, column_index] = max(matrix[row_index, column_index], float(hit.score))
    finally:
        temp_path.unlink(missing_ok=True)
    column_names = [path.stem for path in hmm_paths] + ["length"]
    return matrix, column_names


def _normalize_matrix(matrix: np.ndarray) -> np.ndarray:
    if matrix.size == 0:
        return matrix.copy()
    means = matrix.mean(axis=0)
    means[means == 0] = 1.0
    return matrix / means


def _pca_coordinates(matrix: np.ndarray, n_components: int = 3) -> tuple[np.ndarray, np.ndarray]:
    if matrix.size == 0:
        return np.zeros((0, n_components), dtype=float), np.zeros((n_components,), dtype=float)
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    if centered.shape[0] == 1:
        coords = np.zeros((1, n_components), dtype=float)
        explained = np.zeros((n_components,), dtype=float)
        return coords, explained
    u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    available = min(n_components, len(singular_values))
    coords = u[:, :available] * singular_values[:available]
    if available < n_components:
        coords = np.pad(coords, ((0, 0), (0, n_components - available)))
    variance = singular_values**2
    explained = variance / variance.sum() if variance.sum() else variance
    if len(explained) < n_components:
        explained = np.pad(explained, (0, n_components - len(explained)))
    return coords, explained[:n_components]


def _distance_matrix(features: np.ndarray) -> np.ndarray:
    if len(features) == 0:
        return np.zeros((0, 0), dtype=float)
    squared = np.sum(features**2, axis=1, keepdims=True)
    distances = squared + squared.T - (2 * features @ features.T)
    np.maximum(distances, 0, out=distances)
    return np.sqrt(distances)


def _upgma_newick(names: list[str], distances: np.ndarray) -> str:
    if len(names) == 1:
        return f"{sanitize_newick_name(names[0])};"

    def key(node_a: int, node_b: int) -> tuple[int, int]:
        return (node_a, node_b) if node_a < node_b else (node_b, node_a)

    active = list(range(len(names)))
    sizes = {index: 1 for index in active}
    heights = {index: 0.0 for index in active}
    labels = {index: sanitize_newick_name(names[index]) for index in active}
    distance_lookup = {
        key(row, column): float(distances[row, column])
        for row in range(len(names))
        for column in range(row + 1, len(names))
    }
    next_index = len(names)

    while len(active) > 1:
        best_pair = None
        best_distance = math.inf
        for i, node_a in enumerate(active):
            for node_b in active[i + 1 :]:
                current_distance = distance_lookup[key(node_a, node_b)]
                if current_distance < best_distance:
                    best_distance = current_distance
                    best_pair = (node_a, node_b)
        assert best_pair is not None
        node_a, node_b = best_pair
        height = best_distance / 2
        label_a = f"{labels[node_a]}:{max(height - heights[node_a], 0):.6f}"
        label_b = f"{labels[node_b]}:{max(height - heights[node_b], 0):.6f}"
        new_label = f"({label_a},{label_b})"
        new_node = next_index
        next_index += 1
        for node_c in active:
            if node_c in (node_a, node_b):
                continue
            distance_lookup[key(new_node, node_c)] = (
                distance_lookup[key(node_a, node_c)] * sizes[node_a]
                + distance_lookup[key(node_b, node_c)] * sizes[node_b]
            ) / (sizes[node_a] + sizes[node_b])
        active = [node for node in active if node not in (node_a, node_b)]
        active.append(new_node)
        sizes[new_node] = sizes[node_a] + sizes[node_b]
        heights[new_node] = height
        labels[new_node] = new_label

    return labels[active[0]] + ";"


def _color_map(sources: list[str]) -> dict[str, str]:
    palette = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#17becf",
        "#bcbd22",
    ]
    mapping: dict[str, str] = {}
    for index, source in enumerate(sorted(set(sources))):
        mapping[source] = palette[index % len(palette)]
    return mapping


def _load_candidate_annotations(path: Optional[PathLike]) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    annotation_path = Path(path)
    if not annotation_path.exists():
        return {}
    import csv

    annotations: dict[str, dict[str, str]] = {}
    with annotation_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sequence_id = (row.get("sequence_id") or "").strip()
            if sequence_id:
                annotations[sequence_id] = {key: value for key, value in row.items() if value}
    return annotations


def _candidate_group(record: FastaRecord) -> str:
    if record.metadata.get("source") != "candidate":
        return f"ref:{record.metadata.get('source', 'unknown')}"
    if record.metadata.get("predicted_cembrene_like") == "yes":
        return "candidate:cembrene_like"
    if record.metadata.get("is_tps") == "yes":
        return "candidate:tps_positive"
    if record.metadata.get("is_tps") == "no":
        return "candidate:tps_negative"
    return "candidate:unlabeled"


def _display_styles(records: list[FastaRecord]) -> dict[str, dict[str, str]]:
    reference_sources = sorted({record.metadata.get("source", "unknown") for record in records if record.metadata.get("source") != "candidate"})
    reference_colors = _color_map(reference_sources)
    styles: dict[str, dict[str, str]] = {}
    for source in reference_sources:
        styles[f"ref:{source}"] = {
            "label": f"reference:{source}",
            "fill": reference_colors[source],
            "shape": "circle",
            "opacity": "0.42",
            "radius": "3.8",
            "stroke": "none",
        }
    styles["candidate:unlabeled"] = {
        "label": "candidate",
        "fill": "#111827",
        "shape": "circle",
        "opacity": "0.88",
        "radius": "6.0",
        "stroke": "#111827",
    }
    styles["candidate:tps_positive"] = {
        "label": "candidate TPS+",
        "fill": "#0f766e",
        "shape": "diamond",
        "opacity": "0.92",
        "radius": "7.4",
        "stroke": "#134e4a",
    }
    styles["candidate:tps_negative"] = {
        "label": "candidate TPS-",
        "fill": "#dc2626",
        "shape": "square",
        "opacity": "0.92",
        "radius": "7.0",
        "stroke": "#7f1d1d",
    }
    styles["candidate:cembrene_like"] = {
        "label": "candidate cembrene-like",
        "fill": "#f59e0b",
        "shape": "triangle",
        "opacity": "0.96",
        "radius": "8.2",
        "stroke": "#92400e",
    }
    return styles


def _marker_svg(x: float, y: float, *, shape: str, radius: float, fill: str, stroke: str, opacity: float, stroke_width: float = 1.0) -> str:
    if shape == "square":
        size = radius * 2
        return (
            f'<rect x="{x - radius:.2f}" y="{y - radius:.2f}" width="{size:.2f}" height="{size:.2f}" '
            f'fill="{fill}" opacity="{opacity}" stroke="{stroke}" stroke-width="{stroke_width}" />'
        )
    if shape == "diamond":
        points = [
            f"{x:.2f},{y - radius:.2f}",
            f"{x + radius:.2f},{y:.2f}",
            f"{x:.2f},{y + radius:.2f}",
            f"{x - radius:.2f},{y:.2f}",
        ]
        return (
            f'<polygon points="{" ".join(points)}" fill="{fill}" opacity="{opacity}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}" />'
        )
    if shape == "triangle":
        points = [
            f"{x:.2f},{y - radius:.2f}",
            f"{x + radius:.2f},{y + radius:.2f}",
            f"{x - radius:.2f},{y + radius:.2f}",
        ]
        return (
            f'<polygon points="{" ".join(points)}" fill="{fill}" opacity="{opacity}" '
            f'stroke="{stroke}" stroke-width="{stroke_width}" />'
        )
    return (
        f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius:.2f}" fill="{fill}" opacity="{opacity}" '
        f'stroke="{stroke}" stroke-width="{stroke_width}" />'
    )


def _render_scatter(records: list[FastaRecord], coords: np.ndarray, output_path: PathLike) -> Path:
    if len(records) == 0:
        return Path(output_path)
    styles = _display_styles(records)
    groups = [_candidate_group(record) for record in records]
    width = 1200
    height = 800
    margin = 80
    xs = coords[:, 0] if len(coords) else np.zeros((len(records),))
    ys = coords[:, 1] if coords.shape[1] > 1 else np.zeros((len(records),))
    min_x, max_x = float(xs.min()), float(xs.max())
    min_y, max_y = float(ys.min()), float(ys.max())
    span_x = max(max_x - min_x, 1e-9)
    span_y = max(max_y - min_y, 1e-9)

    def scale_x(value: float) -> float:
        return margin + ((value - min_x) / span_x) * (width - (2 * margin))

    def scale_y(value: float) -> float:
        return height - margin - ((value - min_y) / span_y) * (height - (2 * margin))

    legend_entries = []
    ordered_groups = [group for group in styles if group in groups]
    for legend_index, group in enumerate(ordered_groups):
        y = margin + (legend_index * 24)
        style = styles[group]
        legend_entries.append(
            _marker_svg(
                width - 220,
                y,
                shape=style["shape"],
                radius=7,
                fill=style["fill"],
                stroke=style["stroke"],
                opacity=float(style["opacity"]),
                stroke_width=1.0,
            )
            + f'<text x="{width - 205}" y="{y + 5}" font-size="14">{style["label"]}</text>'
        )

    points = []
    for record, coordinate in zip(records, coords):
        x = scale_x(float(coordinate[0]))
        y = scale_y(float(coordinate[1] if len(coordinate) > 1 else 0.0))
        group = _candidate_group(record)
        style = styles[group]
        is_candidate = record.metadata.get("source") == "candidate"
        points.append(
            _marker_svg(
                x,
                y,
                shape=style["shape"],
                radius=float(style["radius"]),
                fill=style["fill"],
                stroke=style["stroke"],
                opacity=float(style["opacity"]),
                stroke_width=1.2 if is_candidate else 0.8,
            )
        )
        if is_candidate:
            points.append(
                f'<text x="{x + 10:.2f}" y="{y - 10:.2f}" font-size="14" '
                f'font-family="monospace">{record.id}</text>'
            )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        f'<rect width="{width}" height="{height}" fill="white" />'
        f'<text x="{margin}" y="40" font-size="24" font-weight="bold">Ariadne TPS feature embedding</text>'
        f'<line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" stroke="#999999" />'
        f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height - margin}" stroke="#999999" />'
        + "".join(points)
        + "".join(legend_entries)
        + "</svg>"
    )
    target = Path(output_path)
    target.write_text(svg)
    return target


def _render_3d_sections(records: list[FastaRecord], coords: np.ndarray, output_path: PathLike) -> Path:
    target = Path(output_path)
    if len(records) == 0:
        target.write_text("<svg xmlns='http://www.w3.org/2000/svg' width='1200' height='500'></svg>")
        return target
    if coords.shape[1] < 3:
        coords = np.pad(coords, ((0, 0), (0, 3 - coords.shape[1])))

    styles = _display_styles(records)
    ordered_groups = [group for group in styles if group in {_candidate_group(record) for record in records}]
    projections = [
        (0, 1, "PC1", "PC2"),
        (0, 2, "PC1", "PC3"),
        (1, 2, "PC2", "PC3"),
    ]
    panel_width = 500
    panel_height = 400
    panel_gap = 24
    margin = 56
    legend_height = 90
    width = (panel_width * 3) + (panel_gap * 2) + (margin * 2)
    height = panel_height + legend_height + (margin * 2)

    def panel_scale(values_x: np.ndarray, values_y: np.ndarray):
        min_x, max_x = float(values_x.min()), float(values_x.max())
        min_y, max_y = float(values_y.min()), float(values_y.max())
        span_x = max(max_x - min_x, 1e-9)
        span_y = max(max_y - min_y, 1e-9)

        def sx(value: float, x0: float) -> float:
            return x0 + ((value - min_x) / span_x) * (panel_width - (2 * margin))

        def sy(value: float, y0: float) -> float:
            return y0 + panel_height - margin - ((value - min_y) / span_y) * (panel_height - (2 * margin))

        return sx, sy

    svg_lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="white" />',
        f'<text x="{margin}" y="34" font-size="24" font-weight="bold">Ariadne 3D cluster context (PCA sections)</text>',
        f'<text x="{margin}" y="58" font-size="14" fill="#475569">Reference clades are shown as faint background points; motif-aware candidate states are highlighted as foreground markers.</text>',
    ]

    for panel_index, (dim_x, dim_y, label_x, label_y) in enumerate(projections):
        panel_x = margin + (panel_index * (panel_width + panel_gap))
        panel_y = margin
        xs = coords[:, dim_x]
        ys = coords[:, dim_y]
        scale_x, scale_y = panel_scale(xs, ys)
        svg_lines.append(
            f'<rect x="{panel_x}" y="{panel_y}" width="{panel_width}" height="{panel_height}" fill="#fbfdff" stroke="#dbe3ee" />'
        )
        svg_lines.append(
            f'<line x1="{panel_x + margin}" y1="{panel_y + panel_height - margin}" x2="{panel_x + panel_width - margin}" y2="{panel_y + panel_height - margin}" stroke="#94a3b8" />'
        )
        svg_lines.append(
            f'<line x1="{panel_x + margin}" y1="{panel_y + margin}" x2="{panel_x + margin}" y2="{panel_y + panel_height - margin}" stroke="#94a3b8" />'
        )
        svg_lines.append(
            f'<text x="{panel_x + margin}" y="{panel_y + 24}" font-size="16" font-weight="600">{label_x} vs {label_y}</text>'
        )
        for record, x_val, y_val in zip(records, xs, ys):
            group = _candidate_group(record)
            style = styles[group]
            is_candidate = record.metadata.get("source") == "candidate"
            cx = scale_x(float(x_val), panel_x)
            cy = scale_y(float(y_val), panel_y)
            svg_lines.append(
                _marker_svg(
                    cx,
                    cy,
                    shape=style["shape"],
                    radius=float(style["radius"]) if is_candidate else 3.6,
                    fill=style["fill"],
                    stroke=style["stroke"],
                    opacity=float(style["opacity"]) if is_candidate else 0.38,
                    stroke_width=1.0 if is_candidate else 0.6,
                )
            )
            if is_candidate:
                svg_lines.append(
                    f'<text x="{cx + 8:.2f}" y="{cy - 8:.2f}" font-size="9.5" font-family="monospace" fill="#111827">{record.id}</text>'
                )

    legend_y = margin + panel_height + 30
    legend_x = margin
    for index, group in enumerate(ordered_groups):
        x = legend_x + (index * 180)
        style = styles[group]
        svg_lines.append(
            _marker_svg(
                x,
                legend_y,
                shape=style["shape"],
                radius=7,
                fill=style["fill"],
                stroke=style["stroke"],
                opacity=float(style["opacity"]),
                stroke_width=1.0,
            )
        )
        svg_lines.append(f'<text x="{x + 12}" y="{legend_y + 5}" font-size="14">{style["label"]}</text>')
    svg_lines.append("</svg>")
    target.write_text("".join(svg_lines))
    return target


def classify_candidates(
    candidate_fasta: PathLike,
    reference_dir: PathLike,
    output_dir: PathLike,
    *,
    hmm_dir: PathLike,
    candidate_annotations: Optional[PathLike] = None,
    top_k: int = 5,
    tree_neighbors: int = 12,
) -> dict[str, Path]:
    references = load_reference_records(reference_dir)
    if not references:
        raise ValueError(f"No reference sequences were found in {reference_dir}.")
    candidates = read_fasta(candidate_fasta)
    annotations = _load_candidate_annotations(candidate_annotations)
    for candidate in candidates:
        candidate.metadata["source"] = "candidate"
        candidate.metadata["label"] = "candidate"
        candidate.metadata.update(annotations.get(candidate.id, {}))
        candidate.metadata["candidate_group"] = _candidate_group(candidate)
    all_records = references + candidates
    hmm_paths = sorted_hmm_paths(hmm_dir)
    raw_matrix, column_names = _score_records_against_hmms(all_records, hmm_paths)
    normalized_matrix = _normalize_matrix(raw_matrix)
    coords, explained = _pca_coordinates(normalized_matrix, n_components=3)
    distances = _distance_matrix(normalized_matrix)

    output_root = ensure_directory(output_dir)
    feature_rows = []
    for index, record in enumerate(all_records):
        row = {
            "sequence_id": record.id,
            "source": record.metadata.get("source", "unknown"),
            "label": record.metadata.get("label", ""),
            "candidate_group": record.metadata.get("candidate_group", ""),
        }
        if len(column_names) != len(raw_matrix[index]):
            raise ValueError("Feature column names and feature values have different lengths.")
        for column_name, value in zip(column_names, raw_matrix[index]):
            row[column_name] = round(float(value), 6)
        feature_rows.append(row)
    features_path = write_tsv(feature_rows, output_root / "tps_features.tsv")

    embedding_rows = []
    for index, record in enumerate(all_records):
        embedding_rows.append(
            {
                "sequence_id": record.id,
                "source": record.metadata.get("source", "unknown"),
                "candidate_group": record.metadata.get("candidate_group", ""),
                "pc1": round(float(coords[index, 0]), 6),
                "pc2": round(float(coords[index, 1]), 6),
                "pc3": round(float(coords[index, 2]), 6),
            }
        )
    embedding_path = write_tsv(embedding_rows, output_root / "embedding.tsv")
    variance_path = write_tsv(
        [
            {"component": f"PC{component_index + 1}", "explained_variance": round(float(value), 6)}
            for component_index, value in enumerate(explained, start=0)
        ],
        output_root / "embedding_variance.tsv",
    )

    reference_count = len(references)
    candidate_rows: list[dict[str, object]] = []
    neighbor_rows: list[dict[str, object]] = []
    cluster_context_rows: list[dict[str, object]] = []
    tree_dir = ensure_directory(output_root / "trees")

    for candidate_index in range(reference_count, len(all_records)):
        candidate = all_records[candidate_index]
        reference_distances = [
            (reference_index, float(distances[candidate_index, reference_index]))
            for reference_index in range(reference_count)
        ]
        reference_distances.sort(key=lambda item: item[1])
        top_neighbors = reference_distances[: max(1, min(top_k, len(reference_distances)))]
        predicted_source = Counter(references[neighbor_index].metadata.get("source", "unknown") for neighbor_index, _ in top_neighbors).most_common(1)[0][0]
        confidence = sum(
            1 for neighbor_index, _ in top_neighbors if references[neighbor_index].metadata.get("source", "unknown") == predicted_source
        ) / len(top_neighbors)
        nearest_reference_index, nearest_distance = top_neighbors[0]
        nearest_reference = references[nearest_reference_index]
        candidate_rows.append(
            {
                "sequence_id": candidate.id,
                "candidate_group": candidate.metadata.get("candidate_group", "candidate:unlabeled"),
                "is_tps": candidate.metadata.get("is_tps", ""),
                "predicted_cembrene_like": candidate.metadata.get("predicted_cembrene_like", ""),
                "predicted_source": predicted_source,
                "confidence": round(confidence, 4),
                "nearest_reference": nearest_reference.id,
                "nearest_reference_source": nearest_reference.metadata.get("source", "unknown"),
                "nearest_distance": round(nearest_distance, 6),
            }
        )
        for rank, (neighbor_index, distance) in enumerate(top_neighbors, start=1):
            neighbor = references[neighbor_index]
            neighbor_rows.append(
                {
                    "candidate_id": candidate.id,
                    "rank": rank,
                    "neighbor_id": neighbor.id,
                    "neighbor_source": neighbor.metadata.get("source", "unknown"),
                    "neighbor_label": neighbor.metadata.get("label", ""),
                    "distance": round(distance, 6),
                }
            )
        cluster_context_rows.append(
            {
                "sequence_id": candidate.id,
                "candidate_group": candidate.metadata.get("candidate_group", "candidate:unlabeled"),
                "predicted_source": predicted_source,
                "nearest_reference_source": nearest_reference.metadata.get("source", "unknown"),
                "nearest_distance": round(nearest_distance, 6),
                "pc1": round(float(coords[candidate_index, 0]), 6),
                "pc2": round(float(coords[candidate_index, 1]), 6),
                "pc3": round(float(coords[candidate_index, 2]), 6),
            }
        )

        local_neighbor_indices = [neighbor_index for neighbor_index, _ in reference_distances[: max(2, tree_neighbors)]]
        subset_indices = local_neighbor_indices + [candidate_index]
        subset_names = [all_records[index].id for index in subset_indices]
        subset_distances = distances[np.ix_(subset_indices, subset_indices)]
        tree_path = tree_dir / f"{candidate.id}.nwk"
        tree_path.write_text(_upgma_newick(subset_names, subset_distances) + "\n")

    classification_path = write_tsv(candidate_rows, output_root / "classification.tsv")
    neighbors_path = write_tsv(neighbor_rows, output_root / "nearest_neighbors.tsv")
    cluster_context_path = write_tsv(cluster_context_rows, output_root / "candidate_cluster_context.tsv")
    scatter_path = _render_scatter(all_records, coords, output_root / "embedding.svg")
    sections_path = _render_3d_sections(all_records, coords, output_root / "embedding_3d_sections.svg")
    global_tree_path = output_root / "global_context_tree.nwk"
    global_tree_path.write_text(_upgma_newick([record.id for record in all_records], distances) + "\n")
    assignment_summary_rows = []
    grouped: dict[str, list[float]] = {}
    for row in candidate_rows:
        source = str(row["predicted_source"])
        grouped.setdefault(source, []).append(float(row["confidence"]))
    for source in sorted(grouped):
        confidences = grouped[source]
        assignment_summary_rows.append(
            {
                "predicted_source": source,
                "count": len(confidences),
                "mean_confidence": round(sum(confidences) / len(confidences), 6),
            }
        )
    assignment_summary_path = write_tsv(assignment_summary_rows, output_root / "assignment_summary.tsv")
    return {
        "features": features_path,
        "embedding": embedding_path,
        "variance": variance_path,
        "classification": classification_path,
        "neighbors": neighbors_path,
        "candidate_cluster_context": cluster_context_path,
        "embedding_svg": scatter_path,
        "embedding_3d_sections": sections_path,
        "global_tree": global_tree_path,
        "assignment_summary": assignment_summary_path,
        "tree_dir": tree_dir,
    }
