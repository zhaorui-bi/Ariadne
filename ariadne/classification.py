"""Stage 3: classify candidates in TPS HMM feature space.

Each sequence is converted into a vector of HMM scores, normalised, embedded,
and compared against reference sequences to infer a likely source group such as
coral, insect, plant, fungal, or bacterial TPS collections.
"""

from __future__ import annotations

import logging
import math
import tempfile
from collections import Counter
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pyhmmer
from sklearn.cluster import KMeans
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from ariadne.fasta_utils import FastaRecord, ensure_directory, read_fasta, sanitize_newick_name, write_fasta, write_tsv
from ariadne.references import load_reference_records
from ariadne.esm_type import DEFAULT_ESM_MODEL_NAME

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


def sorted_hmm_paths(hmm_dir: PathLike) -> list[Path]:
    """Return HMM paths in a stable order for reproducible feature columns."""
    directory = Path(hmm_dir)
    hmm_paths = list(directory.glob("*.hmm"))
    if not hmm_paths:
        raise FileNotFoundError(
            f"No HMM profiles were found in {directory}. "
            "Please check whether the bundled ariadne/tps_hmm directory exists, or pass --tps-hmm-dir explicitly."
        )
    return sorted(hmm_paths, key=lambda path: (0, int(path.stem)) if path.stem.isdigit() else (1, path.stem))


def _score_records_against_hmms(records: list[FastaRecord], hmm_paths: list[Path], *, missing_score: float = 5.0) -> tuple[np.ndarray, list[str]]:
    """Score all records against all HMMs and build the raw feature matrix."""
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
    """Normalise each feature column by its mean value."""
    if matrix.size == 0:
        return matrix.copy()
    means = matrix.mean(axis=0)
    means[means == 0] = 1.0
    return matrix / means


def _zscore_matrix(matrix: np.ndarray) -> np.ndarray:
    """Z-score each column so embedding axes are not dominated by one feature."""
    if matrix.size == 0:
        return matrix.copy()
    means = matrix.mean(axis=0, keepdims=True)
    stds = matrix.std(axis=0, keepdims=True)
    stds[stds == 0] = 1.0
    return (matrix - means) / stds


def _pca_coordinates(matrix: np.ndarray, n_components: int = 3) -> tuple[np.ndarray, np.ndarray]:
    """Compute PCA coordinates with NumPy's SVD implementation."""
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


def _embedding_group_label(record: FastaRecord) -> str:
    """Return the label used by supervised embedding to separate display groups."""
    if record.metadata.get("source") == "candidate":
        return _candidate_group(record)
    return f"ref:{record.metadata.get('source', 'unknown')}"


def _cluster_label_count(size: int) -> int:
    """Choose a conservative subcluster count for one large reference subset."""
    if size >= 220:
        return 5
    if size >= 140:
        return 4
    if size >= 70:
        return 3
    if size >= 24:
        return 2
    return 1


def _stable_kmeans_labels(matrix: np.ndarray, *, prefix: str) -> list[str]:
    """Partition a subset into a small number of deterministic KMeans groups."""
    if len(matrix) == 0:
        return []
    n_clusters = _cluster_label_count(len(matrix))
    if n_clusters <= 1:
        return [f"{prefix}_1"] * len(matrix)
    kmeans = KMeans(n_clusters=n_clusters, random_state=0, n_init=20)
    assignments = kmeans.fit_predict(matrix)
    order = sorted(
        range(n_clusters),
        key=lambda cluster_index: tuple(np.round(kmeans.cluster_centers_[cluster_index][: min(3, matrix.shape[1])], 6)),
    )
    remap = {cluster_index: rank + 1 for rank, cluster_index in enumerate(order)}
    return [f"{prefix}_{remap[int(cluster)]}" for cluster in assignments]


def _embedding_group_labels(matrix: np.ndarray, records: list[FastaRecord]) -> list[str]:
    """Build subclade-aware labels for the embedding stage.

    The main goal is to stop the large coral reference collection from behaving
    like one monolithic class. We therefore split coral references into small,
    deterministic feature-space subclusters while keeping candidates grouped
    together as one foreground layer.
    """
    labels = [_embedding_group_label(record) for record in records]
    coral_indices: list[int] = []
    for index, record in enumerate(records):
        if record.metadata.get("source") == "candidate":
            continue
        if record.metadata.get("source") != "coral":
            continue
        coral_indices.append(index)
    if coral_indices:
        coral_matrix = matrix[coral_indices]
        coral_labels = _stable_kmeans_labels(coral_matrix, prefix="ref:coral:tps_subclade")
        for index, label in zip(coral_indices, coral_labels):
            labels[index] = label
    return labels


def _deterministic_direction(index: int, dimensions: int) -> np.ndarray:
    """Return a deterministic unit vector used when two centroids fully overlap."""
    base = np.zeros((dimensions,), dtype=float)
    base[index % dimensions] = 1.0
    if dimensions > 1:
        base[(index + 1) % dimensions] = 0.5
    norm = np.linalg.norm(base)
    return base if norm == 0 else base / norm


def _spread_group_coords(coords: np.ndarray, labels: list[str], *, within_scale: float = 0.55) -> np.ndarray:
    """Spread group centroids apart while preserving local within-group structure.

    This post-processing step is used only for the visual embedding so that
    major reference clades and candidate groups no longer collapse onto the
    same centroid in the final SVG.
    """
    if len(coords) == 0:
        return coords.copy()

    dimensions = coords.shape[1]
    group_to_indices: dict[str, list[int]] = {}
    for index, label in enumerate(labels):
        group_to_indices.setdefault(label, []).append(index)
    if len(group_to_indices) < 2:
        return coords.copy()

    centroids = {
        label: coords[indices].mean(axis=0)
        for label, indices in group_to_indices.items()
    }
    radii = {}
    for label, indices in group_to_indices.items():
        residuals = coords[indices] - centroids[label]
        norms = np.linalg.norm(residuals, axis=1) if len(indices) else np.zeros((0,), dtype=float)
        radii[label] = float(np.percentile(norms, 75)) if len(norms) else 0.0
    median_radius = float(np.median(list(radii.values()))) if radii else 0.0
    target_gap = max(median_radius * 4.0, 1.5)

    shifted = {label: centroid.copy() for label, centroid in centroids.items()}
    labels_sorted = sorted(group_to_indices)
    for _ in range(48):
        moved = False
        for pair_index, label_a in enumerate(labels_sorted):
            for label_b in labels_sorted[pair_index + 1 :]:
                delta = shifted[label_a] - shifted[label_b]
                distance = float(np.linalg.norm(delta))
                min_distance = max(target_gap, radii[label_a] + radii[label_b] + (median_radius * 1.5))
                if distance >= min_distance:
                    continue
                if distance < 1e-9:
                    direction = _deterministic_direction(pair_index, dimensions)
                else:
                    direction = delta / distance
                push = (min_distance - distance) * 0.5
                shifted[label_a] = shifted[label_a] + (direction * push)
                shifted[label_b] = shifted[label_b] - (direction * push)
                moved = True
        for label in labels_sorted:
            shifted[label] = shifted[label] * 0.98 + centroids[label] * 0.02
        if not moved:
            break

    adjusted = np.zeros_like(coords)
    for label, indices in group_to_indices.items():
        residuals = coords[indices] - centroids[label]
        adjusted[indices] = shifted[label] + (residuals * within_scale)
    return adjusted


def _embedding_coordinates(
    matrix: np.ndarray,
    records: list[FastaRecord],
    *,
    n_components: int = 3,
) -> tuple[np.ndarray, np.ndarray, list[str], str]:
    """Compute embedding coordinates, preferring class-aware LDA over PCA.

    LDA is used to better separate reference clades and candidate groups in the
    visual embedding. When the label configuration is too small or unstable for
    LDA, the function falls back to PCA.
    """
    if matrix.size == 0:
        return (
            np.zeros((0, n_components), dtype=float),
            np.zeros((n_components,), dtype=float),
            [f"Axis{i + 1}" for i in range(n_components)],
            "empty",
        )

    scaled_matrix = _zscore_matrix(matrix)
    embedding_labels = _embedding_group_labels(scaled_matrix, records)
    unique_labels = sorted(set(embedding_labels))
    if len(unique_labels) >= 2:
        max_components = min(n_components, len(unique_labels) - 1, scaled_matrix.shape[1])
        if max_components >= 2:
            try:
                lda = LinearDiscriminantAnalysis(
                    solver="eigen",
                    shrinkage="auto",
                    n_components=max_components,
                )
                coords = lda.fit_transform(scaled_matrix, embedding_labels)
                explained = np.asarray(getattr(lda, "explained_variance_ratio_", np.zeros((coords.shape[1],), dtype=float)))
                if coords.shape[1] < n_components:
                    coords = np.pad(coords, ((0, 0), (0, n_components - coords.shape[1])))
                if explained.shape[0] < n_components:
                    explained = np.pad(explained, (0, n_components - explained.shape[0]))
                labels = [f"LD{i + 1}" for i in range(n_components)]
                spread_coords = _spread_group_coords(coords, embedding_labels)
                return spread_coords, explained[:n_components], labels, "lda_reference_subclade_spread"
            except Exception as error:
                logger.warning("LDA embedding failed, falling back to PCA: %s", error)

    coords, explained = _pca_coordinates(scaled_matrix, n_components=n_components)
    spread_coords = _spread_group_coords(coords, embedding_labels)
    return spread_coords, explained, [f"PC{i + 1}" for i in range(n_components)], "pca_reference_subclade_spread"


def _distance_matrix(features: np.ndarray) -> np.ndarray:
    """Compute Euclidean distances between all rows in the feature matrix."""
    if len(features) == 0:
        return np.zeros((0, 0), dtype=float)
    squared = np.sum(features**2, axis=1, keepdims=True)
    distances = squared + squared.T - (2 * features @ features.T)
    np.maximum(distances, 0, out=distances)
    return np.sqrt(distances)


def _upgma_newick(names: list[str], distances: np.ndarray) -> str:
    """Build a simple UPGMA tree from a pairwise distance matrix."""
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
    """Assign deterministic colours to reference sources for plotting."""
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


def _candidate_group(record: FastaRecord) -> str:
    """Map a record to a candidate/reference display category."""
    if record.metadata.get("source") != "candidate":
        return f"ref:{record.metadata.get('source', 'unknown')}"
    if record.metadata.get("is_ceess_candidate") == "yes":
        return "candidate_ceess"
    if record.metadata.get("esm_ceess_label") == "non-CeeSs" or record.metadata.get("is_coral_like") == "yes":
        return "candidate_non_ceess"
    return "candidate"


def _marker_svg(x: float, y: float, *, shape: str, radius: float, fill: str, stroke: str, opacity: float, stroke_width: float = 1.0) -> str:
    """Render one SVG marker in the requested geometric style."""
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


def _candidate_group_overlay(points: list[tuple[float, float]], *, fill: str, stroke: str, label: str) -> str:
    """Render a soft ellipse + label that highlights one candidate subgroup."""
    if not points:
        return ""
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    cx = (min_x + max_x) / 2
    cy = (min_y + max_y) / 2
    rx = max((max_x - min_x) / 2 + 18.0, 24.0)
    ry = max((max_y - min_y) / 2 + 16.0, 22.0)
    label_y = cy - ry - 10.0
    return (
        f'<ellipse cx="{cx:.2f}" cy="{cy:.2f}" rx="{rx:.2f}" ry="{ry:.2f}" '
        f'fill="{fill}" fill-opacity="0.08" stroke="{stroke}" stroke-opacity="0.6" '
        f'stroke-width="1.5" stroke-dasharray="6,4" />'
        f'<text x="{cx:.2f}" y="{label_y:.2f}" text-anchor="middle" font-size="11" '
        f'font-weight="700" fill="{stroke}">{label}</text>'
    )


def _render_scatter(
    records: list[FastaRecord],
    coords: np.ndarray,
    output_path: PathLike,
    *,
    explained: Optional[np.ndarray] = None,
    component_labels: Optional[list[str]] = None,
    method: str = "embedding",
) -> Path:
    """Render a readable 2D embedding scatter plot for candidate classification."""
    target = Path(output_path)
    if len(records) == 0:
        target.write_text("<svg xmlns='http://www.w3.org/2000/svg' width='1200' height='820'></svg>")
        return target
    if coords.shape[1] < 2:
        coords = np.pad(coords, ((0, 0), (0, 2 - coords.shape[1])))

    width = 1200
    height = 820
    outer_left = 78
    outer_top = 42
    outer_right = 46
    outer_bottom = 124
    plot_left = outer_left
    plot_top = outer_top + 64
    plot_width = width - plot_left - outer_right - 250
    plot_height = height - plot_top - outer_bottom
    plot_right = plot_left + plot_width
    plot_bottom = plot_top + plot_height
    font = "Arial, Helvetica, sans-serif"
    axis_color = "#374151"
    grid_color = "#E5E7EB"
    xs = coords[:, 0]
    ys = coords[:, 1]
    min_x, max_x = float(xs.min()), float(xs.max())
    min_y, max_y = float(ys.min()), float(ys.max())
    pad_x = (max_x - min_x) * 0.05 or 0.5
    pad_y = (max_y - min_y) * 0.05 or 0.5
    min_x -= pad_x
    max_x += pad_x
    min_y -= pad_y
    max_y += pad_y
    scale_x, scale_y = _make_scalers(min_x, max_x, min_y, max_y, plot_left, plot_top, plot_bottom, plot_width, plot_height)

    ref_sources: list[str] = []
    for record in records:
        if record.metadata.get("source") != "candidate":
            source = record.metadata.get("source", "unknown")
            if source not in ref_sources:
                ref_sources.append(source)
    ref_colour = {source: _ref_source_fill(source, index) for index, source in enumerate(ref_sources)}
    candidate_groups = [
        group for group in _CANDIDATE_CFG
        if any(_candidate_group(record) == group for record in records if record.metadata.get("source") == "candidate")
    ]

    axis_names = component_labels or [f"Axis{i + 1}" for i in range(max(2, coords.shape[1]))]
    candidate_overlay_points: dict[str, list[tuple[float, float]]] = {group: [] for group in candidate_groups}

    def axis_label(dim: int) -> str:
        name = axis_names[dim] if dim < len(axis_names) else f"Axis{dim + 1}"
        if explained is not None and len(explained) > dim:
            return f"{name} ({100.0 * float(explained[dim]):.1f}%)"
        return name

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" font-family="{font}">',
        f'<rect width="{width}" height="{height}" fill="white" />',
        f'<text x="{outer_left}" y="{outer_top + 24}" font-size="22" font-weight="700" fill="#111827">Ariadne TPS Feature Embedding</text>',
        f'<text x="{outer_left}" y="{outer_top + 46}" font-size="11" fill="#6B7280" font-style="italic">Reference clades define the background feature space; candidate markers show discovered sequences. Projection: {method.upper()}.</text>',
        f'<rect x="{plot_left}" y="{plot_top}" width="{plot_width}" height="{plot_height}" fill="white" stroke="{grid_color}" stroke-width="1"/>',
    ]

    for tick in _nice_ticks(min_x, max_x):
        tx = scale_x(tick)
        svg.append(f'<line x1="{tx:.2f}" y1="{plot_top}" x2="{tx:.2f}" y2="{plot_bottom}" stroke="{grid_color}" stroke-width="0.7" stroke-dasharray="3,3"/>')
        svg.append(f'<line x1="{tx:.2f}" y1="{plot_bottom}" x2="{tx:.2f}" y2="{plot_bottom + 5}" stroke="{axis_color}" stroke-width="1.2"/>')
        svg.append(f'<text x="{tx:.2f}" y="{plot_bottom + 18}" text-anchor="middle" font-size="10" fill="#6B7280">{_fmt_tick(tick)}</text>')
    for tick in _nice_ticks(min_y, max_y):
        ty = scale_y(tick)
        svg.append(f'<line x1="{plot_left}" y1="{ty:.2f}" x2="{plot_right}" y2="{ty:.2f}" stroke="{grid_color}" stroke-width="0.7" stroke-dasharray="3,3"/>')
        svg.append(f'<line x1="{plot_left - 5}" y1="{ty:.2f}" x2="{plot_left}" y2="{ty:.2f}" stroke="{axis_color}" stroke-width="1.2"/>')
        svg.append(f'<text x="{plot_left - 10}" y="{ty + 4:.2f}" text-anchor="end" font-size="10" fill="#6B7280">{_fmt_tick(tick)}</text>')

    svg.append(f'<line x1="{plot_left}" y1="{plot_bottom}" x2="{plot_right}" y2="{plot_bottom}" stroke="{axis_color}" stroke-width="1.5"/>')
    svg.append(f'<line x1="{plot_left}" y1="{plot_top}" x2="{plot_left}" y2="{plot_bottom}" stroke="{axis_color}" stroke-width="1.5"/>')
    svg.append(f'<text x="{plot_left + (plot_width / 2):.2f}" y="{plot_bottom + 40}" text-anchor="middle" font-size="13" font-weight="600" fill="{axis_color}">{axis_label(0)}</text>')
    y_label_x = outer_left - 54
    y_label_y = plot_top + (plot_height / 2)
    svg.append(
        f'<text transform="rotate(-90,{y_label_x:.2f},{y_label_y:.2f})" x="{y_label_x:.2f}" y="{y_label_y:.2f}" text-anchor="middle" font-size="13" font-weight="600" fill="{axis_color}">{axis_label(1)}</text>'
    )

    for record, x_value, y_value in zip(records, xs, ys):
        x = max(plot_left, min(plot_right, scale_x(float(x_value))))
        y = max(plot_top, min(plot_bottom, scale_y(float(y_value))))
        if record.metadata.get("source") != "candidate":
            source = record.metadata.get("source", "unknown")
            svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.6" fill="{ref_colour.get(source, "#9CA3AF")}" fill-opacity="0.58" stroke="none"/>')
            continue
        group = _candidate_group(record)
        cfg = _CANDIDATE_CFG.get(group, _CANDIDATE_CFG["candidate"])
        if group in candidate_overlay_points:
            candidate_overlay_points[group].append((x, y))
        svg.append(
            _marker_svg(
                x,
                y,
                shape=cfg["shape"],
                radius=cfg["r"],
                fill=cfg["fill"],
                stroke=cfg["stroke"],
                opacity=cfg["op"],
                stroke_width=1.3,
            )
        )

    for group in candidate_groups:
        cfg = _CANDIDATE_CFG[group]
        overlay = _candidate_group_overlay(
            candidate_overlay_points.get(group, []),
            fill=cfg["fill"],
            stroke=cfg["stroke"],
            label=cfg["label"],
        )
        if overlay:
            svg.insert(5, overlay)

    legend_x = plot_right + 36
    legend_y = plot_top + 24
    svg.append(f'<text x="{legend_x}" y="{legend_y}" font-size="12" font-weight="600" fill="#374151">Reference clades</text>')
    current_y = legend_y + 22
    for source in ref_sources:
        svg.append(f'<circle cx="{legend_x + 6}" cy="{current_y - 4}" r="5.5" fill="{ref_colour[source]}" fill-opacity="0.7" stroke="none"/>')
        svg.append(f'<text x="{legend_x + 18}" y="{current_y}" font-size="12" fill="#374151">{source}</text>')
        current_y += 24
    current_y += 10
    svg.append(f'<text x="{legend_x}" y="{current_y}" font-size="12" font-weight="600" fill="#374151">Candidate groups</text>')
    current_y += 22
    for group in candidate_groups:
        cfg = _CANDIDATE_CFG[group]
        svg.append(_marker_svg(legend_x + 7, current_y - 5, shape=cfg["shape"], radius=7, fill=cfg["fill"], stroke=cfg["stroke"], opacity=cfg["op"], stroke_width=1.0))
        svg.append(f'<text x="{legend_x + 20}" y="{current_y}" font-size="12" fill="#374151">{cfg["label"]}</text>')
        current_y += 24
    svg.append(f'<text x="{outer_left}" y="{height - 24}" font-size="10" fill="#9CA3AF" font-style="italic">Candidate IDs are omitted from embedding.svg on purpose so the SVG stays readable. Use classification.tsv and embedding.tsv for per-sequence details.</text>')
    svg.append("</svg>")
    target.write_text("".join(svg))
    return target


# ---------------------------------------------------------------------------
# Publication-quality palette — AFPK_finder colour scheme
# (script.r: scale_color_manual + theme_bw + no grid lines)
# ---------------------------------------------------------------------------

# Named source → fill colour (mirrors AFPK_finder 12-colour palette)
_SOURCE_FILL: dict[str, str] = {
    "coral":    "#EE6A50",   # coral-red
    "insect":   "#7B68EE",   # medium-slate-blue
    "bacteria": "#9ACD32",   # yellow-green
    "fungal":   "#49C3C3",   # teal
    "plant":    "#3A7D44",   # dark-green
}
# Fallback colours for unlisted sources (remaining AFPK palette entries)
_FILL_EXTRA = ["#87CEFA", "#FFC0CB", "#800080", "#191970", "#FF7F50", "#808080", "#FFD700"]

# Candidate-group display configuration (shape + colour, NO text labels)
_CANDIDATE_CFG: dict[str, dict] = {
    "candidate_ceess": {
        "label": "Candidate CeeSs",
        "fill": "#D97706", "stroke": "#7C2D12", "shape": "diamond", "r": 8.2, "op": 0.94,
    },
    "candidate_non_ceess": {
        "label": "Candidate non-CeeSs",
        "fill": "#0F766E", "stroke": "#134E4A", "shape": "square", "r": 7.6, "op": 0.9,
    },
    "candidate": {
        "label": "Candidate other",
        "fill": "#374151", "stroke": "#111827", "shape": "circle",   "r": 7.2, "op": 0.84,
    },
}


def _ref_source_fill(source: str, index: int) -> str:
    """Return a fallback fill colour for one reference source."""
    return _SOURCE_FILL.get(source, _FILL_EXTRA[index % len(_FILL_EXTRA)])


def _nice_ticks(lo: float, hi: float, target_n: int = 5) -> list[float]:
    """Return ≈ *target_n* human-readable tick values spanning [lo, hi]."""
    span = hi - lo
    if span < 1e-12:
        return [round(lo, 4)]
    raw = span / target_n
    mag = 10 ** math.floor(math.log10(abs(raw)))
    step = mag
    for c in (1, 2, 2.5, 5, 10):
        if c * mag >= raw:
            step = c * mag
            break
    start = math.ceil(lo / step - 1e-9) * step
    ticks: list[float] = []
    v = start
    while v <= hi + step * 1e-6:
        ticks.append(round(v, 10))
        v += step
    return ticks


def _fmt_tick(v: float) -> str:
    """Format an axis tick label with compact precision."""
    if abs(v) < 1e-12:
        return "0"
    if v == int(v):
        return str(int(v))
    return f"{v:.2g}"


def _make_scalers(
    dxlo: float, dxhi: float,
    dylo: float, dyhi: float,
    plot_x0: float, plot_y0: float, plot_y1: float,
    plot_w: float, plot_h: float,
):
    """Return (svgx, svgy) callables that map data coordinates to SVG pixels."""
    span_x = max(dxhi - dxlo, 1e-9)
    span_y = max(dyhi - dylo, 1e-9)

    def svgx(v: float) -> float:
        return plot_x0 + (v - dxlo) / span_x * plot_w

    def svgy(v: float) -> float:
        return plot_y1 - (v - dylo) / span_y * plot_h

    return svgx, svgy


def _render_3d_sections(
    records: list[FastaRecord],
    coords: np.ndarray,
    output_path: PathLike,
    *,
    explained: Optional[np.ndarray] = None,
    component_labels: Optional[list[str]] = None,
    method: str = "embedding",
) -> Path:
    """Render a publication-quality embedding-section figure (3 orthogonal projections).

    Visual style mirrors AFPK_finder (Lin *et al.*, Nat. Commun. 2024):

    * Reference sequences coloured by clade (``theme_bw`` analogue — white
      background, no grid, light axis lines).
    * Candidate markers rendered as distinct shapes **without text labels**,
      eliminating the visual clutter reported when sequence IDs overlap.
    * Axis labels include the fraction of separation/variance explained by each axis.
    * A two-row legend (reference clades / candidate groups) is placed below
      the panels.
    """
    target = Path(output_path)
    if len(records) == 0:
        target.write_text(
            "<svg xmlns='http://www.w3.org/2000/svg' width='1580' height='680'></svg>"
        )
        return target
    if coords.shape[1] < 3:
        coords = np.pad(coords, ((0, 0), (0, 3 - coords.shape[1])))

    # ── reference sources in appearance order ─────────────────────────────
    ref_sources: list[str] = []
    for rec in records:
        if rec.metadata.get("source") != "candidate":
            src = rec.metadata.get("source", "unknown")
            if src not in ref_sources:
                ref_sources.append(src)
    ref_colour = {src: _ref_source_fill(src, i) for i, src in enumerate(ref_sources)}

    # candidate groups actually present in this dataset
    cand_groups_present = [
        g for g in _CANDIDATE_CFG
        if any(
            _candidate_group(r) == g
            for r in records
            if r.metadata.get("source") == "candidate"
        )
    ]

    # ── layout constants ──────────────────────────────────────────────────
    PW, PH   = 460, 420          # panel width, height
    PGAP     = 50                # gap between panels
    OL, OT   = 45, 40            # outer left / top margin
    OB       = 130               # outer bottom margin (legend area)
    PAD_L    = 64                # within-panel left  (Y-axis labels)
    PAD_R    = 18                # within-panel right (breathing room)
    PAD_T    = 36                # within-panel top   (sub-label)
    PAD_B    = 54                # within-panel bottom (X-axis label)
    PLOT_W   = PW - PAD_L - PAD_R
    PLOT_H   = PH - PAD_T - PAD_B
    TICK_LEN = 5
    FONT     = "Arial, Helvetica, sans-serif"
    AX_COL   = "#374151"
    GRID_COL = "#E5E7EB"

    TOTAL_W  = OL + 3 * PW + 2 * PGAP + OL
    TOTAL_H  = OT + 58 + PH + OB      # 58 = title block height
    PANELS_Y = OT + 58                 # top edge of the panel row

    axis_names = component_labels or [f"Axis{i + 1}" for i in range(max(3, coords.shape[1]))]
    projections = [
        (0, 1),
        (0, 2),
        (1, 2),
    ]

    def axis_label(dim: int) -> str:
        name = axis_names[dim] if dim < len(axis_names) else f"Axis{dim + 1}"
        if explained is not None and len(explained) > dim:
            return f"{name} ({100.0 * float(explained[dim]):.1f}%)"
        return name

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg"'
        f' width="{TOTAL_W}" height="{TOTAL_H}"'
        f' viewBox="0 0 {TOTAL_W} {TOTAL_H}"'
        f' font-family="{FONT}">',
        f'<rect width="{TOTAL_W}" height="{TOTAL_H}" fill="white"/>',
        # figure title
        f'<text x="{OL}" y="{OT + 26}"'
        f' font-size="19" font-weight="700" fill="#111827">'
        f'Ariadne TPS Feature Space — {method.upper()} Embedding</text>',
        f'<text x="{OL}" y="{OT + 48}"'
        f' font-size="11" fill="#6B7280" font-style="italic">'
        f'Background: reference clades.  Foreground: discovered candidates. '
        f'Sequence labels omitted for visual clarity — see classification.tsv.</text>',
    ]

    for pi, (dx, dy) in enumerate(projections):
        px0 = OL + pi * (PW + PGAP)     # panel left edge (SVG x)
        py0 = PANELS_Y                   # panel top  edge (SVG y)

        # plot-area corners
        plot_x0 = px0 + PAD_L
        plot_y0  = py0 + PAD_T
        plot_x1  = plot_x0 + PLOT_W
        plot_y1  = plot_y0 + PLOT_H
        mid_x    = (plot_x0 + plot_x1) / 2
        mid_y    = (plot_y0 + plot_y1) / 2

        xs = coords[:, dx]
        ys = coords[:, dy]
        dxlo, dxhi = float(xs.min()), float(xs.max())
        dylo, dyhi = float(ys.min()), float(ys.max())

        # 5 % padding around data extent
        xpad = (dxhi - dxlo) * 0.05 or 0.5
        ypad = (dyhi - dylo) * 0.05 or 0.5
        dxlo -= xpad;  dxhi += xpad
        dylo -= ypad;  dyhi += ypad

        svgx, svgy = _make_scalers(dxlo, dxhi, dylo, dyhi,
                                   plot_x0, plot_y0, plot_y1, PLOT_W, PLOT_H)

        # panel background (very light grey — theme_bw analogue)
        svg.append(
            f'<rect x="{px0}" y="{py0}" width="{PW}" height="{PH}"'
            f' fill="#F9FAFB" rx="4"/>',
        )
        # inner plot area (white, thin border)
        svg.append(
            f'<rect x="{plot_x0}" y="{plot_y0}" width="{PLOT_W}" height="{PLOT_H}"'
            f' fill="white" stroke="{GRID_COL}" stroke-width="1"/>',
        )

        # faint dashed grid lines (panel.grid analogue — very subtle)
        for tv in _nice_ticks(dxlo, dxhi):
            if dxlo <= tv <= dxhi:
                tx = svgx(tv)
                svg.append(
                    f'<line x1="{tx:.2f}" y1="{plot_y0}" x2="{tx:.2f}" y2="{plot_y1}"'
                    f' stroke="{GRID_COL}" stroke-width="0.7" stroke-dasharray="3,3"/>',
                )
        for tv in _nice_ticks(dylo, dyhi):
            if dylo <= tv <= dyhi:
                ty = svgy(tv)
                svg.append(
                    f'<line x1="{plot_x0}" y1="{ty:.2f}" x2="{plot_x1}" y2="{ty:.2f}"'
                    f' stroke="{GRID_COL}" stroke-width="0.7" stroke-dasharray="3,3"/>',
                )

        # X-axis ticks + numeric labels
        for tv in _nice_ticks(dxlo, dxhi):
            if dxlo <= tv <= dxhi:
                tx = svgx(tv)
                svg.append(
                    f'<line x1="{tx:.2f}" y1="{plot_y1}"'
                    f' x2="{tx:.2f}" y2="{plot_y1 + TICK_LEN}"'
                    f' stroke="{AX_COL}" stroke-width="1.2"/>',
                )
                svg.append(
                    f'<text x="{tx:.2f}" y="{plot_y1 + TICK_LEN + 13}"'
                    f' text-anchor="middle" font-size="10" fill="#6B7280">'
                    f'{_fmt_tick(tv)}</text>',
                )

        # Y-axis ticks + numeric labels
        for tv in _nice_ticks(dylo, dyhi):
            if dylo <= tv <= dyhi:
                ty = svgy(tv)
                svg.append(
                    f'<line x1="{plot_x0 - TICK_LEN}" y1="{ty:.2f}"'
                    f' x2="{plot_x0}" y2="{ty:.2f}"'
                    f' stroke="{AX_COL}" stroke-width="1.2"/>',
                )
                svg.append(
                    f'<text x="{plot_x0 - TICK_LEN - 4}" y="{ty + 4:.2f}"'
                    f' text-anchor="end" font-size="10" fill="#6B7280">'
                    f'{_fmt_tick(tv)}</text>',
                )

        # solid axis lines (X bottom, Y left)
        svg.append(
            f'<line x1="{plot_x0}" y1="{plot_y1}" x2="{plot_x1}" y2="{plot_y1}"'
            f' stroke="{AX_COL}" stroke-width="1.5"/>',
        )
        svg.append(
            f'<line x1="{plot_x0}" y1="{plot_y0}" x2="{plot_x0}" y2="{plot_y1}"'
            f' stroke="{AX_COL}" stroke-width="1.5"/>',
        )

        # axis labels with variance explained
        svg.append(
            f'<text x="{mid_x:.2f}" y="{plot_y1 + TICK_LEN + 30}"'
            f' text-anchor="middle" font-size="13" font-weight="600" fill="{AX_COL}">'
            f'{axis_label(dx)}</text>',
        )
        yl_x = px0 + 13
        svg.append(
            f'<text transform="rotate(-90,{yl_x:.2f},{mid_y:.2f})"'
            f' x="{yl_x:.2f}" y="{mid_y:.2f}"'
            f' text-anchor="middle" font-size="13" font-weight="600" fill="{AX_COL}">'
            f'{axis_label(dy)}</text>',
        )

        # panel sub-label: a / b / c  (top-left, publication convention)
        svg.append(
            f'<text x="{px0 + 8}" y="{py0 + 24}"'
            f' font-size="16" font-weight="700" fill="#111827">'
            f'{chr(ord("a") + pi)}</text>',
        )

        # ── data points ─────────────────────────────────────────────────────
        # paint references (background) before candidates (foreground)
        ref_pts: list[str] = []
        cand_pts: list[str] = []
        cand_overlay_points: dict[str, list[tuple[float, float]]] = {group: [] for group in cand_groups_present}

        for rec, xv, yv in zip(records, xs, ys):
            cx = max(plot_x0, min(plot_x1, svgx(float(xv))))
            cy = max(plot_y0, min(plot_y1, svgy(float(yv))))
            src = rec.metadata.get("source", "unknown")

            if src != "candidate":
                # reference point: small filled circle, coloured by clade
                col = ref_colour.get(src, "#9CA3AF")
                ref_pts.append(
                    f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="3.5"'
                    f' fill="{col}" fill-opacity="0.55" stroke="none"/>',
                )
            else:
                # candidate: larger shape marker, NO text label
                grp = _candidate_group(rec)
                cfg = _CANDIDATE_CFG.get(grp, _CANDIDATE_CFG["candidate"])
                if grp in cand_overlay_points:
                    cand_overlay_points[grp].append((cx, cy))
                cand_pts.append(
                    _marker_svg(
                        cx, cy,
                        shape=cfg["shape"],
                        radius=cfg["r"],
                        fill=cfg["fill"],
                        stroke=cfg["stroke"],
                        opacity=cfg["op"],
                        stroke_width=1.3,
                    )
                )

        for grp in cand_groups_present:
            cfg = _CANDIDATE_CFG[grp]
            overlay = _candidate_group_overlay(
                cand_overlay_points.get(grp, []),
                fill=cfg["fill"],
                stroke=cfg["stroke"],
                label=cfg["label"],
            )
            if overlay:
                svg.append(overlay)
        svg.extend(ref_pts)
        svg.extend(cand_pts)

    # ── two-row legend ────────────────────────────────────────────────────
    leg_y = PANELS_Y + PH + 24

    # row 1 — reference clades
    svg.append(
        f'<text x="{OL}" y="{leg_y}"'
        f' font-size="12" font-weight="600" fill="#374151">'
        f'Reference clades:</text>',
    )
    ix = OL + 128
    for src in ref_sources:
        col = ref_colour[src]
        svg.append(
            f'<circle cx="{ix + 6}" cy="{leg_y - 5}" r="5.5"'
            f' fill="{col}" fill-opacity="0.7" stroke="none"/>',
        )
        svg.append(
            f'<text x="{ix + 16}" y="{leg_y}"'
            f' font-size="12" fill="#374151">{src}</text>',
        )
        ix += max(66, 9 * len(src) + 28)

    # row 2 — candidate groups
    leg_y2 = leg_y + 30
    svg.append(
        f'<text x="{OL}" y="{leg_y2}"'
        f' font-size="12" font-weight="600" fill="#374151">'
        f'Candidate groups:</text>',
    )
    ix2 = OL + 128
    for grp in cand_groups_present:
        cfg = _CANDIDATE_CFG[grp]
        svg.append(
            _marker_svg(ix2 + 7, leg_y2 - 5,
                        shape=cfg["shape"], radius=7,
                        fill=cfg["fill"], stroke=cfg["stroke"],
                        opacity=cfg["op"], stroke_width=1.0)
        )
        svg.append(
            f'<text x="{ix2 + 19}" y="{leg_y2}"'
            f' font-size="12" fill="#374151">{cfg["label"]}</text>',
        )
        ix2 += max(82, 9 * len(cfg["label"]) + 32)

    # footnote
    fn_y = leg_y2 + 30
    svg.append(
        f'<text x="{OL}" y="{fn_y}"'
        f' font-size="10" fill="#9CA3AF" font-style="italic">'
        f'{method.upper()} embedding computed from HMM-score feature matrix (column-mean normalised, z-scored for display). '
        f'Candidate IDs omitted from plot; see classification.tsv for per-sequence assignments.</text>',
    )

    svg.append("</svg>")
    target.write_text("".join(svg))
    return target


def classify_candidates(
    candidate_fasta: PathLike,
    reference_dir: PathLike,
    output_dir: PathLike,
    *,
    hmm_dir: PathLike,
    top_k: int = 5,
    tree_neighbors: int = 12,
    ceess_xlsx: Optional[PathLike] = None,
    ceess_model_name: str = DEFAULT_ESM_MODEL_NAME,
    ceess_batch_size: int = 4,
    ceess_max_length: int = 2048,
    ceess_device: Optional[str] = None,
    ceess_cv_folds: int = 5,
    ceess_random_state: int = 0,
    ceess_threshold: float = 0.5,
) -> dict[str, Path]:
    """Assign each candidate to the dominant source of its nearest references."""
    references = load_reference_records(reference_dir)
    if not references:
        raise ValueError(f"No reference sequences were found in {reference_dir}.")
    candidates = read_fasta(candidate_fasta)
    for candidate in candidates:
        candidate.metadata["source"] = "candidate"
        candidate.metadata["label"] = "candidate"
        candidate.metadata["candidate_group"] = _candidate_group(candidate)
    all_records = references + candidates
    hmm_paths = sorted_hmm_paths(hmm_dir)
    raw_matrix, column_names = _score_records_against_hmms(all_records, hmm_paths)
    normalized_matrix = _normalize_matrix(raw_matrix)
    coords, explained, component_labels, embedding_method = _embedding_coordinates(
        normalized_matrix,
        all_records,
        n_components=3,
    )
    embedding_groups = _embedding_group_labels(_zscore_matrix(normalized_matrix), all_records)
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

    reference_count = len(references)
    candidate_rows: list[dict[str, object]] = []
    neighbor_rows: list[dict[str, object]] = []
    cluster_context_rows: list[dict[str, object]] = []
    candidate_records_by_id: dict[str, FastaRecord] = {}
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
        is_coral_like = predicted_source == "coral" or nearest_reference.metadata.get("source", "unknown") == "coral"
        candidate.metadata["is_coral_like"] = "yes" if is_coral_like else "no"
        candidate.metadata["is_ceess_candidate"] = "no"
        candidate.metadata["esm_ceess_label"] = ""
        candidate_records_by_id[candidate.id] = candidate
        candidate_rows.append(
            {
                "sequence_id": candidate.id,
                "candidate_group": candidate.metadata.get("candidate_group", "candidate"),
                "predicted_source": predicted_source,
                "confidence": round(confidence, 4),
                "is_coral_like": "yes" if is_coral_like else "no",
                "nearest_reference": nearest_reference.id,
                "nearest_reference_source": nearest_reference.metadata.get("source", "unknown"),
                "nearest_distance": round(nearest_distance, 6),
                "esm_type_prediction": "",
                "esm_ceess_label": "",
                "esm_ceess_subtype": "",
                "esm_ceess_probability": "",
                "esm_non_ceess_probability": "",
                "esm_cembrene_a_probability": "",
                "esm_cembrene_b_probability": "",
                "is_ceess_candidate": "no",
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
                "candidate_group": candidate.metadata.get("candidate_group", "candidate"),
                "predicted_source": predicted_source,
                "nearest_reference_source": nearest_reference.metadata.get("source", "unknown"),
                "nearest_distance": round(nearest_distance, 6),
                "embedding_group": embedding_groups[candidate_index],
                "embedding_method": embedding_method,
                "axis1": round(float(coords[candidate_index, 0]), 6),
                "axis2": round(float(coords[candidate_index, 1]), 6),
                "axis3": round(float(coords[candidate_index, 2]), 6),
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

    ceess_outputs: dict[str, Path] = {}
    if ceess_xlsx is not None:
        ceess_xlsx_path = Path(ceess_xlsx)
        coral_like_ids = [str(row["sequence_id"]) for row in candidate_rows if row["is_coral_like"] == "yes"]
        coral_like_records = [candidate_records_by_id[sequence_id] for sequence_id in coral_like_ids]
        if coral_like_records and ceess_xlsx_path.exists():
            try:
                from ariadne.esm_type import classify_ceess_candidates_with_esm

                ceess_result = classify_ceess_candidates_with_esm(
                    coral_like_records,
                    ceess_xlsx_path,
                    output_root,
                    model_name=ceess_model_name,
                    batch_size=ceess_batch_size,
                    max_length=ceess_max_length,
                    device=ceess_device,
                    cv_folds=ceess_cv_folds,
                    random_state=ceess_random_state,
                    ceess_threshold=ceess_threshold,
                )
                by_id = {str(row["sequence_id"]): row for row in ceess_result.prediction_rows}
                for row in candidate_rows:
                    prediction = by_id.get(str(row["sequence_id"]))
                    if prediction is None:
                        continue
                    candidate_record = candidate_records_by_id.get(str(row["sequence_id"]))
                    if candidate_record is not None:
                        candidate_record.metadata["is_ceess_candidate"] = prediction["is_ceess_candidate"]
                        candidate_record.metadata["esm_ceess_label"] = prediction.get("esm_ceess_label", "")
                    row["esm_type_prediction"] = prediction["esm_type_prediction"]
                    row["esm_ceess_label"] = prediction.get("esm_ceess_label", "")
                    row["esm_ceess_subtype"] = prediction["esm_ceess_subtype"]
                    row["esm_ceess_probability"] = prediction["esm_ceess_probability"]
                    row["esm_non_ceess_probability"] = prediction.get("esm_non_ceess_probability", "")
                    row["esm_cembrene_a_probability"] = prediction["esm_cembrene_a_probability"]
                    row["esm_cembrene_b_probability"] = prediction["esm_cembrene_b_probability"]
                    row["is_ceess_candidate"] = prediction["is_ceess_candidate"]
                classification_path = write_tsv(candidate_rows, output_root / "classification.tsv")
                ceess_outputs = ceess_result.output_paths
            except Exception as error:
                logger.warning("Skipped optional CeeSs ESM stage: %s", error)
        elif not ceess_xlsx_path.exists():
            logger.warning("CeeSs ESM spreadsheet not found, skipped: %s", ceess_xlsx_path)

    coords, explained, component_labels, embedding_method = _embedding_coordinates(
        normalized_matrix,
        all_records,
        n_components=3,
    )
    embedding_groups = _embedding_group_labels(_zscore_matrix(normalized_matrix), all_records)

    embedding_rows = []
    for index, record in enumerate(all_records):
        embedding_rows.append(
            {
                "sequence_id": record.id,
                "source": record.metadata.get("source", "unknown"),
                "candidate_group": _candidate_group(record) if record.metadata.get("source") == "candidate" else record.metadata.get("candidate_group", ""),
                "embedding_group": embedding_groups[index],
                "embedding_method": embedding_method,
                "axis1_label": component_labels[0] if len(component_labels) > 0 else "Axis1",
                "axis2_label": component_labels[1] if len(component_labels) > 1 else "Axis2",
                "axis3_label": component_labels[2] if len(component_labels) > 2 else "Axis3",
                "axis1": round(float(coords[index, 0]), 6),
                "axis2": round(float(coords[index, 1]), 6),
                "axis3": round(float(coords[index, 2]), 6),
                "pc1": round(float(coords[index, 0]), 6),
                "pc2": round(float(coords[index, 1]), 6),
                "pc3": round(float(coords[index, 2]), 6),
            }
        )
    embedding_path = write_tsv(embedding_rows, output_root / "embedding.tsv")
    variance_path = write_tsv(
        [
            {
                "embedding_method": embedding_method,
                "component": component_labels[component_index] if component_index < len(component_labels) else f"Axis{component_index + 1}",
                "explained_variance": round(float(value), 6),
            }
            for component_index, value in enumerate(explained, start=0)
        ],
        output_root / "embedding_variance.tsv",
    )
    cluster_context_rows = []
    for candidate_index in range(reference_count, len(all_records)):
        candidate = all_records[candidate_index]
        nearest_distance = float(
            next(row["nearest_distance"] for row in candidate_rows if row["sequence_id"] == candidate.id)
        )
        nearest_source = str(
            next(row["nearest_reference_source"] for row in candidate_rows if row["sequence_id"] == candidate.id)
        )
        predicted_source = str(
            next(row["predicted_source"] for row in candidate_rows if row["sequence_id"] == candidate.id)
        )
        cluster_context_rows.append(
            {
                "sequence_id": candidate.id,
                "candidate_group": _candidate_group(candidate),
                "predicted_source": predicted_source,
                "nearest_reference_source": nearest_source,
                "nearest_distance": round(nearest_distance, 6),
                "embedding_group": embedding_groups[candidate_index],
                "embedding_method": embedding_method,
                "axis1": round(float(coords[candidate_index, 0]), 6),
                "axis2": round(float(coords[candidate_index, 1]), 6),
                "axis3": round(float(coords[candidate_index, 2]), 6),
                "pc1": round(float(coords[candidate_index, 0]), 6),
                "pc2": round(float(coords[candidate_index, 1]), 6),
                "pc3": round(float(coords[candidate_index, 2]), 6),
            }
        )
    cluster_context_path = write_tsv(cluster_context_rows, output_root / "candidate_cluster_context.tsv")

    scatter_path = _render_scatter(
        all_records,
        coords,
        output_root / "embedding.svg",
        explained=explained,
        component_labels=component_labels,
        method=embedding_method,
    )
    sections_path = _render_3d_sections(
        all_records,
        coords,
        output_root / "embedding_3d_sections.svg",
        explained=explained,
        component_labels=component_labels,
        method=embedding_method,
    )
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
    outputs = {
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
    outputs.update(ceess_outputs)
    return outputs
