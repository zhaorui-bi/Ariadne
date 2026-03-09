from __future__ import annotations

import math
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pyhmmer

from ariadne.fasta_utils import FastaRecord, ensure_directory, read_fasta, sanitize_newick_name, write_fasta, write_tsv
from ariadne.references import load_reference_records


def sorted_hmm_paths(hmm_dir: str | Path) -> list[Path]:
    directory = Path(hmm_dir)
    hmm_paths = list(directory.glob("*.hmm"))
    if not hmm_paths:
        raise FileNotFoundError(f"No HMM profiles were found in {directory}.")
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


def _render_scatter(records: list[FastaRecord], coords: np.ndarray, output_path: str | Path) -> Path:
    if len(records) == 0:
        return Path(output_path)
    sources = [record.metadata.get("source", "unknown") for record in records]
    colors = _color_map(sources)
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
    for legend_index, source in enumerate(sorted(colors)):
        y = margin + (legend_index * 24)
        legend_entries.append(
            f'<circle cx="{width - 220}" cy="{y}" r="7" fill="{colors[source]}" />'
            f'<text x="{width - 205}" y="{y + 5}" font-size="14">{source}</text>'
        )

    points = []
    for record, coordinate in zip(records, coords):
        x = scale_x(float(coordinate[0]))
        y = scale_y(float(coordinate[1] if len(coordinate) > 1 else 0.0))
        source = record.metadata.get("source", "unknown")
        is_candidate = source == "candidate"
        radius = 8 if is_candidate else 4
        opacity = 0.95 if is_candidate else 0.55
        stroke = "#111111" if is_candidate else "none"
        points.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{radius}" fill="{colors[source]}" '
            f'opacity="{opacity}" stroke="{stroke}" stroke-width="1.2" />'
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
        f'<text x="{margin}" y="40" font-size="24" font-weight="bold">Ariadne TPS candidate feature embedding</text>'
        f'<line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" stroke="#999999" />'
        f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height - margin}" stroke="#999999" />'
        + "".join(points)
        + "".join(legend_entries)
        + "</svg>"
    )
    target = Path(output_path)
    target.write_text(svg)
    return target


def classify_candidates(
    candidate_fasta: str | Path,
    reference_dir: str | Path,
    output_dir: str | Path,
    *,
    hmm_dir: str | Path,
    top_k: int = 5,
    tree_neighbors: int = 12,
) -> dict[str, Path]:
    references = load_reference_records(reference_dir)
    if not references:
        raise ValueError(f"No reference sequences were found in {reference_dir}.")
    candidates = read_fasta(candidate_fasta)
    for candidate in candidates:
        candidate.metadata["source"] = "candidate"
        candidate.metadata["label"] = "candidate"
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
        }
        for column_name, value in zip(column_names, raw_matrix[index], strict=True):
            row[column_name] = round(float(value), 6)
        feature_rows.append(row)
    features_path = write_tsv(feature_rows, output_root / "aflp_features.tsv")

    embedding_rows = []
    for index, record in enumerate(all_records):
        embedding_rows.append(
            {
                "sequence_id": record.id,
                "source": record.metadata.get("source", "unknown"),
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

        local_neighbor_indices = [neighbor_index for neighbor_index, _ in reference_distances[: max(2, tree_neighbors)]]
        subset_indices = local_neighbor_indices + [candidate_index]
        subset_names = [all_records[index].id for index in subset_indices]
        subset_distances = distances[np.ix_(subset_indices, subset_indices)]
        tree_path = tree_dir / f"{candidate.id}.nwk"
        tree_path.write_text(_upgma_newick(subset_names, subset_distances) + "\n")

    classification_path = write_tsv(candidate_rows, output_root / "classification.tsv")
    neighbors_path = write_tsv(neighbor_rows, output_root / "nearest_neighbors.tsv")
    scatter_path = _render_scatter(all_records, coords, output_root / "embedding.svg")
    return {
        "features": features_path,
        "embedding": embedding_path,
        "variance": variance_path,
        "classification": classification_path,
        "neighbors": neighbors_path,
        "embedding_svg": scatter_path,
        "tree_dir": tree_dir,
    }
