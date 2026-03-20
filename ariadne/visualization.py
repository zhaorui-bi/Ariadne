from __future__ import annotations

import math
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE

from ariadne.fasta_utils import ensure_directory, write_tsv

PathLike = Union[str, Path]

LEGACY_PALETTE = [
    "#87CEFA",
    "#7B68EE",
    "#808080",
    "#9ACD32",
    "#99FF33",
    "#FFC0CB",
    "#EE6A50",
    "#800080",
    "#8DEEEE",
    "#006400",
    "#FFFF00",
    "#191970",
]


def _load_tps_features_tsv(path: PathLike) -> tuple[list[str], np.ndarray, list[str], list[str]]:
    import csv

    ids: list[str] = []
    clades: list[str] = []
    lengths: list[str] = []
    matrix_rows: list[list[float]] = []
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path} does not contain header fields.")
        non_feature_columns = {"sequence_id", "source", "label", "candidate_group"}
        feature_columns = [name for name in reader.fieldnames if name not in non_feature_columns]
        if not feature_columns:
            raise ValueError(f"{path} does not contain numeric feature columns.")
        for row in reader:
            sequence_id = (row.get("sequence_id") or "").strip()
            if not sequence_id:
                continue
            try:
                feature_row = [float((row.get(name) or "0").strip() or 0) for name in feature_columns]
            except ValueError:
                continue
            ids.append(sequence_id)
            clades.append((row.get("source") or "unknown").strip() or "unknown")
            lengths.append((row.get("length") or "").strip())
            matrix_rows.append(feature_row)
    if not matrix_rows:
        raise ValueError(f"No valid rows were parsed from {path}.")
    return ids, np.array(matrix_rows, dtype=float), clades, lengths


def _load_legacy_hmm_tab(path: PathLike) -> tuple[list[str], np.ndarray, list[str], list[str]]:
    ids: list[str] = []
    clades: list[str] = []
    lengths: list[str] = []
    matrix_rows: list[list[float]] = []
    with Path(path).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                feature_row = [float(value) for value in parts[1:-2]]
            except ValueError:
                # Skip header rows if present.
                continue
            ids.append(parts[0])
            clades.append(parts[-2])
            lengths.append(parts[-1])
            matrix_rows.append(feature_row)
    if not matrix_rows:
        raise ValueError(f"No valid rows were parsed from {path}.")
    return ids, np.array(matrix_rows, dtype=float), clades, lengths


def _load_hmm_score_table(path: PathLike) -> tuple[list[str], np.ndarray, list[str], list[str]]:
    first_line = ""
    with Path(path).open() as handle:
        for raw_line in handle:
            first_line = raw_line.strip()
            if first_line:
                break
    lowered = first_line.lower()
    if "sequence_id" in lowered and ("\t" in first_line):
        return _load_tps_features_tsv(path)
    return _load_legacy_hmm_tab(path)


def _normalize_by_column_mean(matrix: np.ndarray) -> np.ndarray:
    means = matrix.mean(axis=0)
    means[means == 0] = 1.0
    return matrix / means


def _default_perplexities(n_rows: int) -> list[int]:
    if n_rows <= 2:
        return [1]
    scale = round(n_rows / 100) / 10.0
    if scale <= 0:
        return [max(2, min(10, n_rows - 1))]
    raw = [math.ceil(40 * (scale**0.3)), math.ceil(60 * (scale**0.3))]
    return _sanitize_perplexities(raw, n_rows)


def _sanitize_perplexities(values: Sequence[int], n_rows: int) -> list[int]:
    if n_rows <= 2:
        return [1]
    normalized: list[int] = []
    for value in values:
        bounded = max(2, min(int(value), n_rows - 1))
        if bounded not in normalized:
            normalized.append(bounded)
    if normalized:
        return normalized
    return [max(2, min(30, n_rows - 1))]


def _pca_fallback(matrix: np.ndarray) -> np.ndarray:
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    if centered.shape[0] == 1:
        return np.zeros((1, 2), dtype=float)
    u, s, _ = np.linalg.svd(centered, full_matrices=False)
    dims = min(2, len(s))
    coords = u[:, :dims] * s[:dims]
    if dims < 2:
        coords = np.pad(coords, ((0, 0), (0, 2 - dims)))
    return coords


def _cluster_coordinates(coords: np.ndarray, min_points: int) -> np.ndarray:
    if len(coords) < min_points or len(coords) < 3:
        return np.full((len(coords),), -1, dtype=int)
    spread = float(np.std(coords))
    eps = max(0.5, spread * 0.35)
    model = DBSCAN(eps=eps, min_samples=min_points)
    return model.fit_predict(coords)


def _color_map(clades: list[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for index, clade in enumerate(sorted(set(clades))):
        mapping[clade] = LEGACY_PALETTE[index % len(LEGACY_PALETTE)]
    return mapping


def _render_tsne_svg(
    ids: list[str],
    coords: np.ndarray,
    clades: list[str],
    output_path: PathLike,
    *,
    title: str,
) -> Path:
    if len(ids) == 0:
        target = Path(output_path)
        target.write_text("<svg xmlns='http://www.w3.org/2000/svg' width='800' height='400'></svg>")
        return target
    width = 1200
    height = 800
    margin = 80
    xs = coords[:, 0]
    ys = coords[:, 1]
    min_x, max_x = float(xs.min()), float(xs.max())
    min_y, max_y = float(ys.min()), float(ys.max())
    span_x = max(max_x - min_x, 1e-9)
    span_y = max(max_y - min_y, 1e-9)
    colors = _color_map(clades)

    def scale_x(value: float) -> float:
        return margin + ((value - min_x) / span_x) * (width - (2 * margin))

    def scale_y(value: float) -> float:
        return height - margin - ((value - min_y) / span_y) * (height - (2 * margin))

    points = []
    for index, seq_id in enumerate(ids):
        x = scale_x(float(xs[index]))
        y = scale_y(float(ys[index]))
        clade = clades[index]
        points.append(
            f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.5" fill="{colors[clade]}" opacity="0.75" />'
            f'<title>{seq_id} ({clade})</title>'
        )

    legend = []
    for legend_index, clade in enumerate(sorted(colors)):
        y = margin + (legend_index * 24)
        legend.append(
            f'<circle cx="{width - 230}" cy="{y}" r="7" fill="{colors[clade]}" />'
            f'<text x="{width - 215}" y="{y + 5}" font-size="14">{clade}</text>'
        )

    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">'
        f'<rect width="{width}" height="{height}" fill="white" />'
        f'<text x="{margin}" y="40" font-size="24" font-weight="bold">{title}</text>'
        f'<line x1="{margin}" y1="{height - margin}" x2="{width - margin}" y2="{height - margin}" stroke="#999999" />'
        f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height - margin}" stroke="#999999" />'
        + "".join(points)
        + "".join(legend)
        + "</svg>"
    )
    target = Path(output_path)
    target.write_text(svg)
    return target


def visualize_profiles(
    input_table: PathLike,
    output_dir: PathLike,
    *,
    perplexities: Optional[Sequence[int]] = None,
    min_points: int = 20,
) -> dict[str, Path]:
    ids, matrix, clades, lengths = _load_hmm_score_table(input_table)
    normalized = _normalize_by_column_mean(matrix)
    n_rows = normalized.shape[0]
    target_perplexities = _sanitize_perplexities(perplexities, n_rows) if perplexities else _default_perplexities(n_rows)
    root = ensure_directory(output_dir)

    generated_tables: list[Path] = []
    generated_svgs: list[Path] = []
    summary_rows: list[dict[str, object]] = []
    for perplexity in target_perplexities:
        if n_rows < 3:
            coords = _pca_fallback(normalized)
        else:
            tsne = TSNE(
                n_components=2,
                perplexity=float(perplexity),
                random_state=42,
                init="pca",
                learning_rate="auto",
                max_iter=1000,
            )
            coords = tsne.fit_transform(normalized)
        clusters = _cluster_coordinates(coords, min_points=min_points)
        rows = []
        for index, seq_id in enumerate(ids):
            rows.append(
                {
                    "ID": seq_id,
                    "X": round(float(coords[index, 0]), 6),
                    "Y": round(float(coords[index, 1]), 6),
                    "cluster": int(clusters[index]),
                    "clade": clades[index],
                    "length": lengths[index],
                }
            )
        table_path = write_tsv(rows, root / f"tsne-db{perplexity}.tsv")
        svg_path = _render_tsne_svg(ids, coords, clades, root / f"tsne{perplexity}.svg", title=f"Ariadne TPS t-SNE (perplexity={perplexity})")
        generated_tables.append(table_path)
        generated_svgs.append(svg_path)
        summary_rows.append(
            {
                "perplexity": perplexity,
                "rows": len(rows),
                "table": str(table_path),
                "svg": str(svg_path),
            }
        )

    summary_path = write_tsv(summary_rows, root / "visualization_summary.tsv")
    return {
        "summary": summary_path,
        "tables_dir": root,
        "figures_dir": root,
    }
