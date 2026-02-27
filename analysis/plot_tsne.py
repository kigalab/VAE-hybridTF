#!/usr/bin/env python3
"""
t-SNE Plot Generator

Reads a CSV containing embedding vectors (e.g., 32D) and generates a 2D t-SNE plot.

- Points after the first N special rows are plotted as gray circles.
- The first N special rows are highlighted with distinct colors/markers.
  By default:
    row 1 -> "LuxR"
    row 2 -> "LasR"
    rows 3..N -> colored (hsv) markers

Requirements:
  numpy, pandas, matplotlib, scikit-learn

Usage:
  python plot_tsne.py dataset_info.csv --dim_start 1 --dim_end 33 --output tsne.svg

Author: Sota Okuda, Atsushi Minami
License: MIT
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE


DEFAULT_LUXR_COLOR = "#CC79A7"
DEFAULT_LASR_COLOR = "#009E73"


def _validate_slice(dim_start: int, dim_end: int) -> None:
    if dim_start < 0 or dim_end <= dim_start:
        raise ValueError(f"Invalid slice: dim_start={dim_start}, dim_end={dim_end}")


def _load_embedding_matrix(csv_path: Path, dim_start: int, dim_end: int) -> np.ndarray:
    """
    Load embedding vectors from CSV by slicing columns [dim_start:dim_end).
    Example: dim_start=1, dim_end=33 -> columns 1..32 (32D)
    """
    df = pd.read_csv(csv_path)
    if df.shape[1] < dim_end:
        raise ValueError(
            f"CSV has only {df.shape[1]} columns but dim_end={dim_end} was requested."
        )
    X = df.iloc[:, dim_start:dim_end].to_numpy(dtype=float)
    if np.any(~np.isfinite(X)):
        raise ValueError("Embedding matrix contains NaN/Inf values.")
    return X


def _run_tsne(
    X: np.ndarray,
    seed: int,
    perplexity: float,
    learning_rate: float | str,
    n_iter: int,
) -> np.ndarray:
    tsne = TSNE(
        n_components=2,
        perplexity=perplexity,
        random_state=seed,
        learning_rate=learning_rate,
        max_iter=n_iter,
        init="pca",
    )
    return tsne.fit_transform(X)


def _plot_tsne(
    X_2d: np.ndarray,
    output_path: Path,
    special_n: int,
    luxr_label: str,
    lasr_label: str,
    luxr_color: str,
    lasr_color: str,
    background_label: str,
    marker_special: str,
    marker_background: str,
    special_size: int,
    background_size: int,
    alpha_background: float,
    x_ticks: list[float] | None,
    y_ticks: list[float] | None,
) -> None:
    n = X_2d.shape[0]
    special_n = max(0, min(special_n, n))

    fig, ax = plt.subplots(figsize=(8, 6))

    # Background points (after special_n)
    if n > special_n:
        ax.scatter(
            X_2d[special_n:, 0],
            X_2d[special_n:, 1],
            s=background_size,
            alpha=alpha_background,
            marker=marker_background,
            color="gray",
            label=background_label,
        )

    # Colors for rows 3..special_n
    remaining_count = max(0, special_n - 2)
    remaining_colors = (
        plt.cm.hsv(np.linspace(0, 1, remaining_count, endpoint=False))
        if remaining_count > 0
        else []
    )

    for i in range(special_n):
        if i == 0:
            c = luxr_color
            lbl = luxr_label
        elif i == 1:
            c = lasr_color
            lbl = lasr_label
        else:
            c = remaining_colors[i - 2]
            lbl = None  # avoid huge legends for many highlighted points

        ax.scatter(
            X_2d[i, 0],
            X_2d[i, 1],
            s=special_size,
            marker=marker_special,
            color=c,
            label=lbl,
        )

    ax.set_xlabel("t-SNE dim 1")
    ax.set_ylabel("t-SNE dim 2")

    if x_ticks is not None:
        ax.set_xticks(x_ticks)
    if y_ticks is not None:
        ax.set_yticks(y_ticks)

    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, transparent=True)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a t-SNE plot from embeddings in CSV.")
    parser.add_argument("input_csv", type=str, help="Path to CSV containing embeddings.")
    parser.add_argument("--dim_start", type=int, default=1, help="Start column index (0-based, inclusive).")
    parser.add_argument("--dim_end", type=int, default=33, help="End column index (0-based, exclusive).")
    parser.add_argument("--output", type=str, default="tsne.svg", help="Output figure path (.svg or .png).")

    # t-SNE params
    parser.add_argument("--perplexity", type=float, default=50.0, help="t-SNE perplexity.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument(
        "--learning_rate",
        type=str,
        default="auto",
        help="t-SNE learning_rate (float or 'auto').",
    )
    parser.add_argument("--n_iter", type=int, default=1000, help="Max iterations (max_iter).")

    # Plot params
    parser.add_argument("--special_n", type=int, default=24, help="Number of special rows to highlight.")
    parser.add_argument("--luxr_label", type=str, default="LuxR", help="Legend label for first row.")
    parser.add_argument("--lasr_label", type=str, default="LasR", help="Legend label for second row.")
    parser.add_argument("--luxr_color", type=str, default=DEFAULT_LUXR_COLOR, help="Color for LuxR point.")
    parser.add_argument("--lasr_color", type=str, default=DEFAULT_LASR_COLOR, help="Color for LasR point.")
    parser.add_argument("--background_label", type=str, default="data points", help="Legend label for background points.")

    parser.add_argument("--marker_special", type=str, default="x", help="Marker for special points.")
    parser.add_argument("--marker_background", type=str, default="o", help="Marker for background points.")
    parser.add_argument("--special_size", type=int, default=80, help="Marker size for special points.")
    parser.add_argument("--background_size", type=int, default=20, help="Marker size for background points.")
    parser.add_argument("--alpha_background", type=float, default=0.5, help="Alpha for background points.")

    # Optional fixed ticks (match your original behavior)
    parser.add_argument("--x_ticks", type=str, default="-50,0,50", help="Comma-separated x ticks, or empty to disable.")
    parser.add_argument("--y_ticks", type=str, default="-50,0,50", help="Comma-separated y ticks, or empty to disable.")

    args = parser.parse_args()

    input_path = Path(args.input_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    _validate_slice(args.dim_start, args.dim_end)

    # Parse learning_rate
    if args.learning_rate.lower() == "auto":
        learning_rate: float | str = "auto"
    else:
        learning_rate = float(args.learning_rate)

    # Parse ticks
    def parse_ticks(s: str) -> list[float] | None:
        s = s.strip()
        if s == "":
            return None
        return [float(x) for x in s.split(",")]

    x_ticks = parse_ticks(args.x_ticks)
    y_ticks = parse_ticks(args.y_ticks)

    X = _load_embedding_matrix(input_path, args.dim_start, args.dim_end)
    print(f"Loaded embedding matrix: {X.shape} from {input_path}")

    X_2d = _run_tsne(
        X=X,
        seed=args.seed,
        perplexity=args.perplexity,
        learning_rate=learning_rate,
        n_iter=args.n_iter,
    )

    output_path = Path(args.output)
    _plot_tsne(
        X_2d=X_2d,
        output_path=output_path,
        special_n=args.special_n,
        luxr_label=args.luxr_label,
        lasr_label=args.lasr_label,
        luxr_color=args.luxr_color,
        lasr_color=args.lasr_color,
        background_label=args.background_label,
        marker_special=args.marker_special,
        marker_background=args.marker_background,
        special_size=args.special_size,
        background_size=args.background_size,
        alpha_background=args.alpha_background,
        x_ticks=x_ticks,
        y_ticks=y_ticks,
    )

    print(f"Done. Saved: {output_path.resolve()}")


if __name__ == "__main__":
    main()