#!/usr/bin/env python3
"""
Grid Search Result Plotter

Reads a grid_search_result.csv and generates:
1) Sorted model performance (rank vs correlation)
2) Correlation vs latent_dim (log2 x-axis) with jitter
3) Correlation vs coef_kld (log x-axis) with jitter

Requirements:
  - numpy
  - pandas
  - matplotlib

Usage:
  python plot_grid_search_results.py path/to/grid_search_result.csv --output_dir plots

Author: Sota Okuda, Atsushi Minami
License: MIT
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DEFAULT_TARGET_COL = "combination_2_seqs_aa_pair_frequency_coef"


def _ensure_columns(df: pd.DataFrame, required: list[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required column(s): {missing}")


def plot_sorted_rank(
    df: pd.DataFrame,
    target_col: str,
    outpath: Path,
    threshold: float | None = 0.8,
) -> None:
    """Plot model rank (sorted by target_col desc) vs target_col."""
    df_sorted = df.sort_values(by=target_col, ascending=False).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(
        df_sorted.index.to_numpy() + 1,
        df_sorted[target_col].to_numpy(),
        marker="o",
        linestyle="None",
        markersize=3,
        color="black",
        label="Model Performance",
    )

    ax.set_xlabel("Model Rank (highest to lowest)")
    ax.set_ylabel("Correlation Coefficient")
    ax.grid(False)

    if threshold is not None:
        ax.axhline(
            y=threshold,
            color="red",
            linestyle="--",
            alpha=0.5,
            label=f"Threshold ({threshold})",
        )
        ax.legend()

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, transparent=True)
    plt.close(fig)


def plot_vs_latent_dim(
    df: pd.DataFrame,
    target_col: str,
    x_col: str,
    outpath: Path,
    seed: int = 0,
    jitter_strength: float = 0.2,
) -> None:
    """Scatter plot of target_col vs latent_dim (log2 x-axis) with multiplicative jitter."""
    rng = np.random.default_rng(seed)
    jitter = rng.uniform(-0.1, 0.1, size=len(df))
    x = df[x_col].to_numpy(dtype=float)
    x_jittered = x * (1 + jitter * jitter_strength)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(
        x_jittered,
        df[target_col].to_numpy(dtype=float),
        alpha=0.6,
        s=25,
        label="Models",
    )

    ax.set_xlabel("Latent Dimension")
    ax.set_ylabel(f"Correlation Coefficient ({target_col})")

    # Show powers-of-two spacing evenly while labeling actual dims
    ax.set_xscale("log", base=2)
    dims = sorted(df[x_col].unique())
    ax.set_xticks(dims)
    ax.set_xticklabels([str(d) for d in dims])

    ax.grid(False)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, transparent=True)
    plt.close(fig)


def plot_vs_coef_kld(
    df: pd.DataFrame,
    target_col: str,
    x_col: str,
    outpath: Path,
    seed: int = 0,
) -> None:
    """Scatter plot of target_col vs coef_kld (log x-axis) with multiplicative jitter."""
    rng = np.random.default_rng(seed)
    jitter = rng.uniform(0.85, 1.15, size=len(df))

    x = df[x_col].to_numpy(dtype=float)
    x_jittered = x * jitter

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(
        x_jittered,
        df[target_col].to_numpy(dtype=float),
        alpha=0.6,
        s=30,
        label="Models",
    )

    ax.set_xscale("log")
    ticks = sorted(df[x_col].unique())
    ax.set_xticks(ticks)

    # Readable tick labels (powers of 10 where appropriate)
    labels = []
    for t in ticks:
        t = float(t)
        if t < 0.1 or t > 10:
            labels.append(rf"$10^{{{int(np.log10(t))}}}$")
        else:
            labels.append(str(t))
    ax.set_xticklabels(labels)

    ax.minorticks_off()
    ax.set_xlabel("KLD Coefficient (coef_kld)")
    ax.set_ylabel("Correlation Coefficient")
    ax.grid(False)
    ax.legend()

    fig.tight_layout()
    fig.savefig(outpath, dpi=300, transparent=True)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot grid search results from CSV.")
    parser.add_argument("input_csv", type=str, help="Path to grid_search_result.csv")
    parser.add_argument("--output_dir", type=str, default="plots", help="Directory for output figures")
    parser.add_argument("--target_col", type=str, default=DEFAULT_TARGET_COL, help="Target metric column")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for jitter (reproducibility)")
    parser.add_argument("--format", type=str, default="svg", choices=["svg", "png"], help="Output image format")
    parser.add_argument("--threshold", type=float, default=0.8, help="Threshold line for rank plot")
    args = parser.parse_args()

    input_path = Path(args.input_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path)

    required = [args.target_col, "latent_dim", "coef_kld"]
    _ensure_columns(df, required)

    suffix = args.format.lower()

    plot_sorted_rank(
        df=df,
        target_col=args.target_col,
        outpath=outdir / f"sorted_rank.{suffix}",
        threshold=args.threshold,
    )
    plot_vs_latent_dim(
        df=df,
        target_col=args.target_col,
        x_col="latent_dim",
        outpath=outdir / f"correlation_vs_latent_dim.{suffix}",
        seed=args.seed,
    )
    plot_vs_coef_kld(
        df=df,
        target_col=args.target_col,
        x_col="coef_kld",
        outpath=outdir / f"correlation_vs_coef_kld.{suffix}",
        seed=args.seed,
    )

    print(f"Done. Figures saved to: {outdir.resolve()}")


if __name__ == "__main__":
    main()