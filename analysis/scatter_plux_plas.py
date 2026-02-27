#!/usr/bin/env python3
"""
Scatter plot (Plux vs Plas) with mean ± SEM error bars.

Input CSV format (required):

condition,promoter,rep1,rep2,rep3,...
Negative,Plux,,,
Negative,Plas,,,
LuxR,Plux,,,
LuxR,Plas,,,
...

Each condition must contain both:
  - Plux
  - Plas

Requirements:
  numpy
  pandas
  matplotlib

Usage:
  python scatter_plux_plas.py input.csv --output scatter.svg

Author: Atsushi Minami
License: MIT
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator


@dataclass(frozen=True)
class SummaryPoint:
    condition: str
    plux_mean: float
    plas_mean: float
    plux_sem: float
    plas_sem: float


def mean_sem(values: List[float]) -> Tuple[float, float]:
    """Return mean and SEM (sample standard deviation, ddof=1)."""
    arr = np.asarray(values, dtype=float)
    if arr.size < 2:
        return float(arr.mean()), 0.0
    sem = float(arr.std(ddof=1) / np.sqrt(arr.size))
    return float(arr.mean()), sem


def load_from_csv(csv_path: Path) -> Dict[str, Dict[str, List[float]]]:
    """
    Load replicate data from CSV.

    Expected columns:
      condition, promoter, rep1, rep2, ...
    """
    df = pd.read_csv(csv_path)

    required = {"condition", "promoter"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {sorted(required)}")

    rep_cols = [c for c in df.columns if c not in ("condition", "promoter")]
    if len(rep_cols) == 0:
        raise ValueError("CSV must include at least one replicate column.")

    data: Dict[str, Dict[str, List[float]]] = {}

    for _, row in df.iterrows():
        cond = str(row["condition"])
        prom = str(row["promoter"])

        reps = [float(row[c]) for c in rep_cols if pd.notna(row[c])]

        if cond not in data:
            data[cond] = {}

        data[cond][prom] = reps

    return data


def summarize(data: Dict[str, Dict[str, List[float]]]) -> List[SummaryPoint]:
    """Convert replicate data into mean ± SEM summary points."""
    points: List[SummaryPoint] = []

    for condition, by_promoter in data.items():
        if "Plux" not in by_promoter or "Plas" not in by_promoter:
            raise ValueError(
                f"Condition '{condition}' must contain both Plux and Plas."
            )

        plux_mean, plux_sem = mean_sem(by_promoter["Plux"])
        plas_mean, plas_sem = mean_sem(by_promoter["Plas"])

        points.append(
            SummaryPoint(
                condition=condition,
                plux_mean=plux_mean,
                plas_mean=plas_mean,
                plux_sem=plux_sem,
                plas_sem=plas_sem,
            )
        )

    return points


def plot_scatter(
    points: List[SummaryPoint],
    output: Path,
    xlim: Tuple[float, float] | None,
    ylim: Tuple[float, float] | None,
    xlabel: str,
    ylabel: str,
    transparent: bool,
) -> None:
    """Create scatter plot with error bars."""
    fig, ax = plt.subplots(figsize=(6, 5))

    for p in points:
        ax.errorbar(
            p.plux_mean,
            p.plas_mean,
            xerr=p.plux_sem,
            yerr=p.plas_sem,
            fmt="o",
            color="black",
            ecolor="black",
            capsize=4,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    fig.tight_layout()
    fig.savefig(output, transparent=transparent, dpi=300)
    plt.close(fig)


def parse_lim(s: str) -> Tuple[float, float] | None:
    s = s.strip()
    if s == "":
        return None
    a, b = s.split(",")
    return float(a), float(b)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot Plux vs Plas scatter with mean ± SEM from CSV."
    )
    parser.add_argument("input_csv", type=str, help="Input CSV file (required).")
    parser.add_argument("--output", type=str, default="scatter.svg")
    parser.add_argument("--xlim", type=str, default="1000,1000000")
    parser.add_argument("--ylim", type=str, default="1000,100000")
    parser.add_argument("--xlabel", type=str, default="Plux activity")
    parser.add_argument("--ylabel", type=str, default="Plas activity")
    parser.add_argument("--opaque", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    data = load_from_csv(input_path)
    points = summarize(data)

    plot_scatter(
        points=points,
        output=Path(args.output),
        xlim=parse_lim(args.xlim),
        ylim=parse_lim(args.ylim),
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        transparent=not args.opaque,
    )

    print(f"Done. Saved: {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()