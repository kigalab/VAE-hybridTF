#!/usr/bin/env python3
"""
Scatter plot for library screening (Plux vs Plas).

Expected input CSV columns (default):
  - name : variant name (string)
  - plux : x value (float)
  - plas : y value (float)

Color rule (default):
  - if "LuxR" in name -> pink
  - if "LasR" in name -> green
  - if name == "Ctrl" -> gray
  - else -> black

Requirements:
  pandas
  matplotlib

Usage:
  python scatter_library.py input.csv --output scatter.svg

Author: Atsushi Minami
License: MIT
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator


DEFAULT_LUXR_COLOR = "#CC79A7"
DEFAULT_LASR_COLOR = "#009E73"
DEFAULT_CTRL_COLOR = "#B3B3B3"


def get_color(name: str) -> str:
    """Default color rule based on variant name."""
    if "LuxR" in name:
        return DEFAULT_LUXR_COLOR
    if "LasR" in name:
        return DEFAULT_LASR_COLOR
    if name == "Ctrl":
        return DEFAULT_CTRL_COLOR
    return "black"


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot Plux vs Plas scatter from a CSV.")
    parser.add_argument("input_csv", type=str, help="Input CSV path.")
    parser.add_argument("--output", type=str, default="scatter.svg", help="Output figure path (.svg/.png).")
    parser.add_argument("--name_col", type=str, default="name", help="Column name for labels.")
    parser.add_argument("--x_col", type=str, default="plux", help="Column name for x values.")
    parser.add_argument("--y_col", type=str, default="plas", help="Column name for y values.")
    parser.add_argument("--xlim", type=str, default="1e2,1e5", help="x-axis limits 'min,max' (empty to disable).")
    parser.add_argument("--ylim", type=str, default="1e2,1e5", help="y-axis limits 'min,max' (empty to disable).")
    parser.add_argument("--opaque", action="store_true", help="Save with non-transparent background.")
    args = parser.parse_args()

    input_path = Path(args.input_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    df = pd.read_csv(input_path)

    for col in (args.name_col, args.x_col, args.y_col):
        if col not in df.columns:
            raise ValueError(f"Missing required column '{col}' in CSV. Available: {list(df.columns)}")

    fig, ax = plt.subplots(figsize=(6, 5))

    for _, row in df.iterrows():
        name = str(row[args.name_col])
        x = float(row[args.x_col])
        y = float(row[args.y_col])
        ax.scatter(x, y, color=get_color(name))

    ax.set_xscale("log")
    ax.set_yscale("log")

    # Remove minor ticks (match original style)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())

    ax.set_xlabel("Plux activity")
    ax.set_ylabel("Plas activity")

    def parse_lim(s: str):
        s = s.strip()
        if s == "":
            return None
        a, b = s.split(",")
        return float(a), float(b)

    xlim = parse_lim(args.xlim)
    ylim = parse_lim(args.ylim)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    fig.tight_layout()
    fig.savefig(Path(args.output), transparent=not args.opaque, dpi=300)
    plt.close(fig)

    print(f"Done. Saved: {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()