#!/usr/bin/env python3
"""
4D Scatter Plot Generator

Description:
Generates log-scaled 4D scatter plots from activity_summary_*.csv files.

X-axis  : geom_mean_CLasR
Y-axis  : geom_mean_CLuxR
Color   : geom_mean_{condition}
Size    : geom_mean_pET16

Wild-type sequences are highlighted with star markers.

Author: Atsushi Minami
License: MIT
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, LogFormatterSciNotation, NullLocator


# ===============================
# Configuration
# ===============================

CONDITIONS = ["20", "22", "1M"]

WT_LASR = "ACCTATCTCATTTGCTAGTT"
WT_LUXR = "ACCTGTAGGATCGTACAGGT"

MIN_COUNT_THRESHOLD = 20
UNIFIED_MIN = 90


# ===============================
# Utility Functions
# ===============================

def identify_wt(sequence: str) -> str:
    """Classify sequence as LasR(WT), LuxR(WT), or Non-WT."""
    if WT_LASR in str(sequence):
        return "LasR(WT)"
    if WT_LUXR in str(sequence):
        return "LuxR(WT)"
    return "Non-WT"


def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Apply read count filtering."""
    required_columns = [
        "total_CLasR", "total_CLuxR", "total_pET16",
        "total_20", "total_22", "total_1M"
    ]

    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"Missing required column: {col}")

    mask = np.ones(len(df), dtype=bool)
    for col in required_columns:
        mask &= df[col] >= MIN_COUNT_THRESHOLD

    df_filtered = df[mask].copy()
    df_filtered["WT_type"] = df_filtered["Sequence"].apply(identify_wt)

    return df_filtered


def compute_unified_max(df: pd.DataFrame) -> float:
    """Compute unified maximum for axis and color scaling."""
    values = [
        df["geom_mean_CLasR"].max(),
        df["geom_mean_CLuxR"].max(),
        df[["geom_mean_20", "geom_mean_22", "geom_mean_1M"]].max().max()
    ]
    return max(values)


# ===============================
# Plot Function
# ===============================

def create_plot(df, condition, output_path, sample_name):

    unified_max = compute_unified_max(df)

    x = df["geom_mean_CLasR"].values
    y = df["geom_mean_CLuxR"].values
    color_values = df[f"geom_mean_{condition}"].values
    size_raw = df["geom_mean_pET16"].values

    sizes = 40 / np.sqrt(size_raw / size_raw.min())

    mask_non_wt = df["WT_type"] == "Non-WT"
    mask_lasr = df["WT_type"] == "LasR(WT)"
    mask_luxr = df["WT_type"] == "LuxR(WT)"

    fig, ax = plt.subplots(figsize=(10, 9))

    scatter = ax.scatter(
        x[mask_non_wt], y[mask_non_wt],
        c=color_values[mask_non_wt],
        s=sizes[mask_non_wt],
        cmap="viridis",
        norm=LogNorm(vmin=UNIFIED_MIN, vmax=unified_max),
        alpha=0.7,
        edgecolors="black",
        linewidth=0.3,
        marker="o"
    )

    if mask_lasr.any():
        ax.scatter(
            x[mask_lasr], y[mask_lasr],
            c=color_values[mask_lasr],
            s=sizes[mask_lasr],
            cmap="viridis",
            norm=LogNorm(vmin=UNIFIED_MIN, vmax=unified_max),
            marker="*",
            edgecolors="black",
            linewidth=0.5
        )

    if mask_luxr.any():
        ax.scatter(
            x[mask_luxr], y[mask_luxr],
            c=color_values[mask_luxr],
            s=sizes[mask_luxr],
            cmap="viridis",
            norm=LogNorm(vmin=UNIFIED_MIN, vmax=unified_max),
            marker="*",
            edgecolors="black",
            linewidth=0.5
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(UNIFIED_MIN, unified_max)
    ax.set_ylim(UNIFIED_MIN, unified_max)

    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.xaxis.set_major_formatter(LogFormatterSciNotation())
    ax.yaxis.set_major_formatter(LogFormatterSciNotation())
    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())

    ax.set_xlabel("LasR Activity (geom_mean_CLasR)")
    ax.set_ylabel("LuxR Activity (geom_mean_CLuxR)")
    ax.set_title(f"4D Scatter Plot – {sample_name} (Condition: {condition})")

    cbar = plt.colorbar(scatter)
    cbar.set_label(f"geom_mean_{condition}")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


# ===============================
# Main
# ===============================

def main():

    parser = argparse.ArgumentParser(description="Generate 4D scatter plots.")
    parser.add_argument("input_csv", type=str, help="Path to activity_summary_*.csv")
    parser.add_argument("--output_dir", type=str, default="plots", help="Output directory")

    args = parser.parse_args()

    input_path = Path(args.input_csv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    sample_match = re.search(r"activity_summary_(.+)\.csv", input_path.name)
    sample_name = sample_match.group(1) if sample_match else "sample"

    df = pd.read_csv(input_path)
    df_filtered = filter_dataframe(df)

    if len(df_filtered) == 0:
        raise ValueError("No data remaining after filtering.")

    for condition in CONDITIONS:
        output_path = output_dir / f"4Dplot_{sample_name}_{condition}.png"
        create_plot(df_filtered, condition, output_path, sample_name)

    print("Finished successfully.")


if __name__ == "__main__":
    main()