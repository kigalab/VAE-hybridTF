#!/usr/bin/env python3
"""
Panels for randomized promoter libraries.

Panel b:
  Violin plots of activity scores across TFs (log10 scale),
  with group means (horizontal bars),
  and one-way ANOVA followed by Tukey HSD.

Panel c:
  PLS regression (one-hot encoded randomized positions -> TF activities),
  then plot PLS coefficients as 2x2 heatmaps (base × position),
  saved as: PLS_coef_2x2_{library_name}_revaxes.svg

Input CSV requirements:
  - seq column (default: "seq") : promoter sequence string
  - NC column  (default: "NC")  : normalization/count column
  - TF activity columns (default: LuxRWT, LasRWT, TF20, TF22)

Preprocessing (same logic as your original script):
  - filter rows with NC >= nc_max (default 2000)
  - normalize each TF activity by dividing by NC
  - Y for PLS: log1p(Y_raw) then StandardScaler

Dependencies:
  numpy, pandas, matplotlib, seaborn, scipy, statsmodels, scikit-learn

Usage:
  python randomlib.py Plas_random.csv --library_name Plas --outdir figs

Author: Atsushi Minami
License: MIT
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler


BASES = ["A", "C", "G", "T"]


def validate_columns(df: pd.DataFrame, seq_col: str, nc_col: str, tf_cols: List[str]) -> None:
    required = [seq_col, nc_col] + tf_cols
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Found: {list(df.columns)}")


def nc_filter_and_divide(df: pd.DataFrame, nc_col: str, tf_cols: List[str], nc_max: float) -> pd.DataFrame:
    out = df.copy()
    n0 = len(out)
    out = out[out[nc_col] < nc_max].copy()
    print(f"NC filter: removed {n0 - len(out)} rows with {nc_col} >= {nc_max}. Remaining: {len(out)}")

    # Divide-normalize each TF activity by NC
    for c in tf_cols:
        out[c] = out[c].astype(float) / out[nc_col].astype(float)

    return out


def detect_variable_positions(seqs: np.ndarray) -> List[int]:
    lengths = np.array([len(s) for s in seqs])
    if len(set(lengths)) != 1:
        raise ValueError(f"Sequences have different lengths: {set(lengths)}")
    L = int(lengths[0])

    var_pos = []
    for i in range(L):
        if len(set(s[i] for s in seqs)) > 1:
            var_pos.append(i)

    if not var_pos:
        raise ValueError("No variable positions detected. Check sequences / library design.")
    return var_pos


def onehot_encode(seqs: np.ndarray, variable_positions: List[int]) -> np.ndarray:
    base_to_idx: Dict[str, int] = {b: i for i, b in enumerate(BASES)}

    def encode(seq: str) -> np.ndarray:
        vec = np.zeros(len(variable_positions) * 4, dtype=float)
        for pi, pos in enumerate(variable_positions):
            b = seq[pos]
            if b in base_to_idx:
                vec[pi * 4 + base_to_idx[b]] = 1.0
            # non-ACGT remains all-zero for that position
        return vec

    X = np.vstack([encode(s) for s in seqs])
    return X


def panel_b_violin_anova_tukey(df: pd.DataFrame, tf_cols: List[str], outpath: Path, library_name: str) -> None:
    # Long format
    melt = df[tf_cols].melt(var_name="TF", value_name="ActivityNorm")
    melt["log10Activity"] = np.log10(melt["ActivityNorm"] + 1e-9)

    # Plot
    plt.figure(figsize=(7.5, 5.5))
    ax = sns.violinplot(x="TF", y="log10Activity", data=melt, inner=None, cut=0)

    # Mean bars
    means = melt.groupby("TF")["log10Activity"].mean()
    for i, tf in enumerate(tf_cols):
        ax.hlines(means[tf], i - 0.25, i + 0.25, colors="black", linewidth=2)

    ax.set_xlabel("")
    ax.set_ylabel("log10(Activity / NC)")
    ax.set_title(f"Activity distributions ({library_name})")

    plt.tight_layout()
    plt.savefig(outpath, dpi=300, transparent=True)
    plt.close()

    # Stats (ANOVA on log10 values)
    groups = [melt.loc[melt["TF"] == tf, "log10Activity"].to_numpy() for tf in tf_cols]
    anova = f_oneway(*groups)

    # Tukey HSD post hoc
    tukey = pairwise_tukeyhsd(endog=melt["log10Activity"], groups=melt["TF"], alpha=0.05)

    stats_txt = outpath.with_suffix(".stats.txt")
    with open(stats_txt, "w", encoding="utf-8") as f:
        f.write("Panel b statistics\n")
        f.write(f"One-way ANOVA on log10(Activity/NC): F={anova.statistic:.6g}, p={anova.pvalue:.6g}\n\n")
        f.write("Tukey HSD (alpha=0.05)\n")
        f.write(tukey.summary().as_text())
        f.write("\n")

    tukey_csv = outpath.with_suffix(".tukey.csv")
    tukey_df = pd.DataFrame(tukey.summary().data[1:], columns=tukey.summary().data[0])
    tukey_df.to_csv(tukey_csv, index=False)

    print(f"Saved panel b: {outpath}")
    print(f"Saved stats:   {stats_txt}")
    print(f"Saved Tukey:   {tukey_csv}")


def panel_c_pls_coef_heatmap(
    df: pd.DataFrame,
    seq_col: str,
    tf_cols: List[str],
    outpath: Path,
    library_name: str,
    n_components: int = 2,
    symmetric_scale: bool = True,
) -> None:
    seqs = df[seq_col].astype(str).str.upper().to_numpy()

    variable_positions = detect_variable_positions(seqs)
    X = onehot_encode(seqs, variable_positions)

    # Y_raw = normalized activity (already divided by NC)
    Y_raw = df[tf_cols].astype(float).to_numpy()

    # Drop rows with NaN/inf in Y
    mask = np.isfinite(Y_raw).all(axis=1)
    if not mask.all():
        dropped = int((~mask).sum())
        print(f"Warning: dropping {dropped} rows with NaN/inf in Y.")
        X = X[mask]
        Y_raw = Y_raw[mask]

    # log1p + standardize (same idea as original)
    Y_log = np.log1p(Y_raw)
    Y = StandardScaler().fit_transform(Y_log)

    pls = PLSRegression(n_components=n_components)
    pls.fit(X, Y)

    coef = pls.coef_  # could be (n_features, n_targets) or transposed depending on sklearn version
    n_feat = X.shape[1]
    n_targets = len(tf_cols)

    if coef.shape == (n_feat, n_targets):
        coef_feat_first = coef
    elif coef.shape == (n_targets, n_feat):
        coef_feat_first = coef.T
    else:
        raise ValueError(
            f"Unexpected coef shape: {coef.shape}, expected ({n_feat},{n_targets}) or ({n_targets},{n_feat})"
        )

    n_pos = len(variable_positions)

    # color scale
    if symmetric_scale:
        lim = float(np.nanmax(np.abs(coef_feat_first)))
        vmin, vmax = -lim, +lim
    else:
        vmin, vmax = float(np.nanmin(coef_feat_first)), float(np.nanmax(coef_feat_first))

    # 2x2 panel (assumes 4 TFs; if not, remaining panels are hidden)
    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    axes = axes.flatten()

    for t_idx, tf_name in enumerate(tf_cols):
        ax = axes[t_idx]

        # reshape as (base, position): (4, n_pos)
        mat = coef_feat_first[:, t_idx].reshape(n_pos, 4).T

        sns.heatmap(
            mat,
            ax=ax,
            cmap="coolwarm",
            center=0.0 if symmetric_scale else None,
            vmin=vmin,
            vmax=vmax,
            xticklabels=[p + 1 for p in variable_positions],  # 1-based positions
            yticklabels=BASES,
            cbar=True,
            cbar_kws={"label": "PLS coefficient"},
        )
        ax.set_xlabel("Position (1-based)")
        ax.set_ylabel("Base")
        ax.set_title(tf_name)

    for j in range(len(tf_cols), len(axes)):
        axes[j].axis("off")

    plt.suptitle(f"PLS coefficients (base × position) ({library_name})", y=1.03)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, transparent=True, bbox_inches="tight")
    plt.close()

    print(f"Saved panel c: {outpath}")


def main() -> None:
    p = argparse.ArgumentParser(description="Panel b (violin+stats) and panel c (PLS coef heatmaps).")
    p.add_argument("input_csv", type=str, help="Input CSV for randomized promoter library.")
    p.add_argument("--library_name", type=str, default="Plas", help="Used for output names/titles.")
    p.add_argument("--outdir", type=str, default="figs", help="Output directory.")

    p.add_argument("--seq_col", type=str, default="seq")
    p.add_argument("--nc_col", type=str, default="NC")
    p.add_argument("--tf_cols", type=str, default="LuxRWT,LasRWT,TF20,TF22",
                   help="Comma-separated TF columns (must match CSV).")

    p.add_argument("--nc_max", type=float, default=2000.0)
    p.add_argument("--pls_components", type=int, default=2)
    p.add_argument("--no_symmetric_scale", action="store_true",
                   help="If set, use min/max scale instead of symmetric +/- scale for heatmaps.")

    args = p.parse_args()

    input_path = Path(args.input_csv)
    if not input_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    tf_cols = [x.strip() for x in args.tf_cols.split(",") if x.strip()]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path)
    validate_columns(df, args.seq_col, args.nc_col, tf_cols)

    df = nc_filter_and_divide(df, args.nc_col, tf_cols, args.nc_max)

    # Panel b
    b_out = outdir / f"b_violin_log10_activity_{args.library_name}.svg"
    panel_b_violin_anova_tukey(df, tf_cols, b_out, args.library_name)

    # Panel c (this is the one you want)
    c_out = outdir / f"PLS_coef_2x2_{args.library_name}_revaxes.svg"
    panel_c_pls_coef_heatmap(
        df=df,
        seq_col=args.seq_col,
        tf_cols=tf_cols,
        outpath=c_out,
        library_name=args.library_name,
        n_components=args.pls_components,
        symmetric_scale=not args.no_symmetric_scale,
    )

    print(f"Done. Outputs in: {outdir.resolve()}")


if __name__ == "__main__":
    main()