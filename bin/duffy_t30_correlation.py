#!/usr/bin/env python3
"""
Compute and plot the correlation between:
  - Duffy group: Salmon riboseq gene-level TPM (SRR15175557, SRR15175558, SRR15175559)
  - t30 group:   STAR CN sample unstranded read counts, normalized to CPM

Usage:
    duffy_t30_correlation.py --salmon_dir <dir> --star_dir <dir> --output_dir <dir>
"""

import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr

DUFFY_SAMPLES = ["SRR15175557", "SRR15175558", "SRR15175559"]


def load_salmon_sample(quant_path: str) -> pd.Series:
    """
    Read quant.genes.sf, aggregate to true gene level, return gene TPM Series.

    quant.genes.sf may have multiple rows per gene because the Name field
    encodes transcript-level identifiers. We extract the Ensembl gene ID
    (2nd pipe-delimited field), sum NumReads across transcripts per gene,
    and recompute TPM from the summed rates.
    """
    df = pd.read_csv(quant_path, sep="\t")
    # Extract gene ID: "ENST...|ENSG00000123.4|-|..." -> "ENSG00000123.4"
    df["gene_id"] = df["Name"].str.split("|").str[1]
    # Compute per-transcript rate = NumReads / EffectiveLength
    df["rate"] = df["NumReads"] / df["EffectiveLength"].replace(0, float("nan"))
    # Sum rates and NumReads per gene
    gene_df = df.groupby("gene_id", as_index=True).agg(
        gene_rate=("rate", "sum"),
        gene_reads=("NumReads", "sum"),
    )
    # Recompute gene-level TPM
    total_rate = gene_df["gene_rate"].sum()
    gene_df["TPM"] = gene_df["gene_rate"] / total_rate * 1e6
    return gene_df["TPM"]


def load_duffy(salmon_dir: str) -> pd.Series:
    """Average gene TPM across the three Duffy samples."""
    frames = []
    for sample in DUFFY_SAMPLES:
        path = os.path.join(salmon_dir, sample, "quant.genes.sf")
        tpm = load_salmon_sample(path)
        tpm.name = sample
        frames.append(tpm)
    combined = pd.concat(frames, axis=1)
    return combined.mean(axis=1).rename("duffy_avg_tpm")


def load_star_sample(path: str) -> pd.Series:
    """
    Read a STAR ReadsPerGene.out.tab file, return unstranded counts Series.
    Skips the four summary rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous).
    """
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["gene_id", "unstranded", "stranded_fwd", "stranded_rev"],
    )
    df = df[~df["gene_id"].str.startswith("N_")].copy()
    df = df.set_index("gene_id")
    return df["unstranded"]


def load_t30(star_dir: str) -> pd.Series:
    """Average CPM across all CN_* samples."""
    pattern = os.path.join(star_dir, "CN_*ReadsPerGene.out.tab")
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(f"No CN_*ReadsPerGene.out.tab files found in {star_dir}")
    frames = []
    for path in paths:
        counts = load_star_sample(path)
        # Normalize to CPM
        total = counts.sum()
        cpm = counts / total * 1e6
        cpm.name = os.path.basename(path)
        frames.append(cpm)
    combined = pd.concat(frames, axis=1)
    return combined.mean(axis=1).rename("t30_avg_cpm")


def main():
    parser = argparse.ArgumentParser(description="Duffy vs t30 correlation scatter plot")
    parser.add_argument("--salmon_dir", required=True,
                        help="Directory containing salmon riboseq results (one subdir per sample)")
    parser.add_argument("--star_dir", required=True,
                        help="Directory containing STAR ReadsPerGene.out.tab files")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for the PDF plot")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading Duffy (salmon) data...")
    duffy = load_duffy(args.salmon_dir)

    print("Loading t30 (STAR CN) data...")
    t30 = load_t30(args.star_dir)

    # Inner join — only genes present in both datasets
    data = pd.concat([duffy, t30], axis=1, join="inner").dropna()
    print(f"Genes in common: {len(data)}")

    # Log10(x + 1) transform
    import numpy as np
    x = np.log10(data["duffy_avg_tpm"] + 1)
    y = np.log10(data["t30_avg_cpm"] + 1)

    # Spearman correlation
    rho, pval = spearmanr(x, y)
    print(f"Spearman r = {rho:.4f}, p = {pval:.2e}")

    # Scatter plot
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, s=4, alpha=0.3, linewidths=0, color="#2166ac")
    ax.set_xlabel("Duffy (Salmon riboseq)\nlog$_{10}$(avg TPM + 1)", fontsize=11)
    ax.set_ylabel("t30 CN (STAR RNA-seq)\nlog$_{10}$(avg CPM + 1)", fontsize=11)
    ax.set_title("Duffy vs. t30 expression correlation", fontsize=12)

    # Annotate correlation
    p_str = f"{pval:.2e}" if pval > 0 else "< 1e-300"
    ax.text(
        0.05, 0.95,
        f"$r_s$ = {rho:.3f}\n$p$ = {p_str}",
        transform=ax.transAxes,
        va="top", ha="left",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.8"),
    )

    plt.tight_layout()
    out_path = os.path.join(args.output_dir, "duffy_t30_correlation.pdf")
    fig.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
