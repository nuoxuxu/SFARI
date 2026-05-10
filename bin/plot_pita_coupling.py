#!/usr/bin/env python3
import argparse
from pathlib import Path
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde, ks_2samp, ttest_ind


TISSUE_MAP = {"iPSC": "t00", "NPC": "t04", "CN": "t30"}
TISSUE_ORDER = ["t00", "t04", "t30"]
COLORS = {"dualT": "#6b8cae", "soloT": "#b0b0b0"}


def sig_stars(pval):
    if pval < 1e-16:
        return "**"
    elif pval < 1e-8:
        return "*"
    else:
        return "ns"


def kde_curve(r_vals, n_grid=512, cut=3):
    """KDE using R-equivalent nrd0 bandwidth."""
    r_vals = r_vals[~np.isnan(r_vals)]
    if len(r_vals) < 2:
        return np.array([]), np.array([])
    std = np.std(r_vals, ddof=1)
    if std == 0:
        return np.array([]), np.array([])
    iqr = np.percentile(r_vals, 75) - np.percentile(r_vals, 25)
    bw = 0.9 * min(std, iqr / 1.34) * len(r_vals) ** (-1 / 5)
    x = np.linspace(r_vals.min() - cut * bw, r_vals.max() + cut * bw, n_grid)
    kde = gaussian_kde(r_vals, bw_method=bw / std)
    return x, kde(x)


def plot_A(fig, gs_A, df):
    """3 stacked KDE density plots, one per time point."""
    axes = [fig.add_subplot(gs_A[i]) for i in range(3)]

    for i, (ax, tissue) in enumerate(zip(axes, TISSUE_ORDER)):
        dual_r = df.filter(
            (pl.col("tissue") == tissue) & (pl.col("exon_type") == "dualT")
        )["spearman_r"].to_numpy()
        solo_r = df.filter(
            (pl.col("tissue") == tissue) & (pl.col("exon_type") == "soloT")
        )["spearman_r"].to_numpy()

        x_solo, y_solo = kde_curve(solo_r)
        x_dual, y_dual = kde_curve(dual_r)

        if len(x_solo):
            ax.plot(x_solo, y_solo, color=COLORS["soloT"], linewidth=1.5,
                    label="soloT")
            ax.fill_between(x_solo, y_solo, alpha=0.3, color=COLORS["soloT"])
        if len(x_dual):
            ax.plot(x_dual, y_dual, color=COLORS["dualT"], linewidth=1.5,
                    label="dualT")
            ax.fill_between(x_dual, y_dual, alpha=0.3, color=COLORS["dualT"])

        # K-S test
        stars = "ns"
        if len(dual_r) >= 2 and len(solo_r) >= 2:
            _, pval = ks_2samp(dual_r, solo_r)
            stars = sig_stars(pval)

        # Time point label on left y-axis; asterisk inside top-left
        ax.text(0.02, 0.93, stars, transform=ax.transAxes,
                va="top", ha="left", fontsize=9)
        ax.set_ylabel(tissue, fontsize=9, rotation=0, labelpad=13, va="center")
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        ax.axvline(0, color="black", linewidth=0.5, linestyle="--")

        if i == 0:
            # Panel label and legend on top sub-axis
            ax.text(-0.12, 1.05, "A", transform=ax.transAxes,
                    fontsize=12, fontweight="bold", va="top")
            ax.legend(
                handles=[Patch(facecolor=COLORS[k], alpha=0.6, label=k)
                         for k in ["soloT", "dualT"]],
                frameon=False, fontsize=8, loc="upper right"
            )
            ax.tick_params(labelbottom=False)
        elif i == 1:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel("Spearman's ρ")

    return axes


def plot_B(ax, df):
    """Horizontal bar chart of mean ΔAUC ± SEM per time point."""
    sample_dauc = df.group_by(["sample", "tissue", "exon_type"]).agg(
        pl.col("delta_auc").first()
    )

    y_positions = np.arange(len(TISSUE_ORDER) - 1, -1, -1)
    height = 0.35

    for exon_type, offset in [("soloT", height / 2), ("dualT", -height / 2)]:
        means, sems = [], []
        for tissue in TISSUE_ORDER:
            vals = (
                sample_dauc
                .filter((pl.col("tissue") == tissue) & (pl.col("exon_type") == exon_type))
                ["delta_auc"].drop_nulls().drop_nans().to_numpy()
            )
            if len(vals) == 0:
                means.append(0.0)
                sems.append(0.0)
            else:
                means.append(float(np.mean(vals)))
                sems.append(float(np.std(vals, ddof=1) / np.sqrt(len(vals))))

        ax.barh(
            y_positions + offset, means, height=height * 0.9,
            color=COLORS[exon_type], alpha=0.7, label=exon_type
        )
        ax.errorbar(
            means, y_positions + offset, xerr=sems,
            fmt="none", color="black", linewidth=1, capsize=3
        )

    ax.axvline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(TISSUE_ORDER)
    ax.set_xlabel("ΔAUC")
    ax.legend(frameon=False, fontsize=8)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.text(-0.12, 1.05, "B", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")


def plot_C(ax, df):
    """Box plot of ΔAUC across all samples with Welch's t-test bracket."""
    sample_dauc = df.group_by(["sample", "exon_type"]).agg(
        pl.col("delta_auc").first()
    )

    all_vals = {}
    for i, exon_type in enumerate(["soloT", "dualT"]):
        vals = (
            sample_dauc
            .filter(pl.col("exon_type") == exon_type)
            ["delta_auc"].drop_nulls().drop_nans().to_numpy()
        )
        all_vals[exon_type] = vals
        bp = ax.boxplot(
            vals, positions=[i], widths=0.5,
            patch_artist=True, showfliers=True,
            medianprops=dict(color="black", linewidth=1.5),
            flierprops=dict(
                marker="o", markersize=3, alpha=0.5,
                markerfacecolor=COLORS[exon_type], markeredgewidth=0
            )
        )
        bp["boxes"][0].set_facecolor(COLORS[exon_type])
        bp["boxes"][0].set_alpha(0.7)

    # Welch's t-test significance bracket
    solo_vals = all_vals["soloT"]
    dual_vals = all_vals["dualT"]
    if len(solo_vals) >= 2 and len(dual_vals) >= 2:
        _, pval = ttest_ind(dual_vals, solo_vals, equal_var=False)
        stars = sig_stars(pval)
        y_max = max(np.concatenate([solo_vals, dual_vals]))
        margin = (y_max - min(np.concatenate([solo_vals, dual_vals]))) * 0.08
        bracket_y = y_max + margin
        tick_len = margin * 0.4
        ax.plot([0, 0, 1, 1],
                [bracket_y - tick_len, bracket_y, bracket_y, bracket_y - tick_len],
                color="black", linewidth=1)
        ax.text(0.5, bracket_y + tick_len * 0.5, stars,
                ha="center", va="bottom", fontsize=11)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["soloT", "dualT"])
    ax.set_ylabel("ΔAUC")
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.text(-0.12, 1.05, "C", transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top")


def main():
    parser = argparse.ArgumentParser(description="Plot PITA coupling figures (A, B, C).")
    parser.add_argument("--input", nargs="+", required=True, help="Per-sample TSV files")
    parser.add_argument("--outdir", required=True, help="Output directory for PDF")
    args = parser.parse_args()

    df = pl.concat([pl.read_csv(f, separator="\t") for f in args.input])
    df = df.with_columns(
        pl.col("sample").str.split("_").list.first().replace(TISSUE_MAP).alias("tissue")
    )

    fig = plt.figure(figsize=(15, 5))
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.35)
    gs_A = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[0], hspace=0.15)
    ax_B = fig.add_subplot(gs[1])
    ax_C = fig.add_subplot(gs[2])

    axes_A = plot_A(fig, gs_A, df)
    plot_B(ax_B, df)
    plot_C(ax_C, df)

    # Caption for K-S test asterisk legend, placed below panel A
    # Get the bounding box of the bottom sub-axis in figure coordinates
    print("K-S test, *P-value < 10\u207b\u2078; **P-value < 10\u207b\u00b9\u2076")

    fig.savefig(Path(args.outdir) / "pita_coupling.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {Path(args.outdir) / 'pita_coupling.pdf'}")


if __name__ == "__main__":
    main()
