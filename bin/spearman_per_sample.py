#!/usr/bin/env python3
import argparse
from pathlib import Path
import polars as pl
import pybedtools
from scipy.stats import spearmanr, gaussian_kde
import numpy as np


def read_bed(path, extra_column=None):
    if extra_column is None:
        return pl.read_csv(
            path,
            separator="\t", has_header=False, columns=range(6),
            new_columns=["chrom", "chromStart", "chromEnd", "name", "score", "strand"]
        )
    else:
        return pl.read_csv(
            path,
            separator="\t", has_header=False, columns=range(7),
            new_columns=["chrom", "chromStart", "chromEnd", "name", "score", "strand", extra_column]
        )


def get_intersected_df(a, b, prefix):
    intersected = pybedtools.BedTool.from_dataframe(a.to_pandas())\
        .intersect(
            pybedtools.BedTool.from_dataframe(b.to_pandas()),
            s=True, wo=True
        )
    out = (
        pl.from_pandas(intersected.to_dataframe(
            names=a.columns + ["_".join([prefix, col]) for col in b.columns] + ["overlap"]
        ))
        .drop(["_".join([prefix, col]) for col in b.columns if col != "name"])
        .filter(pl.col("overlap") == pl.col("overlap").max().over("name"))
    )
    return out


def spearman_agg(grp):
    result = spearmanr(grp["chromStart"], grp["chromEnd_right"])
    return pl.DataFrame({
        "gene_name": [grp["gene_name"][0]],
        "spearman_r": [result.statistic],
        "pvalue": [result.pvalue],
    })


def compute_delta_auc(r_vals, n_grid=512, cut=3):
    """Compute ΔAUC = AUC(ρ>0) - AUC(ρ<0) using R-equivalent nrd0 bandwidth."""
    r_vals = r_vals[~np.isnan(r_vals)]
    if len(r_vals) < 10:
        return float("nan")
    std = np.std(r_vals, ddof=1)
    if std == 0:
        return float("nan")
    iqr = np.percentile(r_vals, 75) - np.percentile(r_vals, 25)
    bw = 0.9 * min(std, iqr / 1.34) * len(r_vals) ** (-1 / 5)
    x = np.linspace(r_vals.min() - cut * bw, r_vals.max() + cut * bw, n_grid)
    kde = gaussian_kde(r_vals, bw_method=bw / std)
    density = kde(x)
    auc_pos = np.trapezoid(density[x >= 0], x[x >= 0])
    auc_neg = np.trapezoid(density[x <= 0], x[x <= 0])
    return float(auc_pos - auc_neg)


def main():
    parser = argparse.ArgumentParser(
        description="Compute per-gene Spearman r for dualT and soloT genes."
    )
    parser.add_argument("--bed", required=True, help="7-column BED file (last col = gene_name)")
    parser.add_argument("--cage_peaks", required=True, help="CAGE peaks BED file")
    parser.add_argument("--polya_sites", required=True, help="polyA sites BED file")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    sample = Path(args.bed).stem

    reads_all = read_bed(args.bed, extra_column="gene_name")
    cage_peaks = read_bed(args.cage_peaks)
    polya_sites = read_bed(args.polya_sites)

    first_exons = reads_all\
        .filter(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("chromStart") == pl.col("chromStart").min().over("name"))
            .otherwise(pl.col("chromEnd") == pl.col("chromEnd").max().over("name"))
        )\
        .pipe(get_intersected_df, cage_peaks, "CAGE")

    last_exons = reads_all\
        .filter(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("chromEnd") == pl.col("chromEnd").max().over("name"))
            .otherwise(pl.col("chromStart") == pl.col("chromStart").min().over("name"))
        )\
        .pipe(get_intersected_df, polya_sites, "polyA")

    reads_with_both = (
        first_exons[["chrom", "chromStart", "chromEnd", "name", "strand", "gene_name", "CAGE_name"]]
        .join(last_exons[["chromStart", "chromEnd", "name", "polyA_name"]], on="name", how="inner")
    )

    gene_counts = reads_with_both.group_by("gene_name").agg([
        pl.col("CAGE_name").n_unique().alias("n_CAGE"),
        pl.col("polyA_name").n_unique().alias("n_polyA"),
    ])

    dual_t_genes = gene_counts.filter(
        (pl.col("n_CAGE") > 1) & (pl.col("n_polyA") > 1)
    )["gene_name"]

    solo_t_genes = gene_counts.filter(
        (pl.col("n_CAGE") == 1) | (pl.col("n_polyA") == 1)
    )["gene_name"]

    results = []
    for exon_type, gene_list in [("dualT", dual_t_genes), ("soloT", solo_t_genes)]:
        df = (
            reads_with_both
            .filter(pl.col("gene_name").is_in(gene_list))
            .with_columns(pl.col("chromStart").sort().over("gene_name"))
            .group_by("gene_name").map_groups(spearman_agg)
            .drop_nans("spearman_r")
        )
        r_vals = df["spearman_r"].to_numpy()
        delta_auc = compute_delta_auc(r_vals)
        df = df.with_columns([
            pl.lit(exon_type).alias("exon_type"),
            pl.lit(sample).alias("sample"),
            pl.lit(delta_auc).alias("delta_auc"),
        ])
        results.append(df)

    pl.concat(results).write_csv(args.out, separator="\t")
    print(f"Written: {args.out}")


if __name__ == "__main__":
    main()
