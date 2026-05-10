#!/usr/bin/env python3
"""
For each transcript discovery tool, compute:
  1. % of transcripts with all splice junctions supported by short-read RNA-seq
  2. % of transcripts with 5' end within 100 bp of a CAGE-seq peak
Plot a grouped barplot comparing tools.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent / "src"))
from utils import read_gtf, read_SJ


def expressed_tx(expr_df: pl.DataFrame, min_reads: int = 5, min_samples: int = 2) -> set:
    id_col = expr_df.columns[0]
    sample_cols = expr_df.columns[1:]
    return set(
        expr_df.filter(
            pl.sum_horizontal([pl.col(c) > min_reads for c in sample_cols]) > min_samples
        )[id_col].to_list()
    )


def build_sr_sj_set(sj_files: list[Path]) -> pl.DataFrame:
    sr_sj = (
        pl.concat([read_SJ(f) for f in sj_files])
        .filter((pl.col("unique_reads") + pl.col("multi_reads")) > 0)
        .with_columns(
            pl.col("strand").map_elements(
                lambda s: "+" if s == 1 else ("-" if s == 2 else "."),
                return_dtype=pl.String,
            )
        )
        .select(["chrom", "start", "end", "strand"])
        .unique()
    )
    return sr_sj


def build_cage_lookup(cage_bed: Path) -> dict:
    cage = (
        pl.read_csv(
            cage_bed, separator="\t", has_header=False,
            schema_overrides={"column_1": pl.String},
        )
        .select(
            pl.col("column_1").alias("chrom"),
            pl.col("column_2").alias("start"),
            pl.col("column_3").alias("end"),
            pl.col("column_6").alias("strand"),
        )
    )
    lookup: dict = {}
    for (chrom, strand), grp in cage.group_by(["chrom", "strand"]):
        starts = grp["start"].to_numpy()
        ends = grp["end"].to_numpy()
        order = np.argsort(starts)
        lookup[(chrom, strand)] = (starts[order], ends[order])
    return lookup


def within_100bp(tss_pos: int, starts: np.ndarray, ends: np.ndarray) -> bool:
    idx = np.searchsorted(starts, tss_pos + 101)
    if idx == 0:
        return False
    return bool(np.any(ends[:idx] >= tss_pos - 100))


def compute_polya_support(transcripts: pl.DataFrame, polya: pl.DataFrame) -> tuple[int, float]:
    three_prime_df = (
        transcripts
        .with_columns(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("end"))
            .otherwise(pl.col("start"))
            .alias("pos")
        )
        .select([
            pl.col("transcript_id"),
            pl.col("seqname"),
            pl.col("pos"),
        ])
    )
    n_polya = (
        three_prime_df
        .join_where(
            polya,
            (pl.col("pos") >= pl.col("start")) &
            (pl.col("pos") <= pl.col("end")),
        )
        .filter(pl.col("seqname") == pl.col("chrom"))
        .select("transcript_id")
        .unique()
        .height
    )
    return n_polya, n_polya / transcripts.height * 100


def compute_support(
    gtf_path: Path,
    sr_sj: pl.DataFrame,
    cage_lookup: dict,
    expr_ids: set | None = None,
    polya: pl.DataFrame | None = None,
) -> tuple[float, float, float | None]:
    gtf = read_gtf(gtf_path, attributes=["transcript_id"])
    exons = gtf.filter(pl.col("feature") == "exon")
    transcripts = gtf.filter(pl.col("feature") == "transcript")
    all_ids = transcripts["transcript_id"].unique()
    total_raw = len(all_ids)
    if expr_ids is not None:
        expr_list = list(expr_ids)
        transcripts = transcripts.filter(pl.col("transcript_id").is_in(expr_list))
        exons = exons.filter(pl.col("transcript_id").is_in(expr_list))
        all_ids = transcripts["transcript_id"].unique()
        print(f"  Transcripts in GTF: {total_raw}  →  after expression filter: {len(all_ids)}")
    else:
        print(f"  Transcripts in GTF: {total_raw}")
    total = len(all_ids)

    # ── SR splice junction support ─────────────────────────────────────────────
    lr_junctions = (
        exons
        .sort(["transcript_id", "start"])
        .with_columns(
            pl.col("start").shift(-1).over("transcript_id").alias("next_exon_start"),
        )
        .with_columns(
            (pl.col("next_exon_start") - 1).alias("sj_end"),
            (pl.col("end") + 1).alias("sj_start"),
        )
        .filter(pl.col("next_exon_start").is_not_null())
        .select([
            pl.col("transcript_id"),
            pl.col("seqname").alias("chrom"),
            pl.col("sj_start").alias("start"),
            pl.col("sj_end").alias("end"),
            pl.col("strand"),
        ])
    )

    unsupported_ids = (
        lr_junctions
        .join(sr_sj, on=["chrom", "start", "end", "strand"], how="anti")
        .select("transcript_id")
        .unique()
    )
    n_sr = total - len(
        all_ids.filter(all_ids.is_in(unsupported_ids["transcript_id"].to_list()))
    )
    pct_sr = n_sr / total * 100
    print(f"  SR-supported: {n_sr} ({pct_sr:.1f}%)")

    # ── CAGE support ───────────────────────────────────────────────────────────
    tss_df = (
        transcripts
        .with_columns(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("start"))
            .otherwise(pl.col("end"))
            .alias("tss_pos")
        )
        .select([
            pl.col("transcript_id"),
            pl.col("seqname").alias("chrom"),
            pl.col("strand"),
            pl.col("tss_pos"),
        ])
    )

    cage_flags = []
    for row in tss_df.iter_rows(named=True):
        key = (row["chrom"], row["strand"])
        if key not in cage_lookup:
            cage_flags.append(False)
        else:
            starts, ends = cage_lookup[key]
            cage_flags.append(within_100bp(row["tss_pos"], starts, ends))

    n_cage = sum(cage_flags)
    pct_cage = n_cage / total * 100
    print(f"  CAGE-supported: {n_cage} ({pct_cage:.1f}%)")

    # ── PolyA cluster support ──────────────────────────────────────────────────
    pct_polya = None
    if polya is not None:
        n_polya, pct_polya = compute_polya_support(transcripts, polya)
        print(f"  PolyA-supported: {n_polya} ({pct_polya:.1f}%)")

    return pct_sr, pct_cage, pct_polya


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", nargs=2, metavar=("LABEL", "GTF"),
                        action="append", required=True,
                        help="Tool label and GTF path (repeatable)")
    parser.add_argument("--sj", nargs="+", required=True,
                        help="STAR SJ.out.tab files")
    parser.add_argument("--cage", required=True,
                        help="CAGE BED file (refTSS)")
    parser.add_argument("--bambu-expression", metavar="CSV",
                        help="Bambu expression CSV for filtering")
    parser.add_argument("--isoseq-classification", metavar="TXT",
                        help="IsoSeq filtered_lite_classification.txt for pre-filtering")
    parser.add_argument("--isoseq-expression", metavar="PARQUET",
                        help="IsoSeq full_expression.parquet for expression filtering")
    parser.add_argument("--polya", metavar="BED",
                        help="polyA cluster BED file (atlas.clusters.2.0.GRCh38.96.bed)")
    parser.add_argument("--output", required=True,
                        help="Output PDF path")
    args = parser.parse_args()

    print("Loading SR splice junctions …")
    sr_sj = build_sr_sj_set([Path(f) for f in args.sj])
    print(f"  SR-supported unique junctions: {len(sr_sj)}")

    print("Loading CAGE peaks …")
    cage_lookup = build_cage_lookup(Path(args.cage))

    # Build expression filter sets for tools that need it
    bambu_expr_ids = None
    if args.bambu_expression:
        print("\nLoading Bambu expression …")
        bambu_expr_ids = expressed_tx(pl.read_csv(args.bambu_expression))
        print(f"  Expressed Bambu transcripts: {len(bambu_expr_ids)}")

    isoseq_expr_ids = None
    if args.isoseq_classification and args.isoseq_expression:
        print("\nLoading IsoSeq expression …")
        cls_ids = pl.read_csv(
            args.isoseq_classification, separator="\t", null_values=["NA"]
        )["isoform"].to_list()
        isoseq_expr = (
            pl.read_parquet(args.isoseq_expression)
            .filter(pl.col("isoform").is_in(cls_ids))
        )
        isoseq_expr_ids = expressed_tx(isoseq_expr)
        print(f"  Expressed IsoSeq transcripts: {len(isoseq_expr_ids)}")

    polya = None
    if args.polya:
        print("Loading polyA clusters …")
        polya = pl.read_csv(
            args.polya, separator="\t", has_header=False,
            new_columns=["chrom", "start", "end", "name", "score", "strand"],
            schema_overrides={"chrom": pl.String},
        )
        print(f"  polyA clusters loaded: {len(polya)}")

    expr_ids_by_label = {
        "IsoSeq": isoseq_expr_ids,
        "Bambu": bambu_expr_ids,
    }

    labels, pct_sr_vals, pct_cage_vals, pct_polya_vals = [], [], [], []
    for label, gtf_path in args.tool:
        print(f"\nProcessing {label} ({gtf_path}) …")
        expr_ids = expr_ids_by_label.get(label)
        pct_sr, pct_cage, pct_polya = compute_support(
            Path(gtf_path), sr_sj, cage_lookup, expr_ids, polya
        )
        labels.append(label)
        pct_sr_vals.append(pct_sr)
        pct_cage_vals.append(pct_cage)
        if pct_polya is not None:
            pct_polya_vals.append(pct_polya)

    # ── Plot ──────────────────────────────────────────────────────────────────
    n = len(labels)
    x = np.arange(n)
    has_polya = bool(pct_polya_vals)
    width = 0.25 if has_polya else 0.35
    offsets = [-width, 0, width] if has_polya else [-width / 2, width / 2]

    fig, ax = plt.subplots(figsize=(12 if has_polya else 10, 6))

    def add_bars(vals, offset, label, color):
        bars = ax.bar(x + offset, vals, width, label=label, color=color)
        for bar in bars:
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.8,
                    f"{bar.get_height():.1f}%", ha="center", va="bottom", fontsize=9)

    add_bars(pct_sr_vals,   offsets[0], "SR junction support", "#2166ac")
    add_bars(pct_cage_vals, offsets[1], "CAGE support",        "#d6604d")
    if has_polya:
        add_bars(pct_polya_vals, offsets[2], "PolyA cluster support", "#4dac26")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("% of transcripts", fontsize=11)
    ax.set_ylim(0, 110)
    ax.legend(frameon=False, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.savefig(args.output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved plot to {args.output}")


if __name__ == "__main__":
    main()
