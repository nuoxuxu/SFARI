#!/usr/bin/env python3
"""
Combined LRS transcript calling comparison figure (5 panels).

Row 1 — Venn diagrams comparing Bambu and IsoSeq transcript models by intron chain.
  A: Novel transcripts   B: Known transcripts
Row 2 — Multi-omic validation of expressed long-read transcript models.
  C: SR splice junction support   D: CAGE support   E: PolyA cluster support

Transcripts are compared by intron chain (internal splice junctions only;
mono-exonic and two-exon transcripts excluded).
Expression filter: >5 reads in >2 samples.
Novel/known: Bambu uses novelTranscript flag; IsoSeq uses structural_category.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import polars as pl
from matplotlib_venn import venn2

sys.path.insert(0, str(Path(__file__).parent / "src"))
from utils import read_gtf, read_SJ

FONT_PANEL   = 18
FONT_TITLE   = 16
FONT_AXIS    = 14
FONT_TICK    = 13
FONT_BAR_VAL = 12
FONT_VENN    = 13


# ── Venn helpers ───────────────────────────────────────────────────────────────

def build_exon_group(gtf: pl.DataFrame) -> pl.DataFrame:
    return (
        gtf
        .filter(pl.col("feature") == "exon")
        .filter(pl.col("start").count().over("transcript_id") != 2)
        .filter(
            pl.col("start") != pl.col("start").min().over("transcript_id"),
            pl.col("end")   != pl.col("end").max().over("transcript_id"),
        )
        .group_by("transcript_id")
        .agg(
            pl.col("seqname").first(),
            pl.col("start"),
            pl.col("end"),
            pl.col("strand").first(),
        )
    )


def expressed_tx(expr_df: pl.DataFrame, min_reads: int = 5, min_samples: int = 2) -> set:
    id_col = expr_df.columns[0]
    sample_cols = expr_df.columns[1:]
    return set(
        expr_df.filter(
            pl.sum_horizontal([pl.col(c) > min_reads for c in sample_cols]) > min_samples
        )[id_col].to_list()
    )


def get_chain_set(exon_group_df: pl.DataFrame) -> set:
    keys = set()
    for row in exon_group_df.iter_rows(named=True):
        exons = tuple(sorted(zip(row["start"], row["end"])))
        keys.add((row["seqname"], row["strand"], exons))
    return keys


def print_venn_percentages(A: set, B: set, labels: tuple) -> None:
    only_A = A - B
    only_B = B - A
    AB = A & B
    total = len(A | B)
    la, lb = labels
    print(f"  Total unique transcripts: {total}")
    for label, s in [
        (f"{la} only", only_A),
        (f"{lb} only", only_B),
        (f"{la} & {lb}", AB),
    ]:
        print(f"  {label}: {len(s)} ({100 * len(s) / total:.1f}%)")


def draw_venn(ax, set_a: set, set_b: set, labels: tuple, title: str) -> None:
    v = venn2([set_a, set_b], set_labels=labels, ax=ax)
    for text in (v.set_labels or []):
        if text is not None:
            text.set_fontsize(FONT_VENN)
    for text in (v.subset_labels or []):
        if text is not None:
            text.set_fontsize(FONT_VENN)
    ax.set_title(title, fontsize=FONT_TITLE, pad=10)


# ── SR / CAGE / polyA helpers ──────────────────────────────────────────────────

def build_sr_sj_set(sj_files: list) -> pl.DataFrame:
    return (
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
        ends   = grp["end"].to_numpy()
        order  = np.argsort(starts)
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
        .select([pl.col("transcript_id"), pl.col("seqname"), pl.col("pos")])
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
    polya: pl.DataFrame,
    expr_ids: set | None = None,
) -> tuple[float, float, float]:
    gtf = read_gtf(gtf_path, attributes=["transcript_id"])
    exons       = gtf.filter(pl.col("feature") == "exon")
    transcripts = gtf.filter(pl.col("feature") == "transcript")
    all_ids = transcripts["transcript_id"].unique()
    total_raw = len(all_ids)
    if expr_ids is not None:
        expr_list   = list(expr_ids)
        transcripts = transcripts.filter(pl.col("transcript_id").is_in(expr_list))
        exons       = exons.filter(pl.col("transcript_id").is_in(expr_list))
        all_ids     = transcripts["transcript_id"].unique()
        print(f"  Transcripts: {total_raw} → {len(all_ids)} (after expression filter)")
    else:
        print(f"  Transcripts: {total_raw}")
    total = len(all_ids)

    # SR splice junction support
    lr_junctions = (
        exons
        .sort(["transcript_id", "start"])
        .with_columns(pl.col("start").shift(-1).over("transcript_id").alias("next_exon_start"))
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
    n_sr   = total - len(all_ids.filter(all_ids.is_in(unsupported_ids["transcript_id"].to_list())))
    pct_sr = n_sr / total * 100
    print(f"  SR-supported: {n_sr} ({pct_sr:.1f}%)")

    # CAGE support
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
    n_cage   = sum(cage_flags)
    pct_cage = n_cage / total * 100
    print(f"  CAGE-supported: {n_cage} ({pct_cage:.1f}%)")

    # PolyA cluster support
    n_polya, pct_polya = compute_polya_support(transcripts, polya)
    print(f"  PolyA-supported: {n_polya} ({pct_polya:.1f}%)")

    return pct_sr, pct_cage, pct_polya


def draw_bar(ax, labels: list, vals: list, color: str, title: str) -> None:
    x = np.arange(len(labels))
    width = 0.5
    bars = ax.bar(x, vals, width, color=color)
    for bar in bars:
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.8,
            f"{bar.get_height():.1f}%",
            ha="center", va="bottom", fontsize=FONT_BAR_VAL,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=FONT_TICK)
    ax.tick_params(axis="y", labelsize=FONT_TICK)
    ax.set_ylabel("% of transcripts", fontsize=FONT_AXIS)
    ax.set_ylim(0, 110)
    ax.set_title(title, fontsize=FONT_TITLE, pad=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__)

    # Bambu
    p.add_argument("--bambu-gtf",            required=True, help="supportedTranscriptModels.gtf")
    p.add_argument("--bambu-classification", required=True, help="supportedTxClassification.txt")
    p.add_argument("--bambu-expression",     required=True, help="bambu_expression.csv")

    # IsoSeq
    p.add_argument("--isoseq-gtf",            required=True, help="final_transcripts.gtf")
    p.add_argument("--isoseq-classification", required=True, help="final_classification.parquet")
    p.add_argument("--isoseq-expression",     required=True, help="final_expression.parquet")
    p.add_argument("--isoseq-filtered-lite-classification", required=True,
                   help="filtered_lite_classification.txt (for SR/CAGE/polyA expression filter)")

    # SR / CAGE / polyA
    p.add_argument("--sj",   nargs="+", required=True, help="STAR SJ.out.tab files")
    p.add_argument("--cage", required=True,             help="CAGE BED file (refTSS)")
    p.add_argument("--polya", required=True,            help="PolyA cluster BED file")

    p.add_argument("--output", required=True, help="Output PDF path")

    args = p.parse_args()

    # ── Venn data ──────────────────────────────────────────────────────────────
    print("Loading Bambu …")
    bambu_gtf        = read_gtf(args.bambu_gtf, attributes=["transcript_id"])
    bambu_exon_group = build_exon_group(bambu_gtf)
    novel_bambu_tx   = (
        pl.read_csv(args.bambu_classification, separator="\t")
        .filter(pl.col("novelTranscript") == True)["TXNAME"]
        .to_list()
    )
    bambu_novel   = bambu_exon_group.filter(pl.col("transcript_id").is_in(novel_bambu_tx))
    bambu_known   = bambu_exon_group.filter(~pl.col("transcript_id").is_in(novel_bambu_tx))
    bambu_expr_tx = expressed_tx(pl.read_csv(args.bambu_expression))

    print("Loading IsoSeq …")
    SFARI_gtf            = read_gtf(args.isoseq_gtf, attributes=["transcript_id"])
    SFARI_exon_group     = build_exon_group(SFARI_gtf)
    SFARI_classification = pl.read_parquet(args.isoseq_classification)
    novel_SFARI_tx       = (
        SFARI_classification
        .filter(pl.col("structural_category") != "full-splice_match")["isoform"]
        .to_list()
    )
    SFARI_novel      = SFARI_exon_group.filter(pl.col("transcript_id").is_in(novel_SFARI_tx))
    SFARI_known      = SFARI_exon_group.filter(~pl.col("transcript_id").is_in(novel_SFARI_tx))
    SFARI_expression = pl.read_parquet(args.isoseq_expression)
    isoseq_expr_tx   = expressed_tx(SFARI_expression)

    novel_bambu_chains  = get_chain_set(bambu_novel.filter(pl.col("transcript_id").is_in(bambu_expr_tx)))
    novel_isoseq_chains = get_chain_set(SFARI_novel.filter(pl.col("transcript_id").is_in(isoseq_expr_tx)))
    known_bambu_chains  = get_chain_set(bambu_known.filter(pl.col("transcript_id").is_in(bambu_expr_tx)))
    known_isoseq_chains = get_chain_set(SFARI_known.filter(pl.col("transcript_id").is_in(isoseq_expr_tx)))

    venn_labels = ("Bambu", "IsoSeq")
    print("\nNovel transcripts:")
    print_venn_percentages(novel_bambu_chains, novel_isoseq_chains, venn_labels)
    print("\nKnown transcripts:")
    print_venn_percentages(known_bambu_chains, known_isoseq_chains, venn_labels)

    # ── SR / CAGE / polyA data ─────────────────────────────────────────────────
    print("\nLoading SR splice junctions …")
    sr_sj = build_sr_sj_set([Path(f) for f in args.sj])
    print(f"  Unique SR junctions: {len(sr_sj)}")

    print("Loading CAGE peaks …")
    cage_lookup = build_cage_lookup(Path(args.cage))

    print("Loading polyA clusters …")
    polya = pl.read_csv(
        args.polya, separator="\t", has_header=False,
        new_columns=["chrom", "start", "end", "name", "score", "strand"],
        schema_overrides={"chrom": pl.String},
    )
    print(f"  PolyA clusters loaded: {len(polya)}")

    cls_ids = pl.read_csv(
        args.isoseq_filtered_lite_classification, separator="\t", null_values=["NA"]
    )["isoform"].to_list()
    isoseq_cage_expr = expressed_tx(
        pl.read_parquet(args.isoseq_expression).filter(pl.col("isoform").is_in(cls_ids))
    )

    tool_specs = [
        ("IsoSeq", args.isoseq_gtf, isoseq_cage_expr),
        ("Bambu",  args.bambu_gtf,  bambu_expr_tx),
    ]
    tool_labels, pct_sr_vals, pct_cage_vals, pct_polya_vals = [], [], [], []
    for label, gtf_path, expr_ids in tool_specs:
        print(f"\nProcessing {label} …")
        pct_sr, pct_cage, pct_polya = compute_support(
            Path(gtf_path), sr_sj, cage_lookup, polya, expr_ids
        )
        tool_labels.append(label)
        pct_sr_vals.append(pct_sr)
        pct_cage_vals.append(pct_cage)
        pct_polya_vals.append(pct_polya)

    # ── Plot ───────────────────────────────────────────────────────────────────
    # Row 1: 2 Venns (A, B) each half-width
    # Row 2: 3 bar charts (C, D, E) each one-third width
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 6, figure=fig, hspace=0.45, wspace=0.5)

    ax_a = fig.add_subplot(gs[0, 0:3])
    ax_b = fig.add_subplot(gs[0, 3:6])
    ax_c = fig.add_subplot(gs[1, 0:2])
    ax_d = fig.add_subplot(gs[1, 2:4])
    ax_e = fig.add_subplot(gs[1, 4:6])

    draw_venn(ax_a, novel_bambu_chains, novel_isoseq_chains, venn_labels, "Novel transcripts")
    draw_venn(ax_b, known_bambu_chains, known_isoseq_chains, venn_labels, "Known transcripts")
    draw_bar(ax_c, tool_labels, pct_sr_vals,    "#2166ac", "SR splice junction support")
    draw_bar(ax_d, tool_labels, pct_cage_vals,  "#d6604d", "CAGE support")
    draw_bar(ax_e, tool_labels, pct_polya_vals, "#4dac26", "PolyA cluster support")

    for ax, letter in [(ax_a, "A"), (ax_b, "B"), (ax_c, "C"), (ax_d, "D"), (ax_e, "E")]:
        ax.text(-0.08, 1.08, letter, transform=ax.transAxes,
                fontsize=FONT_PANEL, fontweight="bold", va="top", ha="left", clip_on=False)

    out = Path(args.output)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {out} and {out.with_suffix('.png')}")


if __name__ == "__main__":
    main()
