#!/usr/bin/env python3
"""
Two-way Venn diagram comparison of transcript models from Bambu and IsoSeq.

Transcripts are compared by intron chain: two transcripts are considered equivalent
if they share the same set of internal splice junctions (defined by the positions of
all internal exons, excluding the terminal 5' and 3' exons). Mono-exonic and
two-exon transcripts are excluded.

Transcripts are further filtered by expression: only transcripts with >5 reads in
>2 samples are included in the comparison.

Novel/known classification:
  - Bambu:  novelTranscript flag from supportedTxClassification.txt
  - IsoSeq: structural_category != "full-splice_match" (novel) vs == (known)

Output: lrs_calling_venn.pdf — side-by-side Venn diagrams for novel and known transcripts.
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import polars as pl
from matplotlib_venn import venn2

sys.path.insert(0, str(Path(__file__).parent / "src"))
from utils import read_gtf


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def build_exon_group(gtf: pl.DataFrame) -> pl.DataFrame:
    """Extract internal exons per transcript, grouped by transcript_id.

    Removes mono-exonic and two-exon transcripts, then strips the terminal
    exons (min start, max end), leaving only internal exons.
    """
    return (
        gtf
        .filter(pl.col("feature") == "exon")
        .filter(pl.col("start").count().over("transcript_id") != 2)
        .filter(
            pl.col("start") != pl.col("start").min().over("transcript_id"),
            pl.col("end") != pl.col("end").max().over("transcript_id"),
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
    """Return set of transcript IDs expressed >min_reads in >min_samples samples."""
    id_col = expr_df.columns[0]
    sample_cols = expr_df.columns[1:]
    return set(
        expr_df.filter(
            pl.sum_horizontal([pl.col(c) > min_reads for c in sample_cols]) > min_samples
        )[id_col].to_list()
    )


def get_chain_set(exon_group_df: pl.DataFrame) -> set:
    """Convert an exon group dataframe to a set of intron chain keys.

    Each key is a tuple (seqname, strand, sorted_internal_exon_pairs).
    """
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


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)

    # Bambu inputs
    parser.add_argument("--bambu-gtf", required=True,
                        help="Bambu supportedTranscriptModels.gtf")
    parser.add_argument("--bambu-classification", required=True,
                        help="Bambu supportedTxClassification.txt")
    parser.add_argument("--bambu-expression", required=True,
                        help="Bambu bambu_expression.csv")

    # IsoSeq inputs
    parser.add_argument("--isoseq-gtf", required=True,
                        help="IsoSeq final_transcripts.gtf")
    parser.add_argument("--isoseq-classification", required=True,
                        help="IsoSeq final_classification.parquet")
    parser.add_argument("--isoseq-expression", required=True,
                        help="IsoSeq final_expression.parquet")

    parser.add_argument("--output", required=True,
                        help="Output PDF path")

    args = parser.parse_args()

    # ── Bambu ─────────────────────────────────────────────────────────────────
    print("Loading Bambu …")
    bambu_gtf        = read_gtf(args.bambu_gtf, attributes=["transcript_id"])
    bambu_exon_group = build_exon_group(bambu_gtf)
    novel_bambu_tx   = (
        pl.read_csv(args.bambu_classification, separator="\t")
        .filter(pl.col("novelTranscript") == True)["TXNAME"]
        .to_list()
    )
    bambu_novel      = bambu_exon_group.filter(pl.col("transcript_id").is_in(novel_bambu_tx))
    bambu_known      = bambu_exon_group.filter(~pl.col("transcript_id").is_in(novel_bambu_tx))
    bambu_expression = pl.read_csv(args.bambu_expression)

    # ── IsoSeq ────────────────────────────────────────────────────────────────
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

    # ── Expression filtering ───────────────────────────────────────────────────
    print("Filtering by expression …")
    bambu_expr_tx  = expressed_tx(bambu_expression)
    isoseq_expr_tx = expressed_tx(SFARI_expression)

    novel_bambu_chains  = get_chain_set(bambu_novel.filter(pl.col("transcript_id").is_in(bambu_expr_tx)))
    novel_isoseq_chains = get_chain_set(SFARI_novel.filter(pl.col("transcript_id").is_in(isoseq_expr_tx)))

    known_bambu_chains  = get_chain_set(bambu_known.filter(pl.col("transcript_id").is_in(bambu_expr_tx)))
    known_isoseq_chains = get_chain_set(SFARI_known.filter(pl.col("transcript_id").is_in(isoseq_expr_tx)))

    labels = ("Bambu", "IsoSeq")

    print("\nNovel transcripts:")
    print_venn_percentages(novel_bambu_chains, novel_isoseq_chains, labels)
    print("\nKnown transcripts:")
    print_venn_percentages(known_bambu_chains, known_isoseq_chains, labels)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    venn2([novel_bambu_chains, novel_isoseq_chains], set_labels=labels, ax=ax1)
    ax1.set_title("Novel transcripts")
    ax1.text(0.0, 1.0, "A", transform=ax1.title.get_transform(),
             fontsize=14, fontweight="bold",
             va=ax1.title.get_verticalalignment(), ha="left", clip_on=False)

    venn2([known_bambu_chains, known_isoseq_chains], set_labels=labels, ax=ax2)
    ax2.set_title("Known transcripts")
    ax2.text(0.0, 1.0, "B", transform=ax2.title.get_transform(),
             fontsize=14, fontweight="bold",
             va=ax2.title.get_verticalalignment(), ha="left", clip_on=False)

    fig.tight_layout()
    out = Path(args.output)
    fig.savefig(out, bbox_inches="tight")
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved {out} and {out.with_suffix('.png')}")


if __name__ == "__main__":
    main()
