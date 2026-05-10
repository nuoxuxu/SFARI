"""
Compute the percentage of fusion transcripts whose splice junctions are all
supported by short-read RNA-seq, and plot as a pie chart.
"""

import sys
from pathlib import Path

import polars as pl
from ryp import r, to_r

sys.path.insert(0, "bin")
from src.utils import read_gtf, read_SJ

# ── Paths ──────────────────────────────────────────────────────────────────────
CLASSIFICATION = Path("nextflow_results/classify_and_count/final_classification.parquet")
GTF            = Path("nextflow_results/classify_and_count/final_transcripts.gtf")
STAR_DIR       = Path("nextflow_results/STAR")

# ── Step 1: Load classification; filter to fusion isoforms ────────────────────
classification = pl.read_parquet(CLASSIFICATION)
fusion = classification.filter(pl.col("structural_category") == "fusion")
print(f"Fusion isoforms: {len(fusion)}")

# ── Step 2: Build SR splice junction support set ───────────────────────────────
sj_files = list(STAR_DIR.glob("*.SJ.out.tab"))
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

# ── Step 3: Extract LR splice junctions from GTF for fusion isoforms ──────────
gtf = read_gtf(GTF, attributes=["transcript_id"])
fusion_exons = (
    gtf.filter(pl.col("feature") == "exon")
    .filter(pl.col("transcript_id").is_in(fusion["isoform"].to_list()))
)

lr_junctions = (
    fusion_exons
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

# ── Step 4: Classify isoforms by SR support ────────────────────────────────────
unsupported_transcripts = (
    lr_junctions
    .join(sr_sj, on=["chrom", "start", "end", "strand"], how="anti")
    .select("transcript_id")
    .unique()
)

n_supported   = fusion.filter(~pl.col("isoform").is_in(unsupported_transcripts["transcript_id"])).height
n_unsupported = len(fusion) - n_supported
pct_supported = n_supported / len(fusion) * 100

print(f"SR-supported: {n_supported} / {len(fusion)} ({pct_supported:.1f}%)")

# ── Step 5: Plot ───────────────────────────────────────────────────────────────
pie_data = pl.DataFrame({
    "label": ["SR-supported", "Not SR-supported"],
    "n":     [n_supported, n_unsupported],
})
to_r(pie_data, "pie_data")

r("""
library(ggplot2)
library(dplyr)

pie_data <- pie_data %>%
    mutate(
        pct   = n / sum(n) * 100,
        label = factor(label, levels = c("SR-supported", "Not SR-supported"))
    )

ggplot(pie_data, aes(x = "", y = pct, fill = label)) +
    geom_col(width = 1, colour = "white") +
    coord_polar(theta = "y") +
    geom_text(
        aes(label = paste0(round(pct, 1), "%")),
        position = position_stack(vjust = 0.5),
        size = 5, colour = "white", fontface = "bold"
    ) +
    scale_fill_manual(
        values = c("SR-supported" = "#00BFC4", "Not SR-supported" = "#F8766D"),
        name   = NULL
    ) +
    labs(title = "Fusion isoforms: SR splice junction support") +
    theme_void(base_size = 14) +
    theme(
        plot.title    = element_text(hjust = 0.5, size = 14),
        legend.position = "bottom"
    )

ggsave("figures/supplementary/fusion_sr_support.pdf", width = 4, height = 4.5)
""")
