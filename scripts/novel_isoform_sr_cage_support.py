"""
Compute the percentage of novel isoforms (structural_category != 'full-splice_match')
whose splice junctions are all supported by short-read RNA-seq (STAR SJ.out.tab files),
and whose 5' ends are within 100 bp of a CAGE-seq peak (refTSS).
"""

import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).parent.parent / "bin"))
from utils import read_gtf, read_SJ

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
CLASSIFICATION = BASE / "nextflow_results/classify_and_count/final_classification.parquet"
GTF            = BASE / "nextflow_results/classify_and_count/final_transcripts.gtf"
STAR_DIR       = BASE / "nextflow_results/STAR"
CAGE_BED       = BASE / "data/refTSS_v3.3_human_coordinate.hg38.sorted.bed"

# ── Step 1: Load classification; define novel isoforms ─────────────────────────
print("Loading classification …")
classification = pl.read_parquet(CLASSIFICATION)

print(f"Total isoforms: {len(classification)}")
print("Structural category counts:")
print(classification["structural_category"].value_counts().sort("count", descending=True))

novel = classification.filter(pl.col("structural_category") != "full-splice_match")
print(f"\nNovel isoforms (non-FSM): {len(novel)}")

# ── Step 2: Build SR splice junction support set ───────────────────────────────
print("\nLoading STAR splice junctions …")
sj_files = list(STAR_DIR.glob("*.SJ.out.tab"))
print(f"  Found {len(sj_files)} SJ files")

sr_sj = (
    pl.concat([read_SJ(f) for f in sj_files])
    .filter((pl.col("unique_reads") + pl.col("multi_reads")) > 0)
    .with_columns(
        pl.col("strand").map_elements(
            lambda s: "+" if s == 1 else ("-" if s == 2 else "."),
            return_dtype=pl.String
        )
    )
    .select(["chrom", "start", "end", "strand"])
    .unique()
)
print(f"  SR-supported unique junctions: {len(sr_sj)}")

# ── Step 3: Extract LR splice junctions from GTF for novel isoforms ────────────
print("\nParsing GTF exons …")
gtf = read_gtf(GTF, attributes=["transcript_id"])
gtf_exons = gtf.filter(pl.col("feature") == "exon")
print(f"  Total exon rows: {len(gtf_exons)}")

novel_ids = novel["isoform"]

# Derive splice junctions only for novel isoforms (saves memory)
novel_exons = gtf_exons.filter(pl.col("transcript_id").is_in(novel_ids.to_list()))
print(f"  Novel isoform exon rows: {len(novel_exons)}")

# Derive splice junctions from exon pairs (pure polars, handles mono-exon correctly)
# For each transcript, sorted by start: intron_start = exon_i.end + 1,
# intron_end = exon_{i+1}.start - 1
lr_junctions = (
    novel_exons
    .sort(["transcript_id", "start"])
    .with_columns([
        pl.col("start").shift(-1).over("transcript_id").alias("next_exon_start"),
    ])
    .with_columns(
        (pl.col("next_exon_start") - 1).alias("sj_end"),   # exon_{i+1}.start - 1
        (pl.col("end") + 1).alias("sj_start"),             # exon_i.end + 1
    )
    # Drop the last exon of each transcript (no following exon → null next_exon_start)
    .filter(pl.col("next_exon_start").is_not_null())
    .select([
        pl.col("transcript_id"),
        pl.col("seqname").alias("chrom"),
        pl.col("sj_start").alias("start"),
        pl.col("sj_end").alias("end"),
        pl.col("strand"),
    ])
)
print(f"  Novel isoform splice junctions: {len(lr_junctions)}")

# ── Step 4: Find novel isoforms with all SJs supported ────────────────────────
# Anti-join: junctions in LR that have NO match in SR
unsupported = lr_junctions.join(
    sr_sj,
    on=["chrom", "start", "end", "strand"],
    how="anti"
)
print(f"\n  Unsupported LR junctions: {len(unsupported)}")

# Transcripts that have at least one unsupported junction
transcripts_with_unsupported = unsupported.select("transcript_id").unique()
print(f"  Novel isoforms with ≥1 unsupported junction: {len(transcripts_with_unsupported)}")

# Novel isoforms with ALL junctions supported (including mono-exon, which have no junctions)
novel_sr_supported = novel.filter(
    ~pl.col("isoform").is_in(transcripts_with_unsupported["transcript_id"])
)
pct_sr = len(novel_sr_supported) / len(novel) * 100

# Mono-exon count for reference
novel_multiexon = novel.filter(pl.col("isoform").is_in(lr_junctions["transcript_id"].unique()))
novel_monoexon = novel.filter(~pl.col("isoform").is_in(lr_junctions["transcript_id"].unique()))

# ── Step 5: CAGE-seq support from refTSS BED ──────────────────────────────────
print("\nComputing CAGE-seq support …")

# Get 5' end position for each novel isoform from GTF transcript features
gtf_transcripts = gtf.filter(pl.col("feature") == "transcript")
tss_df = (
    gtf_transcripts
    .filter(pl.col("transcript_id").is_in(novel_ids.to_list()))
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
print(f"  Novel isoform TSS positions: {len(tss_df)}")

# Load CAGE peaks (BED: chrom, start, end, name, score, strand — 0-based coords)
cage = (
    pl.read_csv(
        CAGE_BED, separator="\t", has_header=False,
        schema_overrides={"column_1": pl.String},
    )
    .select(
        pl.col("column_1").alias("chrom"),
        pl.col("column_2").alias("start"),
        pl.col("column_3").alias("end"),
        pl.col("column_6").alias("strand"),
    )
)
print(f"  CAGE peaks loaded: {len(cage)}")

# Build a lookup dict: (chrom, strand) → (sorted_starts, ends) as numpy arrays
cage_lookup: dict = {}
for (chrom, strand), grp in cage.group_by(["chrom", "strand"]):
    starts = grp["start"].to_numpy()
    ends   = grp["end"].to_numpy()
    order  = np.argsort(starts)
    cage_lookup[(chrom, strand)] = (starts[order], ends[order])

def within_100bp(tss_pos: int, starts: np.ndarray, ends: np.ndarray) -> bool:
    """True if tss_pos is within 100 bp of any peak [start, end] (0-based BED coords)."""
    # Binary search: first peak whose start > tss_pos + 100 — only need peaks before that
    idx = np.searchsorted(starts, tss_pos + 101)
    if idx == 0:
        return False
    # Candidates: peaks with start <= tss_pos + 100
    cand_ends = ends[:idx]
    # Among those, check if any end >= tss_pos - 100
    return bool(np.any(cand_ends >= tss_pos - 100))

# Evaluate CAGE support per isoform
cage_flags = []
for row in tss_df.iter_rows(named=True):
    key = (row["chrom"], row["strand"])
    if key not in cage_lookup:
        cage_flags.append(False)
    else:
        starts, ends = cage_lookup[key]
        cage_flags.append(within_100bp(row["tss_pos"], starts, ends))

tss_df = tss_df.with_columns(pl.Series("cage_support", cage_flags))
cage_supported_ids = tss_df.filter(pl.col("cage_support"))["transcript_id"]
print(f"  Novel isoforms with CAGE support: {len(cage_supported_ids)}")

novel_cage_supported = novel.filter(pl.col("isoform").is_in(cage_supported_ids.to_list()))
pct_cage = len(novel_cage_supported) / len(novel) * 100

# ── Step 6: Report ────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("RESULTS")
print("=" * 60)
print(f"Total isoforms:              {len(classification)}")
print(f"Novel isoforms (non-FSM):    {len(novel)}")
print(f"  Mono-exon novel:           {len(novel_monoexon)}")
print(f"  Multi-exon novel:          {len(novel_multiexon)}")
print()
print(f"SR splice junction support:")
print(f"  Supported novel isoforms:  {len(novel_sr_supported)} / {len(novel)}")
print(f"  Percentage:                {pct_sr:.1f}%")
print()
print(f"CAGE-seq support (|dist| ≤ 100 bp):")
print(f"  Supported novel isoforms:  {len(novel_cage_supported)} / {len(novel)}")
print(f"  Percentage:                {pct_cage:.1f}%")

# ── Per-category breakdown ────────────────────────────────────────────────────
print("\n" + "-" * 60)
print("Per structural_category breakdown:")
print("-" * 60)

for cat in novel["structural_category"].unique().sort().to_list():
    cat_df = novel.filter(pl.col("structural_category") == cat)
    n = len(cat_df)
    n_sr = len(novel_sr_supported.filter(pl.col("structural_category") == cat))
    n_cage = len(novel_cage_supported.filter(pl.col("structural_category") == cat))
    print(f"\n{cat} (n={n}):")
    print(f"  SR support:   {n_sr} ({n_sr/n*100:.1f}%)")
    print(f"  CAGE support: {n_cage} ({n_cage/n*100:.1f}%)")
