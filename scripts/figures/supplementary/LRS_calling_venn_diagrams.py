"""
Three-way Venn diagram comparison of transcript models from Bambu, IsoQuant, and Isoseq.

Transcripts are compared by intron chain: two transcripts are considered equivalent
if they share the same set of internal splice junctions (defined by the positions of
all internal exons, excluding the terminal 5' and 3' exons). Mono-exonic and
two-exon transcripts are excluded.

Transcripts are further filtered by expression: only transcripts with >5 reads in
>2 samples are included in the comparison.

Novel/known classification:
  - Bambu:    novelTranscript flag from supportedTxClassification.txt
  - IsoQuant: transcript IDs starting with "transcript" (novel) vs "ENST" (known)
  - Isoseq:   structural_category != "full_splice-match" (novel) vs == (known)

Output: venn_novel_known.pdf — side-by-side Venn diagrams for novel and known transcripts.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import polars as pl
from matplotlib_venn import venn3


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def read_gtf(file, attributes=["transcript_id"], keep_attributes=True):
    df = pl.read_csv(
        file,
        separator="\t",
        comment_prefix="#",
        schema_overrides={"seqname": pl.String},
        has_header=False,
        new_columns=["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"],
    ).with_columns(
        [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
    )
    if not keep_attributes:
        df = df.drop("attributes")
    return df


def build_exon_group(gtf):
    """Extract internal exons per transcript and group by transcript_id.

    Removes mono-exonic and two-exon transcripts, then strips the terminal
    exons (min start, max end), leaving only internal exons. Returns one row
    per transcript with lists of internal exon start/end positions.
    """
    return (
        gtf
        .filter(pl.col("feature") == "exon")
        .filter(pl.col("start").count().over("transcript_id") != 2)
        .filter(
            pl.col("start") != pl.col("start").min().over("transcript_id"),
            pl.col("end") != pl.col("end").max().over("transcript_id"),
        )
        .group_by(pl.col("transcript_id"))
        .agg(
            pl.col("seqname").first(),
            pl.col("start"),
            pl.col("end"),
            pl.col("strand").first(),
        )
    )


def expressed_tx(expr_df, min_reads=5, min_samples=2):
    """Return set of transcript IDs expressed >min_reads in >min_samples samples.

    Assumes the first column is the transcript ID and all remaining columns are
    per-sample read counts.
    """
    id_col = expr_df.columns[0]
    sample_cols = expr_df.columns[1:]
    return set(
        expr_df.filter(
            pl.sum_horizontal([pl.col(c) > min_reads for c in sample_cols]) > min_samples
        )[id_col].to_list()
    )


def get_chain_set(exon_group_df):
    """Convert an exon group dataframe to a set of intron chain keys.

    Each key is a tuple (seqname, strand, sorted_internal_exon_pairs) where
    sorted_internal_exon_pairs is a sorted tuple of (start, end) for all
    internal exons. This uniquely identifies the intron chain regardless of
    TSS/TTS variation.
    """
    keys = set()
    for row in exon_group_df.iter_rows(named=True):
        exons = tuple(sorted(zip(row["start"], row["end"])))
        keys.add((row["seqname"], row["strand"], exons))
    return keys


# ---------------------------------------------------------------------------
# Load data: Bambu
# ---------------------------------------------------------------------------

bambu_gtf        = read_gtf("nextflow_results/comapre_other_LRS_tools/bambu/supportedTranscriptModels.gtf")
bambu_exon_group = build_exon_group(bambu_gtf)
novel_bambu_tx   = (
    pl.read_csv("nextflow_results/comapre_other_LRS_tools/bambu/supportedTxClassification.txt", separator="\t")
    .filter(pl.col("novelTranscript") == True)["TXNAME"]
    .to_list()
)
bambu_novel      = bambu_exon_group.filter(pl.col("transcript_id").is_in(novel_bambu_tx))
bambu_known      = bambu_exon_group.filter(pl.col("transcript_id").is_in(novel_bambu_tx).not_())
bambu_expression = pl.read_csv("nextflow_results/comapre_other_LRS_tools/bambu/bambu_expression.csv")

# ---------------------------------------------------------------------------
# Load data: IsoQuant
# ---------------------------------------------------------------------------

# isoquant_gtf        = read_gtf("nextflow_results/comapre_other_LRS_tools/isoquant/isoquant_out/OUT/OUT.transcript_models.gtf")
isoquant_gtf = read_gtf("nextflow_results/comapre_other_LRS_tools/isoquant/isoquant_out/OUT/OUT.extended_annotation.gtf")
isoquant_exon_group = build_exon_group(isoquant_gtf)
isoquant_novel      = isoquant_exon_group.filter(pl.col("transcript_id").str.starts_with("transcript"))
isoquant_known      = isoquant_exon_group.filter(pl.col("transcript_id").str.starts_with("ENST"))
isoquant_expression_discovered = pl.read_csv(
    "nextflow_results/comapre_other_LRS_tools/isoquant/isoquant_out/OUT/OUT.discovered_transcript_grouped_file_name_counts.tsv",
    separator="\t",
)
isoquant_expression_known = pl.read_csv(
    "nextflow_results/comapre_other_LRS_tools/isoquant/isoquant_out/OUT/OUT.transcript_grouped_file_name_counts.tsv",
    separator="\t",
)

# ---------------------------------------------------------------------------
# Load data: Isoseq (SFARI)
# ---------------------------------------------------------------------------

SFARI_gtf            = read_gtf("nextflow_results/classify_and_count/final_transcripts.gtf")
SFARI_exon_group     = build_exon_group(SFARI_gtf)
SFARI_classification = pl.read_parquet("nextflow_results/classify_and_count/final_classification.parquet")
novel_SFARI_tx       = SFARI_classification.filter(pl.col("structural_category") != "full-splice_match")["isoform"].to_list()
SFARI_novel          = SFARI_exon_group.filter(pl.col("transcript_id").is_in(novel_SFARI_tx))
SFARI_known          = SFARI_exon_group.filter(pl.col("transcript_id").is_in(novel_SFARI_tx).not_())
SFARI_expression     = pl.read_parquet("nextflow_results/classify_and_count/final_expression.parquet")

# ---------------------------------------------------------------------------
# Expression filtering and intron chain sets
# ---------------------------------------------------------------------------

bambu_expr_tx    = expressed_tx(bambu_expression)
isoquant_expr_tx = expressed_tx(isoquant_expression_discovered) | expressed_tx(isoquant_expression_known)
isoseq_expr_tx   = expressed_tx(SFARI_expression)

novel_bambu_chains    = get_chain_set(bambu_novel.filter(pl.col("transcript_id").is_in(bambu_expr_tx)))
novel_isoquant_chains = get_chain_set(isoquant_novel.filter(pl.col("transcript_id").is_in(isoquant_expr_tx)))
novel_isoseq_chains   = get_chain_set(SFARI_novel.filter(pl.col("transcript_id").is_in(isoseq_expr_tx)))

known_bambu_chains    = get_chain_set(bambu_known.filter(pl.col("transcript_id").is_in(bambu_expr_tx)))
known_isoquant_chains = get_chain_set(isoquant_known.filter(pl.col("transcript_id").is_in(isoquant_expr_tx)))
known_isoseq_chains   = get_chain_set(SFARI_known.filter(pl.col("transcript_id").is_in(isoseq_expr_tx)))

labels = ("Bambu", "IsoQuant", "Isoseq")

# ---------------------------------------------------------------------------
# Print set membership percentages
# ---------------------------------------------------------------------------

def print_venn_percentages(A, B, C, labels):
    only_A   = A - B - C
    only_B   = B - A - C
    only_C   = C - A - B
    AB_only  = (A & B) - C
    AC_only  = (A & C) - B
    BC_only  = (B & C) - A
    ABC      = A & B & C
    total    = len(A | B | C)

    la, lb, lc = labels
    print(f"  Total unique transcripts: {total}")
    for label, s in [
        (f"{la} only",             only_A),
        (f"{lb} only",             only_B),
        (f"{lc} only",             only_C),
        (f"{la} & {lb} only",      AB_only),
        (f"{la} & {lc} only",      AC_only),
        (f"{lb} & {lc} only",      BC_only),
        (f"{la} & {lb} & {lc}",    ABC),
    ]:
        print(f"  {label}: {len(s)} ({100 * len(s) / total:.1f}%)")

print("Novel transcripts:")
print_venn_percentages(novel_bambu_chains, novel_isoquant_chains, novel_isoseq_chains, labels)
print("\nKnown transcripts:")
print_venn_percentages(known_bambu_chains, known_isoquant_chains, known_isoseq_chains, labels)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

venn3([novel_bambu_chains, novel_isoquant_chains, novel_isoseq_chains],
      set_labels=labels, ax=ax1)
ax1.set_title("Novel transcripts")
ax1.text(0.0, 1.0, "A", transform=ax1.title.get_transform(), fontsize=14, fontweight="bold", va=ax1.title.get_verticalalignment(), ha="left", clip_on=False)

venn3([known_bambu_chains, known_isoquant_chains, known_isoseq_chains],
      set_labels=labels, ax=ax2)
ax2.set_title("Known transcripts")
ax2.text(0.0, 1.0, "B", transform=ax2.title.get_transform(), fontsize=14, fontweight="bold", va=ax2.title.get_verticalalignment(), ha="left", clip_on=False)

fig.tight_layout()
fig.savefig("figures/supplementary/LRS_calling_venn_diagrams.pdf")
fig.savefig("figures/supplementary/LRS_calling_venn_diagrams.png", dpi=300)
plt.close(fig)
print("Saved LRS_calling_venn_diagrams.pdf and LRS_calling_venn_diagrams.png")
