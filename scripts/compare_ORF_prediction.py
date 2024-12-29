import polars as pl
from src.utils import read_gtf
from src.single_cell import SingleCell
import seaborn as sns
import matplotlib.pyplot as plt
import seaborn as sns

# Load transcript annotations
lr_bulk = SingleCell("results/long_read/pbid_filtered.h5ad")
pbid_to_transcript = lr_bulk.var\
    .filter(
        pl.col("associated_transcript")!="novel"
    )\
    .filter(
        pl.col("structural_category2") == "full-splice_match"
    )["pbid", "associated_transcript"]

# Load GTFs
GENCODE = read_gtf("/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.annotation.gtf",keep_attributes=False,attributes = ["transcript_id", "transcript_type"])\
    .filter(
        pl.col("feature") == "CDS"
    )\
    .sort(["transcript_id", "start"])\
    .group_by("transcript_id")\
    .agg(
        pl.col("transcript_type").unique().map_elements(lambda x: x[0], return_dtype=pl.String),
        orf_start = pl.col("start").first(),
        orf_end = pl.col("end").last()
    )

LRP_CDS = read_gtf("nextflow_results/LRP_CDS.gtf", keep_attributes=False)\
    .with_columns(
        pl.col("transcript_id").str.split("|").map_elements(lambda x: x[1], return_dtype=pl.String)
    )\
    .filter(
        pl.col("feature") == "CDS"
    )\
    .rename({"transcript_id": "pbid"})\
    .join(
        pbid_to_transcript, on = "pbid", how = "left"
    )\
    .drop_nulls()\
    .rename({"associated_transcript": "transcript_id"})\
    .sort(["transcript_id", "start"])\
    .group_by("transcript_id")\
    .agg(
        orf_start = pl.col("start").first(),
        orf_end = pl.col("end").last()
    )\
    .with_columns(
        pl.col("transcript_id").cast(pl.String)
    )

TransDecoder_CDS = read_gtf("full_nt.fasta.transdecoder.genome_updated.gff3", keep_attributes=False)\
    .filter(
        pl.col("feature") == "CDS"
    )\
    .rename({"transcript_id": "pbid"})\
    .join(
        pbid_to_transcript, on = "pbid", how = "left"
    )\
    .drop_nulls()\
    .rename({"associated_transcript": "transcript_id"})\
    .sort(["transcript_id", "start"])\
    .group_by("transcript_id")\
    .agg(
        orf_start = pl.col("start").first(),
        orf_end = pl.col("end").last()
    )\
    .with_columns(
        pl.col("transcript_id").cast(pl.String)
    )

# Get diff between predicted and annotated ORF starts and ends
LRP_diff = LRP_CDS\
    .join(GENCODE, on="transcript_id", how="left")\
    .drop_nulls()\
    .with_columns(
        diff_start = (pl.col("orf_start_right") - pl.col("orf_start")).abs(),
        diff_end = (pl.col("orf_end_right") - pl.col("orf_end")).abs()
    )

TD_diff = TransDecoder_CDS\
    .join(GENCODE, on="transcript_id", how="left")\
    .drop_nulls()\
    .with_columns(
        diff_start = (pl.col("orf_start_right") - pl.col("orf_start")).abs(),
        diff_end = (pl.col("orf_end_right") - pl.col("orf_end")).abs()
    )

LRP_start_pie = LRP_diff["diff_start"].hist([0, 10, 100, 1000]).with_columns(percentage = pl.col("count") / LRP_diff.shape[0])
LRP_end_pie = LRP_diff["diff_end"].hist([0, 10, 100, 1000]).with_columns(percentage = pl.col("count") / LRP_diff.shape[0])
TD_start_pie = TD_diff["diff_start"].hist([0, 10, 100, 1000]).with_columns(percentage = pl.col("count") / TD_diff.shape[0])
TD_end_pie = TD_diff["diff_end"].hist([0, 10, 100, 1000]).with_columns(percentage = pl.col("count") / TD_diff.shape[0])

# Plot TransDecoder prediction
fig, axs = plt.subplots(2, 2, figsize=(10, 10))
axs[0, 0].hist(TD_diff["diff_start"], range=(0, 5), label = "diff_start")
axs[0, 1].hist(TD_diff["diff_end"], range=(0, 5), label = "diff_end")
axs[0, 0].set_title("Distance between predicted and \nannotated ORF start coordinate (bp)")
axs[0, 1].set_title("Distance between predicted and \nannotated ORF end coordinate (bp)")
axs[1, 0].pie(TD_start_pie["percentage"], labels=TD_start_pie["category"].cast(pl.String))
axs[1, 1].pie(TD_end_pie["percentage"], labels=TD_end_pie["category"].cast(pl.String))
fig.suptitle("TransDecoder prediction", fontsize=20, fontweight='bold')
plt.savefig("figures/TD_prediction.png")

# Plot CPAT prediction
fig, axs = plt.subplots(2, 2, figsize=(10, 10))
axs[0, 0].hist(LRP_diff["diff_start"], range=(0, 5), label = "diff_start")
axs[0, 1].hist(LRP_diff["diff_end"], range=(0, 5), label = "diff_end")
axs[0, 0].set_title("Distance between predicted and \nannotated ORF start coordinate (bp)")
axs[0, 1].set_title("Distance between predicted and \nannotated ORF end coordinate (bp)")
axs[1, 0].pie(LRP_start_pie["percentage"], labels=LRP_start_pie["category"].cast(pl.String))
axs[1, 1].pie(LRP_end_pie["percentage"], labels=LRP_end_pie["category"].cast(pl.String))
fig.suptitle("CPAT prediction", fontsize=20, fontweight='bold')
plt.savefig("figures/CPAT_prediction.png")


read_gtf("nextflow_results/LRP_CDS.gtf", keep_attributes=False)\
    .with_columns(
        pl.col("transcript_id").str.split("|").map_elements(lambda x: x[1], return_dtype=pl.String)
    )\
    .with_columns(
        gene_id = pl.col("transcript_id").str.split(".").map_elements(lambda x: "".join([x[0], ".", x[1]]), return_dtype=pl.String)
    )\
    .with_columns(
        attributes = pl.lit('gene_id "') + pl.col("gene_id") + pl.lit('"; transcript_id "') + pl.col("transcript_id") + pl.lit('";')
    ).drop(["transcript_id", "gene_id"])\
    .write_csv("proc/LRP_CDS_reformatted.gtf", separator="\t", quote_style="never", include_header=False)