import polars as pl
from src.utils import read_gtf
import numpy as np
import polars.selectors as cs
import os

# Define functions

def read_refmap(file):
    out = pl.read_csv(file, separator="\t")\
    .filter(pl.col("class_code")=="=")\
    .with_columns(
        pl.col("qry_id_list").str.split(",").map_elements(lambda s: [e.split("|")[1] for e in s], return_dtype=pl.List(pl.String)).alias("qry_id_list"))\
    .explode("qry_id_list")\
    .select(["qry_id_list", "ref_id"])\
    .rename({"qry_id_list": "transcript_id"})
    return out

# Read input files

gencode_path = "".join([os.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf"])
gencode_gtf = read_gtf(gencode_path, ["gene_name", "transcript_id", "gene_id"])\
    .filter(pl.col("feature")=="exon")

TALON_gtf = read_gtf("nextflow_results/V47/compare/cp_vz_0.75_min_7_recovery_talon_hg38.gtf", ["gene_name", "transcript_id", "gene_id"])\
    .filter(
        pl.col("feature")=="exon",
        pl.col("gene_name").str.starts_with("TALON").not_()
    )

classification = pl.read_parquet("nextflow_results/V47/final_classification.parquet")

final_transcripts_gtf = read_gtf(
    "nextflow_results/V47/final_transcripts.gtf", 
    ["transcript_id"])\
    .rename({"transcript_id": "isoform"})\
    .filter(pl.col("feature")=="exon")\
    .join(classification["isoform", "associated_transcript", "associated_gene"], on="isoform", how="left")\
    .filter(pl.col("associated_gene").cast(pl.String).str.starts_with("novelGene").not_())\
    .with_columns(
        pl.col("associated_transcript").cast(pl.String)
    )

TALON_ID_to_GENCODE_V47 = read_refmap("nextflow_results/V47/compare/TALON_GENCODE_V39.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap")

TALON_SFARI = read_refmap("nextflow_results/V47/compare/TALON_SFARI.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap")

# Start here

n_transcripts_removed = TALON_ID_to_GENCODE_V47\
    .filter(pl.col("transcript_id").is_in(TALON_ID_to_GENCODE_V47.group_by("transcript_id").len().filter(pl.col("len")>1)["transcript_id"]))\
    .with_columns(
        pl.col("transcript_id").map_elements(lambda s: s.split("_")[0].split(".")[0], return_dtype=pl.String),
        pl.col("ref_id").map_elements(lambda s: s.split(".")[0], return_dtype=pl.String)
    )\
    .with_columns(
        is_same = (pl.col("transcript_id") == pl.col("ref_id"))
    )\
    .with_columns(
        pl.col("is_same").cast(pl.Int16)
    )\
    .group_by("transcript_id").agg(
        pl.col("is_same").sum().alias("is_same")
    )\
    .filter(pl.col("is_same")==0).shape[0]

print(f"{n_transcripts_removed} transcripts will be removed because they are not uniquely mapped to a ref_id and none of the matching ref_id matches the transcript_id")    

total_SFARI_tx = final_transcripts_gtf.unique("isoform").shape[0]

matched_SFARI_tx = final_transcripts_gtf\
    .unique("isoform")\
    .filter(
        pl.col("isoform").is_in(TALON_SFARI["ref_id"].unique())
    ).shape[0]
