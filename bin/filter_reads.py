#!/usr/bin/env python3
import polars as pl
import sys
import argparse

def main():
    parser = argparse.ArgumentParser("Filter for reads for PITA analysis")
    parser.add_argument("--final_classification", action="store", type=str, required=True)
    parser.add_argument("--orfanage_gtf", action="store", type=str, required=True)
    parser.add_argument("--read_stat", action="store", type=str, required=True)
    parser.add_argument("--out", action="store", type=str, required=True)
    params = parser.parse_args()

    classification = pl.read_parquet(params.final_classification)
    classification = classification\
        .filter(pl.col("ref_exons")!=1)\
        .filter(pl.col("within_CAGE_peak"))

    orfanage = pl.read_csv(params.orfanage_gtf, separator="\t", comment_prefix="#", has_header=False)\
        .with_columns(pl.col("column_9").str.extract(rf'transcript_id "([^;]*)";').alias("transcript_id"))
    orfanage_transcript_id = orfanage["transcript_id"].unique()
    tx_to_keep = classification\
        .filter(
            pl.col("structural_category").is_in(["incomplete-splice_match", "full-splice_match", "novel_not_in_catalog", "novel_in_catalog"]),
            pl.col("isoform").is_in(orfanage_transcript_id.to_list())
        )["isoform"].unique()

    merged_collapsed_read_stat = pl.read_csv(params.read_stat, separator="\t")\
        .filter(pl.col("pbid").is_in(tx_to_keep.to_list()))

    CN_1_3 = pl.read_csv(sys.stdin.buffer, separator="\t", has_header=False, new_columns=["chr", "start", "end", "id", "score", "strand"])

    CN_1_3_filtered = CN_1_3\
        .join(merged_collapsed_read_stat, on="id", how="inner")\
        .join(classification["isoform", "associated_gene"], left_on="pbid", right_on="isoform", how="left")\
        .drop("pbid")

    CN_1_3_filtered.write_csv(params.out, include_header=False, separator="\t")

if __name__ == "__main__":
    main()