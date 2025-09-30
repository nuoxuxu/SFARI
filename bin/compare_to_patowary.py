#!/usr/bin/env python3
import polars as pl
from src.utils import read_gtf
import numpy as np
import polars.selectors as cs
import os
import argparse

def read_refmap(file):
    out = pl.read_csv(file, separator="\t")\
    .filter(pl.col("class_code")=="=")\
    .with_columns(
        pl.col("qry_id_list").str.split(",").map_elements(lambda s: [e.split("|")[1] for e in s], return_dtype=pl.List(pl.String)).alias("qry_id_list"))\
    .explode("qry_id_list")\
    .select(["qry_id_list", "ref_id"])\
    .rename({"qry_id_list": "transcript_id"})
    return out

def main():
    parser = argparse.ArgumentParser(description="Compare final classification to Patowary et al. 2024")
    parser.add_argument("--classification", type=str, required=True, help="Path to final classification file")
    parser.add_argument("--TALON_SFARI_refmap", type=str, required=True, help="Path to TALON to SFARI refmap file")
    parser.add_argument("--output", type=str, required=False, help="Output path", default="nextflow_results/V47/")
    args = parser.parse_args()

    classification = pl.read_parquet(args.classification)

    TALON_SFARI = read_refmap(args.TALON_SFARI_refmap)

    classification = classification\
        .join(
            TALON_SFARI.unique("ref_id"), 
            left_on="isoform",
            right_on="ref_id",
            how = "left"
        )\
        .rename({"transcript_id": "patowary_isoform_id"})

    classification.write_csv(args.output, separator="\t", null_value="NA")

    n_matched = classification\
        .filter(pl.col("structural_category").is_in(["novel_in_catalog", "novel_not_in_catalog"]))\
        .filter(pl.col("patowary_isoform_id").is_null().not_())\
        .shape[0]

    n_total = classification\
        .filter(pl.col("structural_category").is_in(["novel_in_catalog", "novel_not_in_catalog"]))\
        .shape[0]

    print(f"{(n_matched/n_total)*100}% of the novel transcripts in SFARI are also in Patowary")

if __name__ == "__main__":
    main()