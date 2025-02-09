#!/usr/bin/env python3
from src.utils import read_gtf
import polars as pl
import argparse

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--orfanage_gtf", type=str, required=True)
    argparser.add_argument("--output", type=str, required=True)

    params = argparser.parse_args()

    orfanage_output = read_gtf(params.orfanage_gtf, attributes=["gene_id", "transcript_id"])

    orfanage_output = orfanage_output\
        .with_columns(
            gene_id = pl.when(pl.col("gene_id").is_null())\
                .then(pl.col("transcript_id").str.extract(r'(PB\.\d+)\.\d+'))\
                .otherwise(pl.col("gene_id"))
        )

    orfanage_output = orfanage_output\
        .with_columns(
            attributes = pl.when(pl.col("feature")!="transcript")\
                .then(pl.lit('gene_id "')+pl.col("gene_id")+pl.lit('"; transcript_id "')+pl.col("transcript_id")+pl.lit('";'))\
                .otherwise(pl.col("attributes"))
        )

    orfanage_output\
        .drop("gene_id", "transcript_id")\
        .write_csv(params.output, separator="\t", quote_style="never", include_header=False)
    
if __name__ == "__main__":
    main()