#!/usr/bin/env python3

import argparse
from src.utils import read_gtf
import polars as pl

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orfanage_gtf", action="store", dest="orfanage_gtf")
    parser.add_argument("--output", action="store", dest="output")

    args = parser.parse_args()

    gtf = read_gtf(args.orfanage_gtf)
    transcripts_to_keep = gtf.filter(pl.col("feature")=="CDS").unique("transcript_id").select("transcript_id")
    
    gtf\
        .filter(
            pl.col("transcript_id").is_in(transcripts_to_keep)
        )\
        .drop("transcript_id")\
        .write_csv(args.output, include_header=False, quote_style="never", separator="\t")
    
if __name__ == "__main__":
    main()