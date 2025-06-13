#!/usr/bin/env python
from src.utils import read_gtf
import polars as pl

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Filter peptide GTF file for UCSC tracks')
    parser.add_argument('--annot_peptides_gtf', action='store', type=str, required=True)
    parser.add_argument('--output', action='store', type=str, required=True)
    params = parser.parse_args()

    annot_peptides_hybrid = read_gtf(params.annot_peptides_gtf, attributes=["detected", "transcript_id"])

    annot_peptides_hybrid\
        .filter(pl.col("detected") == "True")\
        .with_columns(
            score=pl.lit("."),
            frame=pl.lit("."),
        )\
        .drop("detected", "transcript_id")\
        .write_csv(params.output, include_header=False, separator="\t", quote_style="never")
    
if __name__ == "__main__":
    main()