#!/usr/bin/env python
import polars as pl

def read_gtf(file, attributes=["transcript_id"], keep_attributes=True):
    if keep_attributes:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                )
    else:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                ).drop("attributes")
    
def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Filter orfanage GTF file for UCSC tracks')
    parser.add_argument('--orfanage_gtf', action='store', type=str, required=True)
    parser.add_argument('--output', action='store', type=str, required=True)
    params = parser.parse_args()
    
    orfanage_gtf = read_gtf(params.orfanage_gtf)
    
    orfanage_gtf\
        .filter(
            pl.col("feature").is_in(["exon", "CDS"])
        )\
        .drop("transcript_id")\
        .write_csv(params.output, separator="\t", include_header=False, quote_style="never")
    
if __name__ == "__main__":
    main()