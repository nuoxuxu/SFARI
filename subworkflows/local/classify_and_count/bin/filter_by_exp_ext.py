#!/usr/bin/env python
import polars as pl
import polars.selectors as cs
import argparse

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

def get_classification(classification_path):
    print("Reading isoform classification file...")
    classification = pl.read_csv(classification_path, separator="\t", null_values=["NA"])
    print(f"{classification.shape[0]} isoforms are in classification file, these are the isoforms that passes pigeon filter")
    return classification

def main():
    parser = argparse.ArgumentParser(description='Get single_cell object from long read data')
    parser.add_argument('--classification', action='store', dest='classification', type=str, required=True)
    parser.add_argument('--filtered_gff', action='store', dest='filtered_gff', type=str, required=True)
    parser.add_argument('--full_expression', action='store', dest='full_expression', type=str, required=True)
    parser.add_argument('--min_reads', action='store', dest='min_reads', type=float, default=5)
    parser.add_argument('--min_n_sample', action='store', dest='min_n_sample', type=int, default=2)
    parser.add_argument('--polyA_site', action='store', dest='min_n_sample', type=str, required=True)
    parser.add_argument('--refTSS', action='store', dest='min_n_sample', type=str, required=True)
    
    params = parser.parse_args()

    expression = pl.read_parquet(params.read_stat, params.id_to_sample)
    classification = get_classification(params.classification)

    print(f"Filter out {expression.shape[0]-classification.shape[0]} that did not pass pigeon filter, remaining {classification.shape[0]} isoforms")
    
    expression = classification\
        .join(expression, on = "isoform", how = "left")\
        .select(expression.columns)
    
    classification = classification\
        .join(expression, on = "isoform", how = "left")\
        .select(classification.columns)
    
    print(f"Keeping isofroms that have at least {params.min_reads} reads in at least {params.min_n_sample} samples...")
    
    before = expression.shape[0]

    isoform = expression\
        .with_columns(
            cs.numeric() > 5
        )\
        .filter(
            pl.sum_horizontal(cs.boolean()) > 2
        )\
        .select("isoform")
    
    expression = expression.filter(pl.col("isoform").is_in(isoform["isoform"]))
    
    after = expression.shape[0]
    
    print(f"Filtered out {before - after} isoforms, remaining {expression.shape[0]} isoforms")

    print("Filtering the classification file...")

    classification = classification.filter(pl.col("isoform").is_in(expression["isoform"]))
    
    # Write classification file
    
    print(f"Writing classification file to {params.classification_output}...")
    classification.write_parquet(params.classification_output)

    # Write gtf file

    print(f"Writing gtf file to {params.gtf_output}...")
    read_gtf(params.filtered_gff)\
        .filter(
            pl.col("transcript_id").is_in(classification["isoform"])
        )\
        .drop("transcript_id")\
        .write_csv(params.gtf_output, separator="\t", include_header=False, quote_style="never")
    
    # Write expression file

    print(f"Writing expression file to {params.expression_output}...")
    expression\
        .filter(pl.col("isoform").is_in(classification["isoform"]))\
        .write_parquet(params.expression_output)
    
if __name__ == "__main__":
    main()