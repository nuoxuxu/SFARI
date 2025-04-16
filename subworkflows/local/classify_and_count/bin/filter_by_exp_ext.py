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

def remove_ISM_no_polyA(classification, gtf, polyA_site_path):

    print(f"{classification.shape[0]} isoforms are in classification file")

    polyA_site = pl.read_csv(
        polyA_site_path, 
        separator="\t", new_columns=["chrom", "start", "end", "name", "score", "strand"], 
        schema_overrides={"chrom": pl.String})\
        .with_columns(
            pl.col("chrom").map_elements(lambda x: "".join(["chr", x]), return_dtype=pl.String).alias("chrom")
        )

    validated_pbids = gtf\
    .filter(pl.col("feature")=="transcript")\
    .select(
        pl.col("seqname"),
        pl.col("transcript_id"),
        pos = pl.when(pl.col("strand")=="+")\
            .then(pl.col("end"))\
            .otherwise(pl.col("start"))
    )\
    .join_where(
        polyA_site,
        (pl.col("pos") >= pl.col("start")) &
        (pl.col("pos") <= pl.col("end"))
    )\
    .filter(
        pl.col("seqname")==pl.col("chrom")
    )\
    ["transcript_id"].to_list()
    
    classification =  classification\
        .filter(
            (pl.col("structural_category") != "incomplete-splice_match") | (pl.col("isoform").is_in(validated_pbids))
        )
    print(f"After filtering out ISM transcripts that don't have polyA site support, there remaining {classification.shape[0]} isoforms")
    return classification

def remove_ISM_no_refTSS(classification, gtf, refTSS_path):

    print(f"{classification.shape[0]} isoforms are in classification file")

    reftss = pl.read_csv(
        refTSS_path, 
        separator="\t", has_header = False, 
        new_columns=["chrom", "start", "end", "name", "score", "strand"]
        )
    
    validated_pbids = gtf\
    .filter(pl.col("feature")=="transcript")\
    .select(
        pl.col("seqname"),
        pl.col("transcript_id"),
        pos = pl.when(pl.col("strand")=="+")\
            .then(pl.col("start"))\
            .otherwise(pl.col("end"))
    )\
    .join_where(
        reftss,
        (pl.col("pos") >= (pl.col("start")-100)) &
        (pl.col("pos") <= (pl.col("end")+100))
    )\
    .filter(pl.col("seqname") == pl.col("chrom"))\
    ["transcript_id"].to_list()
    
    classification =  classification\
        .filter(
            (pl.col("structural_category") != "incomplete-splice_match") | (pl.col("isoform").is_in(validated_pbids))
        )
    print(f"After filtering out ISM that don't have rerfTSS support, there remaining {classification.shape[0]} isoforms")
    return classification    

def main():
    parser = argparse.ArgumentParser(description='Get single_cell object from long read data')
    parser.add_argument('--classification', action='store', dest='classification', type=str, required=True)
    parser.add_argument('--filtered_gff', action='store', dest='filtered_gff', type=str, required=True)
    parser.add_argument('--full_expression', action='store', dest='full_expression', type=str, required=True)
    parser.add_argument('--min_reads', action='store', dest='min_reads', type=float, default=5)
    parser.add_argument('--min_n_sample', action='store', dest='min_n_sample', type=int, default=2)
    parser.add_argument('--polyA_site', action='store', dest='polyA_site', type=str, required=True)
    parser.add_argument('--refTSS', action='store', dest='refTSS', type=str, required=True)
    
    params = parser.parse_args()

    classification = get_classification(params.classification)
    filtered_gff = read_gtf(params.filtered_gff)
    expression = pl.read_parquet(params.full_expression)
    
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

    print("Filtering the classification file based on the expression file...")

    classification = classification.filter(pl.col("isoform").is_in(expression["isoform"]))

    print("Filtering the filtered gff file based on the expression file...")

    filtered_gff = filtered_gff.filter(pl.col("transcript_id").is_in(expression["isoform"]))
    
    classification = remove_ISM_no_polyA(classification, filtered_gff, params.polyA_site)

    classification = remove_ISM_no_refTSS(classification, filtered_gff, params.refTSS)

    print("Filtering the filtered gff file based on the classification file accordingly...")
    filtered_gff = filtered_gff\
        .filter(pl.col("transcript_id").is_in(classification["isoform"]))
    
    print("Filtering the expression file based on the classification file accordingly...")
    expression = expression\
        .filter(pl.col("isoform").is_in(classification["isoform"]))
    
    # Write classification file
    
    print(f"Writing classification file to final_classification.parquet...")
    classification.write_parquet("final_classification.parquet")

    # Write gtf file

    print(f"Writing gtf file to final_transcripts.gtf...")
    filtered_gff\
        .drop("transcript_id")\
        .write_csv("final_transcripts.gtf", separator="\t", include_header=False, quote_style="never")
    
    # Write expression file

    print(f"Writing expression file to final_expression.parquet...")
    expression\
        .filter(pl.col("isoform").is_in(classification["isoform"]))\
        .write_parquet("final_expression.parquet")
    
if __name__ == "__main__":
    main()