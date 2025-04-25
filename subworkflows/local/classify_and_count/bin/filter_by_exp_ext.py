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

def get_polyA_tx(gtf, polyA_site):
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
    return validated_pbids
    
def get_CAGE_tx(gtf, reftss):
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
        .filter(
            (pl.col("seqname") == pl.col("chrom"))
        )\
        .unique("transcript_id")\
        ["transcript_id"].to_list()
    return validated_pbids

def get_exp_tx(expression, min_reads=5, min_n_sample=2):
    return expression\
        .with_columns(
            cs.numeric() > min_reads
        )\
        .filter(
            pl.sum_horizontal(cs.boolean()) > min_n_sample
        )\
        ["isoform"].to_list()

def main():
    parser = argparse.ArgumentParser(description='Filter isoforms based on expression and external support')
    parser.add_argument('--classification', action='store', dest='classification', type=str, required=True)
    parser.add_argument('--filtered_gff', action='store', dest='filtered_gff', type=str, required=True)
    parser.add_argument('--full_expression', action='store', dest='full_expression', type=str, required=True)
    parser.add_argument('--min_reads', action='store', dest='min_reads', type=float, default=5)
    parser.add_argument('--min_n_sample', action='store', dest='min_n_sample', type=int, default=2)
    parser.add_argument('--polyA_site', action='store', dest='polyA_site', type=str, required=True)
    parser.add_argument('--refTSS', action='store', dest='refTSS', type=str, required=True)
    
    params = parser.parse_args()

    # Read in datasets
    classification = get_classification(params.classification)
    filtered_gff = read_gtf(params.filtered_gff)
    expression = pl.read_parquet(params.full_expression)
    reftss = pl.read_csv(
            params.refTSS, 
            separator="\t", has_header = False, 
            new_columns=["chrom", "start", "end", "name", "score", "strand"]
        )
    polyA_site = pl.read_csv(
            params.polyA_site, 
            separator="\t", new_columns=["chrom", "start", "end", "name", "score", "strand"], 
            schema_overrides={"chrom": pl.String}
        )\
        .with_columns(
            pl.col("chrom").map_elements(lambda x: "".join(["chr", x]), return_dtype=pl.String).alias("chrom")
        )

    # Filtering isoforms based on expression
    expression = expression\
        .filter(
            pl.col("isoform").is_in(classification["isoform"])
        )

    expression = expression\
        .filter(
            pl.col("isoform").is_in(get_exp_tx(expression, 5, 2))
        )

    print(f"We have {classification.shape[0]} isoforms in the classification file")

    classification = classification\
        .filter(
            pl.col("isoform").is_in(expression["isoform"])
        )

    print(f"Keeping isofroms that have at least {5} reads in at least {2} samples, remaining {classification.shape[0]} isoforms")

    # For 5prime_fragment, keep those with polyA site support, for 3prime_fragment, keep those with CAGE peak support

    classification = classification\
        .with_columns(
            within_CAGE_peak = pl.col("isoform").is_in(get_CAGE_tx(filtered_gff, reftss)),
            within_polyA_site = pl.col("isoform").is_in(get_polyA_tx(filtered_gff, polyA_site))
        )
    
    five_prime_fragment_with_polyA = classification\
        .filter(
            pl.col("subcategory") == "5prime_fragment",
            pl.col("within_polyA_site")
        )["isoform"].to_list()
        
    three_prime_fragment_with_CAGE = classification\
        .filter(
            pl.col("subcategory") == "3prime_fragment",
            pl.col("within_CAGE_peak")
        )["isoform"].to_list()

    ISM_to_keep = five_prime_fragment_with_polyA + three_prime_fragment_with_CAGE

    classification = classification\
        .filter(
                (pl.col("structural_category") != "incomplete-splice_match") | (pl.col("isoform").is_in(ISM_to_keep))
            )

    print(f"Keeping ISMs without external support, remaining {classification.shape[0]} isoforms")

    # Write classification file
    classification.write_parquet("final_classification.parquet")

    print(f"Writing classification file to final_classification.parquet...")

    # Write gtf file

    print(f"Writing gtf file to final_transcripts.gtf...")

    filtered_gff\
        .filter(
            pl.col("transcript_id").is_in(classification["isoform"])
        )\
        .drop("transcript_id")\
        .write_csv("final_transcripts.gtf", separator="\t", include_header=False, quote_style="never")

    # Write expression file

    print(f"Writing expression file to final_expression.parquet...")

    expression\
        .filter(pl.col("isoform").is_in(classification["isoform"]))\
        .write_parquet("final_expression.parquet")
    
if __name__ == "__main__":
    main()