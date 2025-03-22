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

def get_expression(read_stat_path, id_to_sample_path):
    
    print("Getting mapping from PB IDs to sample names...")

    id_to_sample = pl.read_csv(id_to_sample_path, separator = "\t", has_header = False, new_columns = ["id", "sample"])\
        .with_columns(
            pl.col("id").map_elements(lambda s: s.split("/")[0], return_dtype = pl.String),
            pl.col("sample").map_elements(lambda s: s.rsplit("/")[9].rsplit("_", 2)[0], return_dtype = pl.String)
        )\
        .to_pandas().set_index("id").to_dict()["sample"]    

    assert len(id_to_sample.items()) == 15, "There should be 15 samples in the id_to_sample file"

    print("Computing expression count matrix...")

    expression = pl.read_csv(read_stat_path, separator = "\t")\
        .rename({"pbid": "isoform"})\
        .with_columns(
            pl.col("id").map_elements(lambda s: s.split("/")[0], return_dtype = pl.String),pl.lit(1).alias("count")
        )\
        .group_by("id", "isoform")\
        .agg(pl.sum("count"))\
        .sort("isoform", descending=False)\
        .pivot("id", index = "isoform", values = "count")\
        .fill_null(0)\
        .rename(id_to_sample)\
        ["isoform", "iPSC_1", "iPSC_2", "iPSC_3", "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3", "CN_1_2", "CN_1_3", "CN_2_1", "CN_2_2", "CN_3_1", "CN_3_2"]

    print(f"There are {expression.shape[0]} collapsed isoforms in the expression count matrix")

    return expression

def get_classification(classification_path):

    print("Reading isoform classification file...")

    classification = pl.read_csv(classification_path, separator="\t", null_values=["NA"])

    print(f"{classification.shape[0]} isoforms are in classification file, these are the isoforms that passes pigeon filter")

    return classification

def main():
    parser = argparse.ArgumentParser(description='Get single_cell object from long read data')
    parser.add_argument('--annotation_gtf', action='store', dest='annotation_gtf', type=str, required=True)
    parser.add_argument('--id_to_sample', action='store', dest='id_to_sample', type=str, required=True)
    parser.add_argument('--classification', action='store', dest='classification', type=str, required=True)
    parser.add_argument('--read_stat', action='store', dest='read_stat', type=str, required=True)
    parser.add_argument('--filtered_gff', action='store', dest='filtered_gff', type=str, required=True)
    parser.add_argument('--min_reads', action='store', dest='min_reads', type=float, default=5)
    parser.add_argument('--min_n_sample', action='store', dest='min_n_sample', type=int, default=2)
    parser.add_argument('--protein_coding_gene', action='store_true')
    parser.add_argument('--classification_output', action='store', dest='classification_output', type=str, required=True)
    parser.add_argument('--gtf_output', action='store', dest='gtf_output', type=str, required=True)
    parser.add_argument('--expression_output', action='store', dest='expression_output', type=str, required=True)
    
    params = parser.parse_args()

    expression = get_expression(params.read_stat, params.id_to_sample)
    classification = get_classification(params.classification)

    expression.write_parquet("fulll_expression.parquet")

    print(f"Filter out {expression.shape[0]-classification.shape[0]} that did not pass pigeon filter, remaining {classification.shape[0]} isoforms")
    
    expression = classification\
        .join(expression, on = "isoform", how = "left")\
        .select(expression.columns)

    print("Matching isoform in classification and read_stat...")

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

    # Filter by protein-coding genes

    if params.protein_coding_gene:
        print("Keeping transcripts associated with protein-coding genes...")
        gencode_gtf = read_gtf(params.annotation_gtf, attributes = ["gene_name", "gene_type"]).rename({"gene_name" : "associated_gene"}).unique("associated_gene")
        before = classification.shape[0]
        classification = classification\
            .join(
                gencode_gtf, 
                on = "associated_gene",
                how = "left"
            )\
            .filter(
                pl.col("gene_type") == "protein_coding"
            )\
            .drop("gene_type")
        after = classification.shape[0]
        print(f"Filtered out {before - after} isoforms that are not associated with protein-coding genes, remaining {classification.shape[0]} isoforms")
    
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