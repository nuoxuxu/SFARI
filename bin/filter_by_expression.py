#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
import polars as pl
from src.utils import read_gtf
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Get single_cell object from long read data')
    parser.add_argument('--id_to_sample', action='store', dest='id_to_sample', type=str, required=True)
    parser.add_argument('--classification', action='store', dest='classification', type=str, required=True)
    parser.add_argument('--read_stat', action='store', dest='read_stat', type=str, required=True)
    parser.add_argument('--filtered_gff', action='store', dest='filtered_gff', type=str, required=True)
    parser.add_argument('--min_reads', action='store', dest='min_reads', type=float, default=5)
    parser.add_argument('--min_n_sample', action='store', dest='min_n_sample', type=int, default=2)
    parser.add_argument('--protein_coding_gene', action='store_true')
    parser.add_argument('--classification_output', action='store', dest='classification_output', type=str, required=True)
    parser.add_argument('--gtf_output', action='store', dest='gtf_output', type=str, required=True)
    params = parser.parse_args()

    id_to_sample = pl.read_csv(params.id_to_sample, separator = "\t", has_header = False, new_columns = ["id", "sample"])\
        .with_columns(
            pl.col("id").map_elements(lambda s: s.split("/")[0], return_dtype = pl.String),
            pl.col("sample").map_elements(lambda s: s.rsplit("/")[9].rsplit("_", 2)[0], return_dtype = pl.String)
        )\
        .to_pandas().set_index("id").to_dict()["sample"]
    
    classification = pl.read_csv(params.classification, separator="\t", null_values=["NA"])

    read_stat = pl.read_csv(params.read_stat, separator = "\t")\
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

    # Matching isoform in classification and read_stat

    read_stat = classification\
        .join(read_stat, on = "isoform", how = "left")\
        .select(read_stat.columns)

    classification = classification\
        .join(read_stat, on = "isoform", how = "left")\
        .select(classification.columns)
    
    # Write classification file
    
    classification.write_parquet(params.classification_output)

    # Write gtf file

    read_gtf(params.filtered_gff)\
        .filter(
            pl.col("transcript_id").is_in(classification["isoform"])
        )\
        .drop("transcript_id")\
        .write_csv(params.gtf_output, separator="\t", include_header=False, quote_style="never")
    
if __name__ == "__main__":
    main()