#!/usr/bin/env python
import polars as pl
import argparse

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
        .pivot(on = "id", index = "isoform", values = "count")\
        .fill_null(0)\
        .rename(id_to_sample)\
        ["isoform", "iPSC_1", "iPSC_2", "iPSC_3", "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3", "CN_1_2", "CN_1_3", "CN_2_1", "CN_2_2", "CN_3_1", "CN_3_2"]

    print(f"There are {expression.shape[0]} collapsed isoforms in the expression count matrix")

    return expression

def main():
    parser = argparse.ArgumentParser(description='Get full expression matrix from long read data')
    parser.add_argument('--read_stat', action='store', dest='read_stat', type=str, required=True)
    parser.add_argument('--id_to_sample', action='store', dest='id_to_sample', type=str, required=True)

    params = parser.parse_args()

    expression = get_expression(params.read_stat, params.id_to_sample)

    expression.write_parquet("full_expression.parquet")

if __name__ == "__main__":
    main()    