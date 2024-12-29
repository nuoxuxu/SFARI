#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
import polars as pl
from src.single_cell import SingleCell
import argparse

def main():
    parser = argparse.ArgumentParser(description="Add peptide support to single cell object")
    parser.add_argument("--h5ad_file", action="store", type=str, required=True)
    parser.add_argument("--peptide_file", action="store", type=str, required=True)
    parser.add_argument("--output", action="store", type=str, required=True)

lr_bulk = SingleCell(parser.h5ad_file)

percolator = pl.read_csv("nextflow_results/percolator.tsv", has_header=True, separator="\t")
percolator = percolator\
    .with_columns(
       proteinIds = percolator["proteinIds"].map_elements(lambda s: s.split(","))
    ).explode("proteinIds")\
    .filter(
        pl.col("q-value") < 0.05
    )\
    .with_columns(
        pl.col("proteinIds").str.split("|").map_elements(lambda s: s[1], return_dtype=pl.String)
    )

known_transcripts = lr_bulk.var.cast({"associated_transcript": pl.String}).filter(pl.col("associated_transcript").str.starts_with("ENST"))["isoform"].unique()

PSMId_mapped_to_novel_transcripts_uniquely = percolator\
    .with_columns(
        is_known = pl.col("proteinIds").str.starts_with("ENST").cast(pl.Boolean)
        )\
    .group_by("PSMId")\
    .agg(
        pl.col("is_known").sum()
    )\
    .filter(
        pl.col("is_known")==0
    )\
    .select("PSMId")

validated_pbids = percolator\
    .filter(
        pl.col("PSMId").is_in(PSMId_mapped_to_novel_transcripts_uniquely)
    )\
    .select("proteinIds").unique()

lr_bulk.var = lr_bulk.var\
    .with_columns(
        in_any_PSMs = pl.col("associated_transcript").is_in(percolator["proteinIds"].unique()) | pl.col("isoform").is_in(percolator["proteinIds"].unique()),
        in_novel_only_PSMs = pl.col("isoform").is_in(validated_pbids)
    )\
    .with_columns(
        validated_proteomics = pl.when(pl.col("structural_category").is_in(["novel_not_in_catalog", "novel_in_catalog"]))\
            .then(pl.col("in_novel_only_PSMs")).otherwise(pl.col("in_any_PSMs"))
    )

if __name__ == "__main__":
    main()