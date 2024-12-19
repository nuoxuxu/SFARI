#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
import polars as pl
from src.single_cell import SingleCell
import argparse

def main():
    parser = argparse.ArgumentParser(description='Get the list of novel transcripts')
    parser.add_argument('--h5ad_file_orf', action='store', type=str, required=True)
    parser.add_argument('--genome_gff3', action='store', type=str, required=True)
    parser.add_argument('--output', action='store', type=str, required=True)
    params = parser.parse_args()

    lr_bulk = SingleCell(params.h5ad_file_orf)

    lr_bulk.var.filter(
        pl.col("structural_category").is_in(["novel_not_in_catalog", "novel_in_catalog"]),
        pl.col("isoform").is_in(lr_bulk.varm["ORF_annot"].filter(pl.col("predicted_orf"))["isoform"])
        )["isoform"]\
        .to_frame().write_csv(params.output, include_header=False)

if __name__ == "__main__":
    main()        