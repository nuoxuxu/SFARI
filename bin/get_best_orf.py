#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python

from src.single_cell import SingleCell
import polars as pl
from src.utils import read_gff
import argparse

parser = argparse.ArgumentParser("Test")
parser.add_argument("--h5ad_file", action='store', type=str, required=True)
parser.add_argument("--pred_orfs_gff3", action='store', type=str, required=True)
parser.add_argument("--output", action='store', type=str, required=True)
params = parser.parse_args()

lr_bulk = SingleCell(params.h5ad_file)

best_orf = lr_bulk.var.filter(pl.col("predicted_orf"))["isoform", "length"]\
    .with_columns(
        orf_frame = 1
    )\
    .rename({"isoform": "pb_acc", "length": "len"})

transdecoder_gff3 = read_gff(params.pred_orfs_gff3)\
    .drop_nulls("seqname")\
    .filter(pl.col("feature")== "CDS")\
    ["seqname", "start", "end"]\
    .rename({"seqname": "pb_acc", "start": "orf_start", "end": "orf_end"})

best_orf = best_orf\
    .join(
        transdecoder_gff3, on = "pb_acc", how = "left"
    )\
    .with_columns(
        orf_len = pl.col("orf_end") - pl.col("orf_start") + 1
    )

best_orf.write_csv(params.output, separator="\t")