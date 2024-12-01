#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
import polars as pl
from src.single_cell import SingleCell
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser(description='Filter long read data')
    parser.add_argument('--filtered_classification', action='store', dest='filtered_classification', type=str, required=True)
    parser.add_argument('--h5ad_file', action='store', dest='h5ad_file', type=str, required=True)
    parser.add_argument('--output', action='store', dest='output', type=str, required=True)
    parser.add_argument('--min_reads', action='store', dest='min_reads', type=float, default=5)
    parser.add_argument('--min_n_sample', action='store', dest='min_n_sample', type=int, default=2)
    params = parser.parse_args()

    lr_bulk = SingleCell(params.h5ad_file)
    classification = pl.read_csv(params.filtered_classification, separator="\t", null_values="NA")

    # Filter pacbio id count matrix to keep only transcripts in `classification`
    # Filter Settings:
    #   -a,--polya-percent     FLOAT  Adenine percentage at genomic 3' end to flag an isoform as intra-priming. [0.6]
    #   -r,--polya-run-length  INT    Continuous run-A length at genomic 3' end to flag an isoform as intra-priming. [6]
    #   -m,--max-distance      INT    Maximum distance to an annotated 3' end to preserve as a valid 3' end. [50]
    #   -c,--min-cov           INT    Minimum junction coverage for each isoform. Only used if min_cov field is not 'NA'. [3]
    #   --mono-exon                   Filter out all mono-exonic transcripts.
    #   --skip-junctions              Skip junctions.txt filtering.

    lr_bulk = lr_bulk.filter_var( pl.col("pbid").is_in(classification["isoform"]))
    lr_bulk = lr_bulk[:, np.where(np.sum(lr_bulk.X > params.min_reads, axis = 0) > params.min_n_sample)[0]]

    lr_bulk.save(params.output)

if __name__ == "__main__":
    main()