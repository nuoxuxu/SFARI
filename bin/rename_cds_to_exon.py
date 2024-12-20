#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python

import polars as pl
from src.utils import read_gtf
from src.single_cell import SingleCell

lr_bulk = SingleCell("nextflow_results/pbid_orf.h5ad")

isoforms_to_keep = lr_bulk.var\
    .filter(pl.col("ORF_type") == "complete")\
    ["base_isoform"].drop_nulls().unique()

genes_to_keep = isoforms_to_keep.cast(pl.String).str.split(".").map_elements(lambda s: "".join(s[0]+"."+s[1]))

sample_gtf = read_gtf("nextflow_results/transcripts_filtered.fasta.transdecoder.genome.gff3.gtf")

sample_gtf\
    .with_columns(
        pl.when(pl.col("feature")=="mRNA")\
            .then(pl.col("feature").map_elements(lambda s: "transcript", return_dtype=pl.String))\
            .otherwise(pl.col("feature"))
    )\
    .filter(
        pl.col("feature") == "mRNA"
    )