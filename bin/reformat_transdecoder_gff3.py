#!/usr/bin/env python

from src.utils import read_gff
import polars as pl
import sys

transdecoder_genome_gff3 = sys.argv[1]
out = sys.argv[2]

GFF3 = read_gff(transdecoder_genome_gff3, attributes=["Parent"])

GFF3 = GFF3\
    .with_columns(
        pl.when(pl.col("type")=="mRNA").then(pl.col("Parent").str.split("^").map_elements(lambda s: s[0])).otherwise(pl.col("Parent").str.extract("^(.*)\.[^.]+$"))
    )\
    .with_columns(
        pl.when(pl.col("type")!= "mRNA").then(pl.col("Parent").str.extract("^(.*)\.[^.]+$")).otherwise(pl.col("Parent")).alias("gene_id"),
        pl.when(pl.col("type")=="mRNA").then(None).otherwise(pl.col("Parent")).alias("transcript_id")
    ).drop("Parent").drop_nulls("seqid")

GFF3\
    .with_columns(
        attributes = pl.when(pl.col("type")=="mRNA").then(pl.lit('gene_id "') + pl.col("gene_id") + pl.lit('";')).otherwise(pl.lit('gene_id "') + pl.col("gene_id") + pl.lit('";') + pl.lit('transcript_id "') + pl.col("transcript_id") + pl.lit('";'))
    )\
    .select(["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])\
    .write_csv(out, quote_style = "never", separator = "\t", include_header = False)