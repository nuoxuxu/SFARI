#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
from src.utils import read_gff
import polars as pl
import argparse

def main():
    parser = argparse.ArgumentParser(description="Reformat TransDecoder genome GFF3 to be used in UCSC Genome Browser")
    parser.add_argument("--genome_gff3", action="store", type=str)
    parser.add_argument("--output", action="store", type=str)
    params = parser.parse_args()

    GFF3 = read_gff(params.genome_gff3, attributes=["Parent"])

    GFF3 = GFF3\
        .with_columns(
            pl.when(pl.col("feature")=="mRNA").then(pl.col("Parent").str.split("^").map_elements(lambda s: s[0], return_dtype=pl.String)).otherwise(pl.col("Parent").str.extract(r"^(.*)\.[^.]+$"))
        )\
        .with_columns(
            pl.when(pl.col("feature")!= "mRNA").then(pl.col("Parent").str.extract(r"^(.*)\.[^.]+$")).otherwise(pl.col("Parent")).alias("gene_id"),
            pl.when(pl.col("feature")=="mRNA").then(None).otherwise(pl.col("Parent")).alias("transcript_id")
        ).drop("Parent").drop_nulls("seqname")

    GFF3\
        .with_columns(
            attributes = pl.when(pl.col("feature")=="mRNA").then(pl.lit('gene_id "') + pl.col("gene_id") + pl.lit('";')).otherwise(pl.lit('gene_id "') + pl.col("gene_id") + pl.lit('";') + pl.lit('transcript_id "') + pl.col("transcript_id") + pl.lit('";'))
        )\
        .select(["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])\
        .write_csv(params.output, quote_style = "never", separator = "\t", include_header = False)
    
if __name__ == "__main__":
    main()