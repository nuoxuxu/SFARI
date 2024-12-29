#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
from src.utils import read_gff
import polars as pl
import argparse

def main():
    parser = argparse.ArgumentParser(description='Reformat TrandDecoder gff3 so that it works with GenomeKit')
    parser.add_argument('--genome_gff3', action='store', type=str, required=True)
    parser.add_argument('--output', action='store', dest='output', type=str, required=True)
    params = parser.parse_args()

    read_gff(params.genome_gff3, attributes=["ID", "Parent"])\
        .with_columns(
            pl.col("feature").replace({"gene": "gene", "five_prime_UTR": "five_prime_UTR", "exon": "exon", "mRNA": "transcript", "CDS": "CDS", "three_prime_UTR": "three_prime_UTR"})
        ).drop_nulls("feature")\
        .with_columns(
            pl.when(pl.col("feature")=="gene")\
                .then(pl.col("ID").str.split("^").map_elements(lambda x: x[0], return_dtype=pl.String))\
                .when(pl.col("feature")=="transcript")\
                .then(pl.col("ID").str.extract(r"^(.*)\.p\d+$"))\
                .otherwise(pl.col("ID")),
            pl.when(pl.col("feature")=="transcript")\
                .then(pl.col("Parent").str.extract(r"^(.*)\^chr"))\
                .when(pl.col("feature").is_in(["three_prime_UTR", "five_prime_UTR", "CDS", "exon"]))\
                .then(pl.col("Parent").str.extract(r"^(.*)\.p\d+$"))\
                .otherwise(pl.col("Parent")),
            exon_number = pl.when(pl.col("feature")=="exon")\
                .then(pl.col("ID").str.extract(r"exon(\d+)"))\
                .otherwise(0)
        )\
        .with_columns(
            gene_id = pl.when(pl.col("feature")=="gene")\
            .then(pl.col("ID"))\
            .when(pl.col("feature")=="transcript")\
            .then(pl.col("Parent"))\
            .otherwise(pl.col("Parent").str.split(".").map_elements(lambda x: ".".join([x[0], x[1]]), return_dtype=pl.String)),

            transcript_id = pl.when(pl.col("feature")=="gene")\
            .then(pl.lit(None))\
            .when(pl.col("feature")=="transcript")\
            .then(pl.col("ID"))\
            .otherwise(pl.col("Parent"))
        )\
        .with_columns(
            attributes = pl.when(pl.col("feature")=="gene")\
                .then(pl.lit("ID=")+pl.col("ID"))\
                .when(pl.col("feature").is_in(["transcript", "five_prime_UTR", "three_prime_UTR", "CDS"]))\
                .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent"))\
                .when(pl.col("feature")=="exon")\
                .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent")+pl.lit(";exon_number=")+pl.col("exon_number").cast(pl.String))
        )\
        .drop("exon_number", "ID", "Parent", "gene_id", "transcript_id")\
        ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']\
        .write_csv(params.output, include_header=False, quote_style="never", separator="\t")
    
if __name__ == "__main__":
    main()    