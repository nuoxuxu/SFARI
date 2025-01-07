#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
from src.utils import read_gff
import polars as pl
import argparse

def main():

    parser = argparse.ArgumentParser(description="Reformat TransDecoder genome GFF3 to be used in UCSC Genome Browser")
    parser.add_argument("--genome_gff3", action="store", type=str)
    parser.add_argument("--output", action="store", type=str)
    params = parser.parse_args()

    GFF3 = read_gff(params.genome_gff3, attributes=["ID", "Parent"]).drop_nulls("seqname")

    # Because we use GENCODE annotations and GENCODE annotations contain "transcript" as a feature,
    # Rename "mRNA" to "transcript" in the third column

    GFF3 = GFF3\
        .with_columns(
            feature = pl.when(pl.col("feature")=="mRNA").then(pl.lit("transcript")).otherwise(pl.col("feature"))
        )

    # Get attributes ID, Parent and exon_number

    GFF3 = GFF3\
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
    )

    # Get attributes gene_id and transcript_id

    GFF3 = GFF3\
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
        )

    GFF3\
        .with_columns(
            attributes = pl.when(pl.col("feature")=="gene")\
                .then(pl.lit('ID=') + pl.col('ID') + pl.lit(';gene_id=') + pl.col('gene_id'))\
                .otherwise(pl.lit('ID=') + pl.col('ID') + pl.lit(';Parent=') + pl.col('Parent') + pl.lit(';gene_id=') + pl.col('gene_id') + pl.lit(';transcript_id') + pl.col('transcript_id'))
        ).drop("ID", "Parent", "exon_number", "gene_id", "transcript_id")\
        .write_csv(params.output, quote_style = "never", separator = "\t", include_header = False)

if __name__ == "__main__":
    main()