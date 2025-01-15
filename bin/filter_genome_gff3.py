#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
import polars as pl
from src.single_cell import SingleCell
import argparse
from src.utils import read_gff

def main():
    parser = argparse.ArgumentParser(description='Filter genome_gffs from TransDecoder with novel transripts')
    parser.add_argument('--mode', action='store', type=str, required=True)
    parser.add_argument('--h5ad_file_orf', action='store', type=str, required=True)
    parser.add_argument('--genome_gff3', action='store', type=str, required=True)
    parser.add_argument('--protein_classification_unfiltered', action='store', type=str, required=True)
    parser.add_argument('--output', action='store', type=str, required=True)
    params = parser.parse_args()

    lr_bulk = SingleCell(params.h5ad_file_orf)

    genome_gff3 = read_gff(params.genome_gff3, attributes=["ID", "Parent"])\
        .drop_nulls("seqname")\
        .with_columns(
            pl.when(pl.col("feature")=="gene")\
                .then(pl.col("ID").str.split("^").map_elements(lambda x: x[0], return_dtype=pl.String))\
                .when(pl.col("feature")=="mRNA")\
                .then(pl.col("ID").str.extract(r"^(.*)\.p\d+$"))\
                .otherwise(pl.col("ID")),
            pl.when(pl.col("feature")=="mRNA")\
                .then(pl.col("Parent").str.extract(r"^(.*)\^chr"))\
                .when(pl.col("feature").is_in(["three_prime_UTR", "five_prime_UTR", "CDS", "exon"]))\
                .then(pl.col("Parent").str.extract(r"^(.*)\.p\d+$"))\
                .otherwise(pl.col("Parent")),
            exon_number = pl.when(pl.col("feature")=="exon")\
                .then(pl.col("ID").str.extract(r"exon(\d+)"))\
                .otherwise(0)
        )\
        .with_columns(
            attributes = pl.when(pl.col("feature")=="gene")\
                .then(pl.lit("ID=")+pl.col("ID"))\
                .when(pl.col("feature").is_in(["mRNA", "five_prime_UTR", "three_prime_UTR", "CDS"]))\
                .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent"))\
                .when(pl.col("feature")=="exon")\
                .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent")+pl.lit(";exon_number=")+pl.col("exon_number").cast(pl.String))
        )\
        .with_columns(
            transcript_id = pl.when(pl.col("feature")=="gene")\
                .then(pl.lit(None))\
                .when(pl.col("feature")=="mRNA")\
                .then(pl.col("ID"))\
                .otherwise(pl.col("Parent"))
        )

    novel_pbids = lr_bulk.var\
        .join(
        pl.read_csv(params.protein_classification_unfiltered, separator = "\t")\
            .rename({"pb": "isoform"})\
            ["isoform", "protein_classification_base"],
        on = "isoform", how = "left"
        )\
        .filter(pl.col("ORF_type")=="complete")\
        .with_columns(
            pl.col("base_isoform").cast(pl.String)
        )\
        .filter(
            pl.col("isoform") == pl.col("base_isoform"), 
            pl.col("protein_classification_base").is_in(["pNNC", "pNIC"])
        )["isoform"]

    novel_gene_ids = novel_pbids.str.split(".").map_elements(lambda x: "".join(x[0]+"." + x[1])).unique()

    base_iso_pbids = lr_bulk.var\
        .filter(pl.col("ORF_type")=="complete")\
        .with_columns(
            pl.col("base_isoform").cast(pl.String)
        )\
        .filter(
            pl.col("isoform") == pl.col("base_isoform")
        )["isoform"]
    
    base_iso_gene_ids = base_iso_pbids.str.split(".").map_elements(lambda x: "".join(x[0]+"." + x[1])).unique()

    if params.mode == "hybrid":
        genome_gff3 = genome_gff3\
            .filter(
                ((pl.col("feature")=="gene") & pl.col("ID").is_in(novel_gene_ids)) | ((pl.col("feature")!="gene") & pl.col("transcript_id").is_in(novel_pbids))
            )
    elif params.mode == "pacbio":
        genome_gff3 = genome_gff3\
            .filter(
                ((pl.col("feature")=="gene") & pl.col("ID").is_in(base_iso_gene_ids)) | ((pl.col("feature")!="gene") & pl.col("transcript_id").is_in(base_iso_pbids))
            )        

    genome_gff3\
        .with_columns(
            attributes = pl.when(pl.col("feature")=="gene")\
                .then(pl.lit("ID=")+pl.col("ID"))\
                .when(pl.col("feature").is_in(["mRNA", "five_prime_UTR", "three_prime_UTR", "CDS"]))\
                .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent"))\
                .when(pl.col("feature")=="exon")\
                .then(pl.lit("ID=")+pl.col("ID")+pl.lit(";Parent=")+pl.col("Parent")+pl.lit(";exon_number=")+pl.col("exon_number"))        
        )\
        .drop(["ID", "Parent", "exon_number", "transcript_id"])\
        .write_csv(params.output, include_header=False, quote_style="never", separator="\t")  


if __name__ == "__main__":
    main()