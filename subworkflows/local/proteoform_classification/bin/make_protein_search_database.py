#!/usr/bin/env python3
import polars as pl
import argparse
from src.utils import read_fasta, write_fasta, read_gtf, collapse_isoforms_to_proteoforms

def main():
    
    parser = argparse.ArgumentParser(description='Filter genome_gffs from TransDecoder with novel transripts')
    parser.add_argument('--mode', action='store', type=str, required=True)
    parser.add_argument('--protein_classification', action='store', type=str, required=True)
    parser.add_argument('--filtered_predicted_cds_gtf', action='store', type=str, required=True)
    parser.add_argument('--peptide_fasta', action='store', type=str, required=True)
    parser.add_argument('--translation_fasta', action='store', type=str, required=True)
    parser.add_argument('--output', action='store', type=str, required=True)
    params = parser.parse_args()

    protein_classification = pl.read_csv(params.protein_classification, separator = "\t")
    predicted_cds_gtf = read_gtf(params.filtered_predicted_cds_gtf)
    isoforms_to_proteoforms = collapse_isoforms_to_proteoforms(predicted_cds_gtf)

    peptide_fasta = read_fasta(params.peptide_fasta)
    translation_fasta = read_fasta(params.translation_fasta)\
        .with_columns(
            pl.col("transcript_id").map_elements(lambda x: x.split("|")[1])
        )

    if params.mode == "hybrid":
        novel_pbids = protein_classification\
            .join(
                isoforms_to_proteoforms.rename({"isoform": "pb"}),
                on = "pb",
                how = "left"
            )\
            .filter(
                pl.col("pb").is_in(predicted_cds_gtf.filter(pl.col("feature")=="transcript")["transcript_id"]),
                pl.col("protein_classification_base").is_in(["pNNC", "pNIC"])
            )\
            .unique("base_isoform")\
            .select("base_isoform")

        novel_fasta = peptide_fasta\
            .filter(
                pl.col("transcript_id").is_in(novel_pbids)
            )
    
        write_fasta(pl.concat([novel_fasta, translation_fasta], how="vertical"), "transcript_id", "seq", params.output)

    elif params.mode == "pacbio":
        
        base_isoforms = protein_classification\
            .join(
                isoforms_to_proteoforms.rename({"transcript_id": "pb"}),
                on = "pb",
                how = "left"
            )\
            .filter(
                pl.col("pb").is_in(predicted_cds_gtf.filter(pl.col("feature")=="transcript")["transcript_id"])
            ).unique("base_isoform").select("base_isoform")
        
        write_fasta(peptide_fasta.filter(pl.col("transcript_id").is_in(base_isoforms)), "transcript_id", "seq", params.output)

if __name__ == "__main__":
    main()