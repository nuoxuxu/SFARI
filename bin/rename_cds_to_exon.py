#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python

import polars as pl
from src.utils import read_gtf
from src.single_cell import SingleCell
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('h5ad_file_orf', type=str, help='Path to the PBID ORF file')
    parser.add_argument('genome_gff3_gtf', type=str, help='Path to the TransDecoder genome GTF file')
    parser.add_argument('reference_gtf', type=str, help='Path to the reference GTF file')
    parser.add_argument('genome_gff3_gtf_renamed', type=str, help='Path to the renamed TransDecoder genome GTF file')
    parser.add_argument('reference_gtf_renamed', type=str, help='Path to the renamed reference GTF file')
    params = parser.parse_args()

    def get_full_gtf_from_exon(CDS_only, output_path):
        transcript_coords = CDS_only\
            .group_by(
                pl.col("transcript_id")
            )\
            .agg(
                pl.col("start").min().alias("start"),
                pl.col("end").max().alias("end")
            )

        transcript_only = CDS_only.unique("transcript_id")\
            .join(transcript_coords, on="transcript_id")\
            .with_columns(
                pl.col("feature").replace("exon", "transcript"),
                start = pl.col("start_right"),
                end = pl.col("end_right"),
            ).drop("start_right", "end_right")

        return pl.concat([transcript_only, CDS_only], how="vertical")\
            .with_columns(
                feature = pl.col("feature").cast(pl.Enum(["transcript", "exon"]))
            )\
            .sort("transcript_id", "feature", "start")\
            .drop("gene_id", "transcript_id")\
            .write_csv(output_path, separator="\t", incldue_header=False, quote_style="never")

    lr_bulk = SingleCell(params.h5ad_file_orf)

    # Filter for isoforms that represent unique ORFs

    isoforms_to_keep = lr_bulk.var\
        .filter(pl.col("ORF_type") == "complete")\
        ["base_isoform"].drop_nulls().unique()

    # genes_to_keep = isoforms_to_keep.cast(pl.String).str.split(".").map_elements(lambda s: "".join(s[0]+"."+s[1]))

    sample_CDS_only = read_gtf(params.genome_gff3_gtf, attributes=["gene_id", "transcript_id"])\
        .filter(
            pl.col("feature") == "CDS",
            pl.col("transcript_id").is_in(isoforms_to_keep)
        )\
        .with_columns(
            pl.col("feature").replace("CDS", "exon")
        )

    reference_CDS_only = read_gtf(params.reference_gtf, attributes=["gene_id", "transcript_id", "transcript_type"])\
        .filter(
            pl.col("feature") == "CDS",
            pl.col("transcript_type") == "protein_coding"
        )\
        .with_columns(
            pl.col("feature").replace("CDS", "exon")
        ).drop("transcript_type")
        

    get_full_gtf_from_exon(reference_CDS_only, params.genome_gff3_gtf_renamed)\
        
    get_full_gtf_from_exon(sample_CDS_only, params.reference_gtf_renamed)