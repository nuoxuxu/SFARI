#!/usr/bin/env python
import polars as pl
from src.utils import gtf_to_SJ, read_gtf
import argparse

def get_CSSs(predicted_cds_gtf, annotation_gtf, out, novel=True, feature="exon"):
    """
    Get novel canonical splice sites present in the predicted GTF but not the GENCODE GTF.
    A splice site is considered novel if it is not present in the GENCODE GTF, whether the 
    matching splice site is novel is not taken into consideration.
    Args:
        predicted_cds_gtf (str): Path to the predicted GTF file.
        out (str): Path to the output file.
        feature (str): Feature to filter on. Default is "exon".
    """
    
    final_transcripts_SJ = read_gtf(predicted_cds_gtf)\
        .filter(pl.col("feature") == feature)\
        .pipe(gtf_to_SJ)\
        .unpivot(
            index=["chrom", "strand", "transcript_id"],
            variable_name="start_or_end",
            value_name="coord"
        )\
        .unique(["chrom", "start_or_end", "coord"]).drop("strand")

    gencode_gtf = annotation_gtf
    GENCODE_SJ = read_gtf(gencode_gtf)\
        .filter(pl.col("feature") == feature)\
        .pipe(gtf_to_SJ)\
        .unpivot(
            index=["chrom", "strand", "transcript_id"],
            variable_name="start_or_end",
            value_name="coord"
        )\
        .unique(["chrom", "start_or_end", "coord"]).drop("strand")

    if novel:
        novel_splice_sites = final_transcripts_SJ\
            .join(
                GENCODE_SJ, 
                on=["chrom", "start_or_end", "coord"], 
                how="left",
                coalesce=True
                )\
            .filter(pl.col("transcript_id_right").is_null())\
            .drop("transcript_id_right")\
            .with_columns(
                pl.when(pl.col("start_or_end") == "start")\
                    .then(pl.col("coord") + 1)\
                    .otherwise(pl.col("coord") - 1).alias("coord_1")
            )

        novel_splice_sites = novel_splice_sites\
            .with_columns(
                pl.when(pl.col("start_or_end") == "start").then(pl.col("coord")).otherwise(pl.col("coord_1")).alias("start"),
                pl.when(pl.col("start_or_end") == "end").then(pl.col("coord")).otherwise(pl.col("coord_1")).alias("end")
            )\
            .drop("coord", "coord_1")
        
        novel_splice_sites.write_csv(out)
    else:
        known_splice_sites = GENCODE_SJ\
            .with_columns(
                pl.when(pl.col("start_or_end") == "start")\
                    .then(pl.col("coord") + 1)\
                    .otherwise(pl.col("coord") - 1).alias("coord_1")
            )\
            .with_columns(
                pl.when(pl.col("start_or_end") == "start").then(pl.col("coord")).otherwise(pl.col("coord_1")).alias("start"),
                pl.when(pl.col("start_or_end") == "end").then(pl.col("coord")).otherwise(pl.col("coord_1")).alias("end")
            )\
            .drop("coord", "coord_1")

        known_splice_sites.write_csv(out)

def main():
    parser = argparse.ArgumentParser(description='Get novel splice sites from predicted GTF')
    parser.add_argument('--annotation_gtf', action='store', type=str, required=True, help='Path to the GENCODE annotationGTF')
    parser.add_argument('--predicted_cds_gtf', action='store', type=str, required=True, help='Path to predicted CDS GTF')

    params = parser.parse_args()

    get_CSSs(params.predicted_cds_gtf, params.annotation_gtf, "novel_splice_sites_all.csv", feature="exon")
    get_CSSs(params.predicted_cds_gtf, params.annotation_gtf, "novel_splice_sites_cds.csv", feature="CDS")
    get_CSSs(params.predicted_cds_gtf, params.annotation_gtf, "known_splice_sites_all.csv", novel=False, feature="exon")
    get_CSSs(params.predicted_cds_gtf, params.annotation_gtf, "known_splice_sites_cds.csv", novel=False, feature="CDS")


if __name__ == "__main__":
    main()