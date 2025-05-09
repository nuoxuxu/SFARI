import polars as pl
from src.utils import gtf_to_SJ, read_gtf
import os

def get_novel_CSSs(predicted_cds_gtf, out, feature="exon"):
    """
    Get novel canonical splice sites present in the predicted GTF but not the GENCODE GTF.
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

    gencode_gtf = os.getenv("GENOMIC_DATA_DIR") + "/GENCODE/gencode.v47.annotation.gtf" # type: ignore
    GENCODE_SJ = read_gtf(gencode_gtf)\
        .filter(pl.col("feature") == feature)\
        .pipe(gtf_to_SJ)\
        .unpivot(
            index=["chrom", "strand", "transcript_id"],
            variable_name="start_or_end",
            value_name="coord"
        )\
        .unique(["chrom", "start_or_end", "coord"]).drop("strand")

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
        .drop("start_or_end")

    novel_splice_sites.write_csv(out)

get_novel_CSSs("nextflow_results/V47/orfanage/orfanage.gtf", "export/variant/novel_splice_sites_all.csv", feature="exon")
get_novel_CSSs("nextflow_results/V47/orfanage/orfanage.gtf", "export/variant/novel_splice_sites_cds.csv", feature="CDS")