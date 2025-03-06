from src.utils import read_gtf
import polars as pl

orfanage_gtf = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")

exon_number = orfanage_gtf\
    .filter(pl.col("feature") == "exon")\
    .with_columns(
        exon_number = pl.col("transcript_id").rank(method="ordinal").over("transcript_id")
    ).select(["transcript_id", "start", "end", "strand", "exon_number"])

orfanage_gtf = orfanage_gtf\
    .join(
        exon_number,
        on = ["transcript_id", "start", "end", "strand"],
        how = "left"
    )\
    .with_columns(
        attributes = pl.when(pl.col("feature") == "exon")\
            .then(pl.col("attributes") + pl.lit(' exon_number ') + pl.col('exon_number').cast(pl.String) + pl.lit(';'))\
            .otherwise(pl.col("attributes"))
        
    )

orfanage_gtf.write_csv("export/orfanage_with_exon_number.gtf", separator="\t", include_header=False, quote_style="never")