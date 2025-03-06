from src.utils import read_gtf
import polars as pl
import polars.selectors as cs

orfanage_gtf = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")

exons = orfanage_gtf\
    .filter(pl.col("feature") == "exon")\
    .with_columns(
        exon_number = pl.col("transcript_id").rank(method="ordinal").over("transcript_id")
    )

other_level2_features = orfanage_gtf\
    .filter(
        pl.col("feature").is_in(["start_codon", "stop_codon", "CDS"])
    )\
    .join_where(
        exons.select(["transcript_id", "start", "end", "strand", "exon_number"]),
        pl.col("start") >= pl.col("start_right"),
        pl.col("end") <= pl.col("end_right"),
        pl.col("strand") == pl.col("strand_right"),
        pl.col("transcript_id") == pl.col("transcript_id_right")
    )\
    .drop(
        cs.ends_with("_right")
    )

orfanage_gtf_with_exon_num = pl.concat([exons, other_level2_features])\
    .sort(["transcript_id", "start"])\
    .with_columns(exon_number = pl.col("exon_number").cast(pl.String))\
    .with_columns(
        attributes = pl.col("attributes") + pl.lit(' exon_number "') + pl.col("exon_number") + pl.lit('";')
    ).drop(["transcript_id", "exon_number"])

orfanage_gtf_with_exon_num.write_csv("export/orfanage_with_exon_number.gtf", separator="\t", include_header=False, quote_style="never")

