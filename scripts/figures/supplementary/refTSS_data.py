import polars as pl
import seaborn as sns
from src.utils import read_gtf

classification = pl.read_parquet("nextflow_results/V47/final_classification.parquet")
reftss = pl.read_csv("data/liftovered_mm39_to_hg38_peaks_overlapped_reftss_hg38_500bp.bed", separator="\t", has_header = False, new_columns=["seqname", "start", "end", "name", "score", "strand"])
gtf = read_gtf("nextflow_results/V47/final_transcripts.gtf")

validated_pbids = gtf\
    .filter(pl.col("feature")=="exon")\
    .group_by("transcript_id")\
    .agg(
        pl.col("seqname").map_elements(lambda x : x[0], return_dtype=pl.String),
        pl.col("strand").map_elements(lambda x : x[0], return_dtype=pl.String),
        pl.col("start").min(),
        pl.col("end").max()
    )\
    .select(
        pl.col("seqname"),
        pl.col("transcript_id"),
        pos = pl.when(pl.col("strand")=="+")\
            .then(pl.col("start"))\
            .otherwise(pl.col("end"))
    )\
    .join(
        reftss,
        on = "seqname",
        how = "inner"
    )\
    .filter(
        (pl.col("start") <= pl.col("pos") + 100) &
        (pl.col("end")   >= pl.col("pos") - 100)
    )\
    .unique("transcript_id")\
    .select("transcript_id")

classification\
    .with_columns(
        CAGE = pl.col("isoform").is_in(validated_pbids),
        structural_category2 = pl.when(pl.col("structural_category").is_in(["full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"]))\
            .then(pl.col("structural_category"))\
            .otherwise(pl.lit("Other"))
    )\
    .with_columns(
        structural_category2 = pl.when(pl.col("subcategory") == "3prime_fragment")\
            .then(pl.lit("3prime_fragment"))\
            .when(pl.col("subcategory") == "5prime_fragment")\
            .then(pl.lit("5prime_fragment"))\
            .otherwise(pl.col("structural_category2"))
    )\
    .group_by(["structural_category2", "CAGE"])\
    .len()\
    .group_by("structural_category2")\
    .agg([
        pl.col("len").filter(pl.col("CAGE") == True).sum().alias("true_len"),
        pl.col("len").sum().alias("total_len")
    ]).with_columns(
        (pl.col("true_len") / pl.col("total_len") * 100).alias("CAGE_percentage")
    ).write_csv("data/structural_category_cage_percentage.csv")