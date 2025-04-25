import polars as pl
from src.utils import read_gtf

# Read in datasets
classification = pl.read_parquet("nextflow_results/V47/final_classification.parquet")
reftss = pl.read_csv("data/liftovered_mm39_to_hg38_peaks_overlapped_reftss_hg38_500bp.bed", separator="\t", has_header = False, new_columns=["chrom", "start", "end", "name", "score", "strand"])
polyA_site = pl.read_csv(
    "data/atlas.clusters.2.0.GRCh38.96.bed", 
    separator="\t", new_columns=["chrom", "start", "end", "name", "score", "strand"], 
    schema_overrides={"chrom": pl.String})\
    .with_columns(
        pl.col("chrom").map_elements(lambda x: "".join(["chr", x]), return_dtype=pl.String).alias("chrom")
    )
gtf = read_gtf("nextflow_results/V47/final_transcripts.gtf")

# Get list of tx validated by CAGE-seq peaks
validated_pbids = gtf\
    .filter(pl.col("feature")=="transcript")\
    .select(
        pl.col("seqname"),
        pl.col("transcript_id"),
        pos = pl.when(pl.col("strand")=="+")\
            .then(pl.col("start"))\
            .otherwise(pl.col("end"))
    )\
    .join_where(
        reftss,
        (pl.col("pos") >= (pl.col("start")-100)) &
        (pl.col("pos") <= (pl.col("end")+100))
    )\
    .filter(
        (pl.col("seqname") == pl.col("chrom"))
    )\
    .unique("transcript_id")\
    ["transcript_id"].to_list()

# Get percentage of isoforms validated by CAGE-seq peaks for visualization
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

# Annotate classification with CAGE peak information
classification = classification\
    .with_columns(
        within_CAGE_peak = pl.col("isoform").is_in(validated_pbids)
    )

# Get list of tx validated by polyA peaks
validated_pbids = gtf\
    .filter(pl.col("feature")=="transcript")\
    .select(
        pl.col("seqname"),
        pl.col("transcript_id"),
        pos = pl.when(pl.col("strand")=="+")\
            .then(pl.col("end"))\
            .otherwise(pl.col("start"))
    )\
    .join_where(
        polyA_site,
        (pl.col("pos") >= (pl.col("start"))) &
        (pl.col("pos") <= (pl.col("end")))
    )\
    .filter(
        pl.col("seqname")==pl.col("chrom")
    ).unique("transcript_id")\
    ["transcript_id"].to_list()

# Get percentage of isoforms validated by polyA peaks for visualization
classification\
    .with_columns(
        polyA_site = pl.col("isoform").is_in(validated_pbids),
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
    .group_by(["structural_category2", "polyA_site"])\
    .len()\
    .group_by("structural_category2")\
    .agg([
        pl.col("len").filter(pl.col("polyA_site") == True).sum().alias("true_len"),
        pl.col("len").sum().alias("total_len")
    ]).with_columns(
        (pl.col("true_len") / pl.col("total_len") * 100).alias("polyA_site_percentage")
    ).write_csv("data/structural_category_polyA_site_percentage.csv")