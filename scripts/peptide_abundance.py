import polars as pl
from src.utils import read_gtf

# Define functions
def remove_peptide_with_conflicted_mapping(peptide_mapping: pl.DataFrame) -> pl.DataFrame:
    """Remove peptides that map to transcripts that include known and novel transcripts."""
    peptide_to_remove = peptide_mapping\
        .group_by("peptide")\
        .agg(
            pl.count("pb"),
            pl.col("GENCODE").unique().len()
        )\
        .filter(pl.col("GENCODE")!=1)\
        .select("peptide")
        
    peptide_mapping = peptide_mapping\
        .filter(pl.col("peptide").is_in(peptide_to_remove["peptide"].to_numpy()).not_())

    return peptide_mapping

def get_logCPM_plus1(final_expression: pl.DataFrame) -> pl.DataFrame:
    """Get log2(CPM + 1) for each isoform."""
    mean_log_CPM_plus1 = final_expression\
        .with_columns([
            ((pl.col(col) / pl.col(col).sum() * 1_000_000) + 1).log(base=2).alias(col)
            for col in final_expression.columns[1:]
        ])\
        .with_columns(
            pl.mean_horizontal(
                pl.exclude("isoform")
            ).alias("mean_expression")
        )\
        .select(["isoform", "mean_expression"])
    return mean_log_CPM_plus1

# Read data
peptide_mapping = pl.read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet")
final_expression = pl.read_parquet("nextflow_results/V47/final_expression.parquet")
annot_peptides_hybrid = read_gtf("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf", attributes=["detected", "transcript_id"])

# mean_log_CPM_plus1 = get_logCPM_plus1(final_expression)
# peptide_mapping = remove_peptide_with_conflicted_mapping(peptide_mapping)
# peptide_detected = annot_peptides_hybrid\
#     .select("detected", "transcript_id")\
#     .unique()

peptide_mapping.unique("peptide")\
    .join(
        annot_peptides_hybrid.unique("transcript_id").select(["transcript_id", "detected"]),
        left_on="peptide",
        right_on="transcript_id",
        how="left"
    )\
    .filter(pl.col("detected").is_null())

peptide_mapping\
    .join(
        mean_log_CPM_plus1,
        left_on="pb",
        right_on="isoform",
        how="left"
    )\
    .group_by("peptide")\
    .agg(
        pl.col("GENCODE").first(),
        pl.col("mean_expression").mean(),
    )\
    .join(
        peptide_detected,
        left_on="peptide",
        right_on="transcript_id",
        how="left"
    )\
    .filter(pl.col("detected").is_null())\
    .write_parquet("nextflow_results/V47/orfanage/peptide_mapping_for_plotting.parquet")