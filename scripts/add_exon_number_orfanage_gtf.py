from src.utils import read_gtf
import polars as pl
from src.utils import collapse_isoforms_to_proteoforms
import polars.selectors as cs

orfanage_gtf = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")

# Add attribute "exon_number" to each exon

exons = orfanage_gtf\
    .filter(pl.col("feature") == "exon")\
    .with_columns(
        exon_number = pl.col("transcript_id").rank(method="ordinal").over("transcript_id")
    )

# Add attribute "exon_number" to each CDS, start_codon, and stop_codon

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

# Combine the two dataframes

orfanage_gtf_with_exon_num = pl.concat([exons, other_level2_features])\
    .sort(["transcript_id", "start"])\
    .with_columns(exon_number = pl.col("exon_number").cast(pl.String))\
    .with_columns(
        attributes = pl.col("attributes") + pl.lit(' exon_number "') + pl.col("exon_number") + pl.lit('";')
    ).drop("exon_number")

# Add protein_id

isoforms_to_proteoforms = collapse_isoforms_to_proteoforms(orfanage_gtf)

orfanage_gtf_with_exon_num_protein_id = orfanage_gtf_with_exon_num\
    .join(
        isoforms_to_proteoforms.rename({"isoform": "transcript_id"}),
        on = "transcript_id",
        how = "left"
    )\
    .rename({"base_isoform": "protein_id"})\
    .with_columns(
        attributes = pl.col("attributes") + pl.lit(' protein_id "') + pl.col("protein_id") + pl.lit('";')
    ).drop(["protein_id", "transcript_id"])

orfanage_gtf_with_exon_num_protein_id.write_csv("export/orfanage_with_exon_number.gtf", separator="\t", include_header=False, quote_style="never")