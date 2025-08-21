from src.utils import read_gtf
import polars as pl

riboseq = read_gtf("nextflow_results/V47/orfanage/UCSC_tracks/riboseq.gtf")

novel_CDS = read_gtf("nextflow_results/V47/orfanage/UCSC_tracks/novel_CDS.gtf")

classification = pl.read_parquet("nextflow_results/V47/final_classification.parquet")

peptides = read_gtf("nextflow_results/V47/orfanage/UCSC_tracks/detected_peptides.gtf")

plus_strand = novel_CDS\
    .filter(pl.col("strand") == "+")\
    .sort("seqname", "start")\
    .group_by("transcript_id", maintain_order=True)\
    .agg(
        pl.col("start").first(),
        pl.col("end").first(),
        pl.col("strand").first(),
        pl.col("seqname").first(),
    )\
    .with_columns(start = pl.col("start")-1)\
    .join(
        riboseq["seqname", "start", "strand", "source"],
        how="left",
        on=["seqname", "start", "strand"],
    )\
    .filter(pl.col("source").is_not_null())\
    .unique("start")\
    .select("transcript_id", "strand", "seqname", "start", "end")

minus_strand = novel_CDS\
    .filter(pl.col("strand") == "-")\
    .sort("seqname", "end", descending=True)\
    .group_by("transcript_id", maintain_order=True)\
    .agg(
        pl.col("start").first(),
        pl.col("end").first(),
        pl.col("strand").first(),
        pl.col("seqname").first(),
    )\
    .join(
        riboseq["seqname", "end", "strand", "source"],
        how="left",
        on=["seqname", "end", "strand"],
    )\
    .filter(pl.col("source").is_not_null())\
    .unique("end")\
    .select("transcript_id", "strand", "seqname", "start", "end")

validated_ncORFs = pl.concat([plus_strand, minus_strand])

validated_ncORFs.write_csv("validated_ncORFs.txt", separator="\t")

validated_ncORFs\
    .join(
        classification[["isoform", "associated_gene"]].rename({"isoform": "transcript_id"}),
        on="transcript_id",
        how="left"
    )