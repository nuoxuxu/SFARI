import polars as pl
from src.utils import gtf_to_SJ, read_gtf, read_vcf
import os
from src.ryp import r, to_r

#--------------------------------Read the CSV with the defined column names--------------------------------
# ClinVar coordinates are 1-based
variants = read_vcf("data/clinvar_20250421.vcf", info=["CLNVC", "GENEINFO", "ALLELEID", "CLNSIG", "MC", "CLNDISDB", "CLNDN"])

variants = variants\
    .with_columns(
        pl.col("chrom").map_elements(lambda x: "".join(["chr", x]), return_dtype=pl.String)
        )\
    .filter(
        pl.col("chrom").str.contains("_").not_()
    )

#--------------------------------Read in final_transcripts and GENCODE SJ--------------------------------
final_transcripts_SJ = read_gtf("nextflow_results/V47/final_transcripts.gtf")\
    .filter(pl.col("feature") == "exon")\
    .pipe(gtf_to_SJ)\
    .unpivot(
        index=["chrom", "strand", "transcript_id"],
        variable_name="start_or_end",
        value_name="coord"
    )\
    .unique(["chrom", "start_or_end", "coord"]).drop("strand")

GENCODE_SJ = read_gtf("".join([os.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf"]))\
    .filter(pl.col("feature") == "exon")\
    .pipe(gtf_to_SJ)\
    .unpivot(
        index=["chrom", "strand", "transcript_id"],
        variable_name="start_or_end",
        value_name="coord"
    )\
    .unique(["chrom", "start_or_end", "coord"]).drop("strand")

#--------------------------------Get novel canonical splice sites--------------------------------
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

variants = variants\
    .filter(
        pl.col("pos").is_in(novel_splice_sites.filter(pl.col("coord").is_in(variants["pos"]))["coord"]) |   
        pl.col("pos").is_in(novel_splice_sites.filter(pl.col("coord_1").is_in(variants["pos"]))["coord_1"])
    )

variants = variants\
    .with_columns(
        pl.col("MC").str.split_exact("|", 1).struct.rename_fields(["SO", "MC"]).alias("MC")
    ).unnest("MC")\
    .with_columns(
        pl.col("MC").str.split_exact(",", 1).struct.rename_fields(["MC", "SO_2"]).alias("MC")
    ).unnest("MC")

variants.drop("info", "SO", "SO_2").write_csv("export/novel_canonical_splice_variants.csv")

#--------------------------------------------Visualization--------------------------------
to_r(variants, "variants")
r(
"""
library(ggplot2)
library(dplyr)

variants %>% 
    ggplot(aes(x=MC)) +
    geom_bar() +
    labs(x="Molecular consequences", y="Count") +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave("figures/test/canonical_splice_variants_MC.png", width=10, height=6)
"""
)

r(
"""
variants %>% 
    ggplot(aes(x=CLNSIG)) +
    geom_bar() +
    labs(x="Clinical significance", y="Count") +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave("figures/test/canonical_splice_variants_CLNSIG.png", width=10, height=6)
"""
)

#--------------------------------------------To-Do--------------------------------

# Are any of the genes SFARI genes?
# Look at the three examples of vatiants that contain genes that are autism related
# Disease ontology enrichment?