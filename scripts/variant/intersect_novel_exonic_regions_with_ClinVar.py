import polars as pl
from src.utils import read_vcf
from src.ryp import r, to_r, to_py

r(
"""
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(VariantAnnotation)

gencode_exons <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
    exonsBy(by = "tx", use.names = TRUE) %>% 
    unlist() %>% 
    unique()

SFARI_exons <- rtracklayer::import("nextflow_results/V47/final_transcripts.gtf") %>%
    subset(type == "exon") %>% 
    unique()

novel_exonic_regions <- GenomicRanges::subtract(SFARI_exons, gencode_exons, ignore.strand=TRUE) %>% unlist() %>% unique() %>% reduce(ignore.strand=TRUE)

clinvar_vcf <- readVcf("data/clinvar_20250421.vcf")
clinvar_vcf_gr <- rowRanges(clinvar_vcf)
seqlevelsStyle(clinvar_vcf_gr) <- "UCSC"

novel_exonic_regions_hits <- findOverlaps(clinvar_vcf_gr, novel_exonic_regions)
gencode_exons_hits <- findOverlaps(clinvar_vcf_gr, gencode_exons)
clinvar_vcf_gr <- clinvar_vcf_gr[setdiff(queryHits(novel_exonic_regions_hits), queryHits(gencode_exons_hits))]

id_list <- as.integer(names(clinvar_vcf_gr))
"""
)

id_list = to_py("id_list")

variants = read_vcf("data/clinvar_20250421.vcf", info=["CLNVC", "GENEINFO", "ALLELEID", "CLNDN", "CLNSIG", "MC"])

variants\
    .drop(["CLNVC", "GENEINFO", "ALLELEID", "CLNDN", "CLNSIG", "MC"])\
    .filter(
        pl.col("id").is_in(id_list)
    )\
    .with_columns(
        pl.col("chrom").map_elements(lambda x: "".join(["chr", x]), return_dtype=pl.String)
        )\
    .filter(
        pl.col("chrom").str.contains("_").not_()
    )\
    .write_csv("export/variant/novel_exonic_regions_clinvar.vcf", include_header=False, separator="\t", quote_style="never")
        
variants = variants\
    .filter(
        pl.col("id").is_in(id_list)
    )\
    .with_columns(
        pl.col("chrom").map_elements(lambda x: "".join(["chr", x]), return_dtype=pl.String)
        )\
    .filter(
        pl.col("chrom").str.contains("_").not_()
    )\
    .with_columns(
        pl.col("MC").str.split_exact("|", 1).struct.rename_fields(["SO", "MC"]).alias("MC")
    ).unnest("MC")\
    .with_columns(
        pl.col("MC").str.split_exact(",", 1).struct.rename_fields(["MC", "SO_2"]).alias("MC")
    ).unnest("MC").drop("SO", "SO_2")

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
ggsave("figures/test/novel_exonic_regions_MC.png", width=10, height=6)
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
ggsave("figures/test/novel_exonic_regions_CLNSIG.png", width=10, height=6)
"""
)