library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(VariantAnnotation)
library(tidyr)
library(parallel)

source("src/utils.R")

novel_exonic_regions <- readRDS("export/variant/novel_exonic_regions.rds")

gencode_exons <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
    exonsBy(by = "tx", use.names = TRUE) %>% 
    unlist() %>% 
    unique()

clinvar_vcf <- readVcf("data/clinvar_20250421.vcf")
clinvar_vcf_gr <- rowRanges(clinvar_vcf)
seqlevelsStyle(clinvar_vcf_gr) <- "UCSC"

novel_exonic_regions_hits <- findOverlaps(clinvar_vcf_gr, novel_exonic_regions)
gencode_exons_hits <- findOverlaps(clinvar_vcf_gr, gencode_exons)
clinvar_vcf_gr <- clinvar_vcf_gr[setdiff(queryHits(novel_exonic_regions_hits), queryHits(gencode_exons_hits))]
clinvar_vcf_gr <- clinvar_vcf_gr[width(clinvar_vcf_gr) == 1]

id_list <- as.integer(names(clinvar_vcf_gr))
ClinVar_df <- read_vcf("data/clinvar_20250421.vcf")
ClinVar_df <- ClinVar_df %>% 
    filter(id %in% id_list)

ClinVar_df %>% dplyr::select(-c(CLNVC, GENEINFO)) %>% 
    write_tsv("export/variant/novel_exonic_regions_ClinVar.vcf")