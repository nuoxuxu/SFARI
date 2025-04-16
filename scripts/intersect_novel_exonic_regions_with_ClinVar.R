library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)
library(VariantAnnotation)

gencode_exons <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
    exonsBy(by = "tx", use.names = TRUE) %>% 
    unlist() %>% 
    unique()

SFARI_exons <- rtracklayer::import("nextflow_results/V47/final_transcripts.gtf") %>%
    subset(type == "exon") %>% 
    unique()

novel_genomic_regions <- GenomicRanges::setdiff(SFARI_exons, gencode_exons)

clinvar_vcf <- readVcf("data/clinvar_20250409.vcf")

clinvar_vcf_gr <- rowRanges(clinvar_vcf)

seqlevelsStyle(clinvar_vcf_gr) <- "UCSC"

GenomicRanges::intersect(novel_genomic_regions, clinvar_vcf_gr, ignore.strand=TRUE)
