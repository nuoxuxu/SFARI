library(GenomicRanges)
library(GenomicFeatures)
library(dplyr)
library(readr)
library(VariantAnnotation)
library(rtracklayer)
library(stringr)
library(purrr)

clinvar_vcf <- readVcf("data/clinvar_20250421.vcf")
clinvar_vcf_gr <- rowRanges(clinvar_vcf)
seqlevelsStyle(clinvar_vcf_gr) <- "UCSC"

orfanage_txdb <- makeTxDbFromGFF("nextflow_results/V47/orfanage/orfanage.gtf")
gencode_txdb <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.sorted.gtf"))

orfanage_splice_site_variants <- locateVariants(clinvar_vcf_gr, orfanage_txdb, SpliceSiteVariants())
orfanage_splice_site_variants <- orfanage_splice_site_variants[!duplicated(names(orfanage_splice_site_variants))]
gencode_splice_site_variants <- locateVariants(clinvar_vcf_gr, gencode_txdb, SpliceSiteVariants())
gencode_splice_site_variants <- gencode_splice_site_variants[!duplicated(names(gencode_splice_site_variants))]

novel_splice_site_variants <- orfanage_splice_site_variants[names(orfanage_splice_site_variants) %in% setdiff(names(orfanage_splice_site_variants), names(gencode_splice_site_variants))]
novel_splice_site_variants <- novel_splice_site_variants[width(novel_splice_site_variants)==1]

vcf_df <- read_tsv("data/clinvar_20250421.vcf", comment="##")

vcf_df %>% 
    dplyr::filter(ID %in% names(novel_splice_site_variants)) %>% 
    write_tsv("export/variant/novel_splice_site_ClinVar.vcf")