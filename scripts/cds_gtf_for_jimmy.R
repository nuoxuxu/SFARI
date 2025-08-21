#!/usr/bin/env Rscript
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(readr)

# Making GTF file for Jimmy's Ribo-seq analysis

annotation_gtf <- paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf")
predicted_cds_gtf <- "nextflow_results/V47/orfanage/orfanage.gtf"
sequence_feature <- "CDS"

SFARI_CDS <- rtracklayer::import(predicted_cds_gtf) %>%
    subset(type == sequence_feature) %>% 
    unique()

SFARI_CDS_nonunique <- rtracklayer::import(predicted_cds_gtf) %>%
    subset(type == sequence_feature)

gencode_CDS <- makeTxDbFromGFF(annotation_gtf) %>%
    cdsBy(by = "tx", use.names = TRUE) %>%
    unlist() %>%
    unique()

novel_exonic_regions_CDS <- GenomicRanges::subtract(SFARI_CDS, gencode_CDS) %>% unlist() %>% unique() %>% reduce()

hits <- findOverlaps(novel_exonic_regions_CDS, SFARI_CDS_nonunique)

novel_exonic_regions_CDS_expanded <- novel_exonic_regions_CDS[queryHits(hits)]

mcols(novel_exonic_regions_CDS_expanded)$transcript_id <- mcols(SFARI_CDS_nonunique[subjectHits(hits)])$transcript_id

novel_exonic_regions_CDS_expanded_grl <- split(novel_exonic_regions_CDS_expanded, novel_exonic_regions_CDS_expanded$transcript_id, drop=TRUE)

rtracklayer::export(unlist(novel_exonic_regions_CDS_expanded_grl), "export/variant/novel_exonic_regions_CDS_expanded_grl.gtf")