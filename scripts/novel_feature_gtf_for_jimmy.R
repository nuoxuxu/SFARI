#!/usr/bin/env Rscript
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(readr)

# This is for Jimmy's Ribo-seq stuff. I think he needs a gtf file that has, for each transcript, novel CDS regions (i.e. not in a GENCODE CDS).
# So this would contain overlapping CDS regions.

annotation_gtf <- "/project/rrg-shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf"
predicted_cds_gtf <- "nextflow_results/V47/orfanage/orfanage.gtf"
sequence_feature <- "three_prime_utr"
expressed_at_t30 <- read_csv("export/variant/expressed_t30.csv")

SFARI_feature <- rtracklayer::import(predicted_cds_gtf) %>%
    subset(type == sequence_feature) %>% 
    unique()

SFARI_feature_nonunique <- rtracklayer::import(predicted_cds_gtf) %>%
    subset(type == sequence_feature)

if (sequence_feature == "CDS") {
    gencode_feature <- makeTxDbFromGFF(annotation_gtf) %>%
        cdsBy(by = "tx", use.names = TRUE) %>%
        unlist() %>%
        unique()
} else if (sequence_feature == "three_prime_utr") {
    gencode_feature <- makeTxDbFromGFF(annotation_gtf) %>%
        threeUTRsByTranscript(use.names = TRUE) %>%
        unlist() %>%
        unique()
} else {
    stop("sequence_feature must be either CDS or three_prime_utr")
}

# Getting unique, non-overlapping novel features
novel_features_unique <- GenomicRanges::subtract(SFARI_feature, gencode_feature) %>% unlist() %>% unique() %>% reduce()

# Since we need each novel feature to be associated with a transcript, we will overlap back with the non-unique SFARI features
hits <- findOverlaps(novel_features_unique, SFARI_feature_nonunique)
novel_features_expanded <- novel_features_unique[queryHits(hits)]

# Label with transcript IDs
mcols(novel_features_expanded)$transcript_id <- mcols(SFARI_feature_nonunique[subjectHits(hits)])$transcript_id
mcols(novel_features_expanded)$type <- sequence_feature
# Split by transcript ID
# novel_features_expanded_grl <- split(novel_features_expanded, novel_features_expanded$transcript_id, drop=TRUE)

# Subset to only transcripts expressed at timepoint 30
subset(novel_features_expanded, mcols(novel_features_expanded)$transcript_id %in% pull(expressed_at_t30, x)) %>% 
    rtracklayer::export(sprintf("export/riboseq/novel_%s_expanded_grl_t30.gtf", sequence_feature))