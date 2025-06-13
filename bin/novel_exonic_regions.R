#!/usr/bin/env Rscript
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(readr)

args <- commandArgs(trailingOnly=TRUE)

annotation_gtf <- args[1]
predicted_cds_gtf <- args[2]

gencode_exons <- makeTxDbFromGFF(annotation_gtf) %>%
    exonsBy(by = "tx", use.names = TRUE) %>% 
    unlist() %>%
    unique()

SFARI_exons <- rtracklayer::import(predicted_cds_gtf) %>%
    subset(type == "exon") %>%
    unique()

novel_exonic_regions <- GenomicRanges::subtract(SFARI_exons, gencode_exons) %>% unlist() %>% unique() %>% reduce()

# In order to show novel exonic regions in UCSC genome browser, we need to add gene_id and transcript_id and
# change the feature name from rtracklayer::export default "sequence_feature" to "exon".

hits <- findOverlaps(novel_exonic_regions, SFARI_exons)

mcols(novel_exonic_regions)$gene_id <- mcols(SFARI_exons[subjectHits(hits[!duplicated(queryHits(hits))])])$gene_id

mcols(novel_exonic_regions)$transcript_id <- mcols(SFARI_exons[subjectHits(hits[!duplicated(queryHits(hits))])])$transcript_id

rtracklayer::export(novel_exonic_regions, "novel_exonic_regions.gtf")

system("sed -i 's/sequence_feature/exon/g' novel_exonic_regions.gtf")