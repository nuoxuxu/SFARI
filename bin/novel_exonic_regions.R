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
sequence_feature <- args[3]
output <- args[4]

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

as_tibble(novel_exonic_regions) %>% 
    mutate(
        source = "ORFanage",
        feature = "exon",
        score = ".",
        frame = ".",
        attributes = paste0(
            'gene_id \"', gene_id, '\"; ',
            'transcript_id \"', transcript_id, '\";'
        )
    ) %>%
    dplyr::select(-c(gene_id, transcript_id)) %>% 
    dplyr::select(seqnames, source, feature, start, end, score, strand, frame, attributes) %>% 
    write_tsv(output, col_names = FALSE, quote="none", escape="none")