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

gencode_CDS <- makeTxDbFromGFF(annotation_gtf) %>%
    cdsBy(by = "tx", use.names = TRUE) %>%
    unlist() %>%
    unique()

SFARI_CDS <- rtracklayer::import(predicted_cds_gtf) %>%
    subset(type == "CDS") %>%
    unique()

novel_CDS <- GenomicRanges::subtract(SFARI_CDS, gencode_CDS) %>% unlist() %>% unique() %>% reduce()

# In order to show novel exonic regions in UCSC genome browser, we need to add gene_id and transcript_id and
# change the feature name from rtracklayer::export default "sequence_feature" to "exon".

hits <- findOverlaps(novel_CDS, SFARI_CDS)

mcols(novel_CDS)$gene_id <- mcols(SFARI_CDS[subjectHits(hits[!duplicated(queryHits(hits))])])$gene_id

mcols(novel_CDS)$transcript_id <- mcols(SFARI_CDS[subjectHits(hits[!duplicated(queryHits(hits))])])$transcript_id

as_tibble(novel_CDS) %>%
    mutate(
        source = "ORFanage",
        feature = "CDS",
        score = ".",
        frame = ".",
        attributes = paste0(
            'gene_id \"', gene_id, '\"; ',
            'transcript_id \"', transcript_id, '\";'
        )
    ) %>% 
    dplyr::select(-c(gene_id, transcript_id)) %>% 
    dplyr::select(seqnames, source, feature, start, end, score, strand, frame, attributes) %>% 
    write_tsv("novel_CDS.gtf", col_names = FALSE, quote="none", escape="none")
