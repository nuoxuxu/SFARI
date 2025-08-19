#!/usr/bin/env Rscript
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(tidyr)
library(readr)

# Read input files

args = commandArgs(trailingOnly=TRUE)
annotation_gtf <- args[1]
predicted_cds_gtf <- args[2]
novel_CDS <- args[3]

# Get novel UTRs

GENCODE_threeUTRs_gr <- makeTxDbFromGFF(annotation_gtf) %>%
    threeUTRsByTranscript(use.names = TRUE) %>%
    unlist() %>%
    unique()

GENCODE_fiveUTRs_gr <- makeTxDbFromGFF(annotation_gtf) %>%
    fiveUTRsByTranscript(use.names = TRUE) %>%
    unlist() %>%
    unique()

GENCODE_UTRs_gr <- c(GENCODE_threeUTRs_gr, GENCODE_fiveUTRs_gr)

SFARI_UTRs_gr <- rtracklayer::import(predicted_cds_gtf) %>%
    subset(type %in% c("five_prime_utr", "three_prime_utr")) %>%
    unique()

novel_UTRs_gr <- GenomicRanges::subtract(SFARI_UTRs_gr, GENCODE_UTRs_gr) %>% unlist() %>% unique() %>% reduce()

# Get novel CDSs

novel_CDS_gr <- rtracklayer::import(novel_CDS)

# Remove UTRs that have overlaps with novel CDSs

novel_UTRs_gr <- novel_UTRs_gr[-subjectHits(findOverlaps(novel_UTRs_gr, novel_CDS_gr)) %>% unique()]

# Formatting for UTRs

novel_UTRs_df <- as_tibble(novel_UTRs_gr) %>%
    mutate(type="novel")

known_UTRs_df <- as_tibble(reduce(GENCODE_UTRs_gr)) %>%
    mutate(type="known")

combined_UTRs_df <- bind_rows(novel_UTRs_df, known_UTRs_df)

# Formatting for CDS

GENCODE_CDS_gr <- makeTxDbFromGFF(annotation_gtf) %>%
    cdsBy("tx", use.names = TRUE) %>%
    unlist() %>%
    unique() %>%
    reduce()

known_CDS_df <- as_tibble(GENCODE_CDS_gr) %>%
    select(c(seqnames, start, end, width, strand)) %>%
    mutate(type = "known")

novel_CDS_df <- as_tibble(novel_CDS_gr) %>%
    select(c(seqnames, start, end, width, strand)) %>%
    mutate(type = "novel")

combined_CDS_df <- bind_rows(known_CDS_df, novel_CDS_df)

# Write files

write_csv(combined_UTRs_df, "UTR_regions.csv")
write_csv(combined_CDS_df, "CDS_regions.csv")