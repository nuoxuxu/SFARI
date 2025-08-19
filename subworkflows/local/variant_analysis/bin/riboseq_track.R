#!/usr/bin/env Rscript
library(readxl)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(stringr)
library(readr)

args <- commandArgs(trailingOnly=TRUE)

riboseq_file <- args[1]

riboseq <- read_excel(riboseq_file, sheet = 2) %>%
    dplyr::rename(
        seqname = chrom,
        start = codon5,
        end = codon3,
        transcript_id = orfID,
    ) %>%
    mutate(
        source = "riboseq",
        feature = "CDS",
        score = ".",
        frame = "."
    ) %>% 
    select(seqname, source, feature, start, end, score, strand, frame, transcript_id)

riboseq %>%
    mutate(gene_id = str_extract(transcript_id, "\\blnckb\\.\\d+(?!\\.\\d)")) %>%
    mutate(transcript_id = str_extract(transcript_id, "\\blnckb\\.\\d+\\.\\d+\\b")) %>%
    mutate(
        attributes = paste0(
            'gene_id \"', gene_id, '\"; ',
            'transcript_id \"', transcript_id, '\";'
        )
    ) %>% 
    select(-c(gene_id, transcript_id)) %>%
    write_tsv("riboseq.gtf", col_names = FALSE, quote="none", escape="none")
