#!/usr/bin/env Rscript
library(dplyr)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(readr)

args <- commandArgs(trailingOnly=TRUE)

annotation_gtf <- args[1]
predicted_cds_gtf <- args[2]
peptides_gtf <- args[3]
peptides_output <- args[4]

# Splice junctions

peptide_SJ <- makeTxDbFromGFF(peptides_gtf) %>%
    intronsByTranscript(use.names=TRUE)
peptide_SJ <- peptide_SJ[lengths(peptide_SJ) != 0] %>%
    unlist()

GENCODE_SJ <- makeTxDbFromGFF(annotation_gtf) %>%
    intronsByTranscript(use.names=TRUE) %>%
    unlist()

peptide_SJ_not_in_GENCODE <- peptide_SJ[!(seq_along(peptide_SJ) %in% queryHits(findOverlaps(peptide_SJ, GENCODE_SJ, type = "equal")))]

SFARI_SJ <- makeTxDbFromGFF(predicted_cds_gtf) %>% 
    intronsByTranscript(use.names=TRUE)
SFARI_SJ <- SFARI_SJ[lengths(SFARI_SJ) != 0] %>%
    unlist()

pbid_containing_novel_SJs <- SFARI_SJ[unique(subjectHits(findOverlaps(peptide_SJ_not_in_GENCODE, SFARI_SJ, type="equal")))] %>%
    names() %>%
    unique()

# Mono-exons

peptide_exon <- makeTxDbFromGFF(peptides_gtf) %>%
    exonsBy(use.names = TRUE, by = "tx")
peptide_exon <- peptide_exon[lengths(peptide_exon) == 1]

GENCODE_exon <- makeTxDbFromGFF(annotation_gtf) %>%
    exonsBy(by = "tx", use.names=TRUE)

peptide_exon_not_in_GENCODE <- peptide_exon[!(seq_along(peptide_exon) %in% queryHits(findSpliceOverlaps(peptide_exon, GENCODE_exon)))] %>% unlist()

SFARI_CDS <- makeTxDbFromGFF(predicted_cds_gtf) %>% 
    cdsBy(by = "tx", use.names=TRUE)

idx <- findSpliceOverlaps(peptide_exon[!(seq_along(peptide_exon) %in% queryHits(findSpliceOverlaps(peptide_exon, GENCODE_exon)))], SFARI_CDS) %>%
    subjectHits()

pbid_containing_novel_exons <- names(SFARI_CDS[idx])

# Formating

## Write SFARI_peptides

SFARI_peptides <- as.data.frame(import.gff(peptides_gtf))
SFARI_peptides <- SFARI_peptides %>%
    mutate(
        type = case_when(
            transcript_id %in% names(peptide_exon) ~ "mono-exonic",
            transcript_id %in% names(peptide_SJ) ~ "splice-junction"
        ),
        novelty = case_when(
            transcript_id %in% c(names(peptide_exon_not_in_GENCODE), names(peptide_SJ_not_in_GENCODE)) ~ "novel",
            .default = "known"
        )
    )

SFARI_peptides %>%
    mutate(
        feature = "exon",
        attribute = paste0('gene_id "', gene_id, '"; transcript_id "', transcript_id, '"; gene_name "', gene_name, '"; detected "', detected, '"; type "', type, '"; novelty "', novelty, '";')) %>%
        dplyr::select(c(seqnames, source, feature, start, end, score, strand, phase, attribute)) %>%
        write.table(peptides_output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)