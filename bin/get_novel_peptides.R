#!/usr/bin/env Rscript
library(dplyr)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(readr)

args = commandArgs(trailingOnly=TRUE)

peptides_gtf <- args[1]
print(peptides_gtf)
reference_gtf <- args[2]
print(args[2])
genome_gff3_gtf <- args[3]
print(args[3])
genome_gff3_gtf_processed <- args[4]
print(args[4])
peptides_output <- args[5]
print(args[5])
genome_gff3_gtf_output <- args[6]
print(args[6])
peptide_SJ <- makeTxDbFromGFF("/gpfs/fs0/scratch/s/shreejoy/nxu/SFARI/work/39/9002e15ad40f6afc2395d311add4c6/SFARI_peptides_hybrid.gtf")
# Splice junctions
peptide_SJ <- makeTxDbFromGFF(peptides_gtf) %>%
    intronsByTranscript(use.names=TRUE)
peptide_SJ <- peptide_SJ[lengths(peptide_SJ) != 0] %>%
    unlist()

GENCODE_SJ <- makeTxDbFromGFF(reference_gtf) %>%
    intronsByTranscript(use.names=TRUE) %>%
    unlist()

peptide_SJ_not_in_GENCODE <- peptide_SJ[!(seq_along(peptide_SJ) %in% queryHits(findOverlaps(peptide_SJ, GENCODE_SJ, type = "equal")))]

SFARI_SJ <- makeTxDbFromGFF(genome_gff3_gtf) %>% 
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

GENCODE_exon <- makeTxDbFromGFF(reference_gtf) %>%
    exonsBy(by = "tx", use.names=TRUE)

peptide_exon_not_in_GENCODE <- peptide_exon[!(seq_along(peptide_exon) %in% queryHits(findSpliceOverlaps(peptide_exon, GENCODE_exon)))] %>% unlist()

SFARI_CDS <- makeTxDbFromGFF(genome_gff3_gtf) %>% 
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
SFARI_peptides[c("seqnames", "start", "end", "strand", "gene_name", "transcript_id", "type", "novelty")] %>%
    write.csv(peptides_output)

## Write genome_gff3_gtf.csv

genome_gff3_gtf <- read_csv(genome_gff3_gtf_processed)
genome_gff3_gtf <- genome_gff3_gtf %>%
    mutate(
        containing_novel_spl = case_when(
            transcript_id %in% c(pbid_containing_novel_SJs, pbid_containing_novel_exons) ~ TRUE,
            .default = FALSE
        )
    )

genome_gff3_gtf %>% write.csv(genome_gff3_gtf_output)
