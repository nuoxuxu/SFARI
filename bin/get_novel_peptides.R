library(dplyr)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(readr)

peptide_gtf <- "all_peptides.gtf"

# Splice junctions
peptide_SJ <- makeTxDbFromGFF(peptide_gtf) %>%
    intronsByTranscript(use.names=TRUE)
peptide_SJ <- peptide_SJ[lengths(peptide_SJ) != 0] %>%
    unlist()

GENCODE_SJ <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v39.annotation.gtf")) %>%
    intronsByTranscript(use.names=TRUE) %>%
    unlist()

peptide_SJ_not_in_GENCODE <- peptide_SJ[!(seq_along(peptide_SJ) %in% queryHits(findOverlaps(peptide_SJ, GENCODE_SJ, type = "equal")))]

SFARI_SJ <- makeTxDbFromGFF("nextflow_results/transcripts_filtered.fasta.transdecoder.genome.gff3.gtf") %>% 
    intronsByTranscript(use.names=TRUE)
SFARI_SJ <- SFARI_SJ[lengths(SFARI_SJ) != 0] %>%
    unlist()

pbid_containing_novel_SJs <- SFARI_SJ[unique(subjectHits(findOverlaps(peptide_SJ_not_in_GENCODE, SFARI_SJ, type="equal")))] %>%
    names() %>%
    unique()

# Mono-exons
peptide_exon <- makeTxDbFromGFF(peptide_gtf) %>%
    exonsBy(use.names = TRUE, by = "tx")
peptide_exon <- peptide_exon[lengths(peptide_exon) == 1]

GENCODE_exon <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v39.annotation.gtf")) %>%
    exonsBy(by = "tx", use.names=TRUE)

peptide_exon_not_in_GENCODE <- peptide_exon[!(seq_along(peptide_exon) %in% queryHits(findSpliceOverlaps(peptide_exon, GENCODE_exon)))] %>% unlist()

SFARI_CDS <- makeTxDbFromGFF("nextflow_results/transcripts_filtered.fasta.transdecoder.genome.gff3.gtf") %>% 
    cdsBy(by = "tx", use.names=TRUE)

idx <- findSpliceOverlaps(peptide_exon[!(seq_along(peptide_exon) %in% queryHits(findSpliceOverlaps(peptide_exon, GENCODE_exon)))], SFARI_CDS) %>%
    subjectHits()

pbid_containing_novel_exons <- names(SFARI_CDS[idx])

# Formating

## Write SFARI_peptides

SFARI_peptides <- as.data.frame(import.gff(peptide_gtf))
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
    write.csv("transcript_vis_app/SFARI_peptides.csv")

## Write genome_gff3_gtf.csv

genome_gff3_gtf <- read_csv("transcript_vis_app/genome_gff3_gtf.csv")
genome_gff3_gtf <- genome_gff3_gtf %>%
    mutate(
        containing_novel_spl = case_when(
            transcript_id %in% c(pbid_containing_novel_SJs, pbid_containing_novel_exons) ~ TRUE,
            .default = FALSE
        )
    )

genome_gff3_gtf %>% write.csv("transcript_vis_app/genome_gff3_gtf.csv")
