library(arrow)
library(dplyr)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(readr)

#-----------------------------------Load Datasets-----------------------------------#
annotation_gtf <- paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")
predicted_cds_gtf <- "nextflow_results/V47/orfanage/orfanage.gtf"
peptides_gtf <- "nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf"

#-----------------------------------Processing-----------------------------------#
# Splice junctions

peptide_SJ <- makeTxDbFromGFF(peptides_gtf) %>%
    intronsByTranscript(use.names=TRUE)
peptide_SJ <- peptide_SJ[lengths(peptide_SJ) != 0] %>%
    unlist()

GENCODE_SJ <- makeTxDbFromGFF(annotation_gtf) %>%
    intronsByTranscript(use.names=TRUE) %>%
    unlist()

SFARI_SJ <- makeTxDbFromGFF(predicted_cds_gtf) %>% 
    intronsByTranscript(use.names=TRUE)
SFARI_SJ <- SFARI_SJ[lengths(SFARI_SJ) != 0] %>%
    unlist()

peptide_SJ_in_GENCODE <- peptide_SJ[seq_along(peptide_SJ) %in% queryHits(findOverlaps(peptide_SJ, GENCODE_SJ, type = "equal"))]
peptide_SJ_not_in_GENCODE <- peptide_SJ[!(seq_along(peptide_SJ) %in% queryHits(findOverlaps(peptide_SJ, GENCODE_SJ, type = "equal")))]    

in_GENCODE_hits <- findOverlaps(peptide_SJ_in_GENCODE, SFARI_SJ, type="equal")
not_in_GENCODE_hits <- findOverlaps(peptide_SJ_not_in_GENCODE, SFARI_SJ, type="equal")

in_GENCODE_df <- tibble(
    peptide_SJ = names(peptide_SJ_in_GENCODE)[queryHits(in_GENCODE_hits)],
    SFARI_SJ = names(SFARI_SJ)[subjectHits(in_GENCODE_hits)],
    GENCODE = rep(TRUE, length(in_GENCODE_hits))
)

not_GENCODE_df <- tibble(
    peptide_SJ = names(peptide_SJ_not_in_GENCODE)[queryHits(hits)],
    SFARI_SJ = names(SFARI_SJ)[subjectHits(hits)],
    GENCODE = rep(FALSE, length(hits))
) 

bind_rows(in_GENCODE_df, not_GENCODE_df)

# Mono-exons

peptide_exon <- makeTxDbFromGFF(peptides_gtf) %>%
    exonsBy(use.names = TRUE, by = "tx")
peptide_exon <- peptide_exon[lengths(peptide_exon) == 1]

GENCODE_exon <- makeTxDbFromGFF(annotation_gtf) %>%
    exonsBy(by = "tx", use.names=TRUE)

SFARI_CDS <- makeTxDbFromGFF(predicted_cds_gtf) %>% 
    cdsBy(by = "tx", use.names=TRUE)

idx <- findSpliceOverlaps(peptide_exon[!(seq_along(peptide_exon) %in% queryHits(findSpliceOverlaps(peptide_exon, GENCODE_exon)))], SFARI_CDS) %>%
    subjectHits()

pbid_containing_novel_exons <- names(SFARI_CDS[idx])

#-----------------------------------Update classification-----------------------------------#
classification <- read_parquet("nextflow_results/V47/final_classification.parquet")

classification <- classification %>%
    mutate(
        containing_novel_spl = case_when(
            isoform %in% c(pbid_containing_novel_SJs, pbid_containing_novel_exons) ~ TRUE,
            .default = FALSE
        )
    )

classification %>% write_parquet("nextflow_results/V47/final_classification.parquet")