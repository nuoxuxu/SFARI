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
    peptide = names(peptide_SJ_in_GENCODE)[queryHits(in_GENCODE_hits)],
    pb = names(SFARI_SJ)[subjectHits(in_GENCODE_hits)],
    GENCODE = rep(TRUE, length(in_GENCODE_hits)),
    type = rep("splice-junction", length(in_GENCODE_hits))
)

not_GENCODE_df <- tibble(
    peptide = names(peptide_SJ_not_in_GENCODE)[queryHits(not_in_GENCODE_hits)],
    pb = names(SFARI_SJ)[subjectHits(not_in_GENCODE_hits)],
    GENCODE = rep(FALSE, length(not_in_GENCODE_hits)),
    type = rep("splice-junction", length(not_in_GENCODE_hits))
) 

peptide_SJ_mapping <- bind_rows(in_GENCODE_df, not_GENCODE_df)

# Mono-exons

peptide_exon <- makeTxDbFromGFF(peptides_gtf) %>%
    exonsBy(use.names = TRUE, by = "tx")
peptide_exon <- peptide_exon[lengths(peptide_exon) == 1]

GENCODE_exon <- makeTxDbFromGFF(annotation_gtf) %>%
    exonsBy(by = "tx", use.names=TRUE)

SFARI_CDS <- makeTxDbFromGFF(predicted_cds_gtf) %>% 
    cdsBy(by = "tx", use.names=TRUE)

peptide_exon_in_GENCODE <- peptide_exon[seq_along(peptide_exon) %in% queryHits(findSpliceOverlaps(peptide_exon, GENCODE_exon))]
peptide_exon_not_in_GENCODE <- peptide_exon[!(seq_along(peptide_exon) %in% queryHits(findSpliceOverlaps(peptide_exon, GENCODE_exon)))]  

in_GENCODE_hits <- findSpliceOverlaps(peptide_exon_in_GENCODE, SFARI_CDS, type="equal")
not_in_GENCODE_hits <- findSpliceOverlaps(peptide_exon_not_in_GENCODE, SFARI_CDS, type="equal")

in_GENCODE_df <- tibble(
    peptide = names(peptide_exon_in_GENCODE)[queryHits(in_GENCODE_hits)],
    pb = names(SFARI_CDS)[subjectHits(in_GENCODE_hits)],
    GENCODE = rep(TRUE, length(in_GENCODE_hits)),
    type = "mono-exonic"

)

not_GENCODE_df <- tibble(
    peptide = names(peptide_exon_not_in_GENCODE)[queryHits(not_in_GENCODE_hits)],
    pb = names(SFARI_CDS)[subjectHits(not_in_GENCODE_hits)],
    GENCODE = rep(FALSE, length(not_in_GENCODE_hits)),
    type = "mono-exonic"
)

peptide_exon_mapping <- bind_rows(in_GENCODE_df, not_GENCODE_df)

# Combine two types of peptides

out <- bind_rows(peptide_SJ_mapping, peptide_exon_mapping)
out %>% write_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet")
