library(arrow)
library(dplyr)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(readr)
library(ggplot2)
library(patchwork)

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

#-----------------------------------Plotting-----------------------------------#

# Add "detected" column to the peptide mapping

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet")

peptides_gtf <- import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>%
    as_tibble() %>% 
    distinct(transcript_id, .keep_all = TRUE)

classification <- read_tsv("nextflow_results/V47/orfanage/SFARI.protein_classification.tsv")    

peptide_mapping <- peptide_mapping %>% 
    left_join(
        peptides_gtf %>% 
            dplyr::select(transcript_id, detected),
        join_by(peptide == transcript_id)
    ) %>% 
    left_join(
        classification %>% dplyr::select(pb, protein_classification_base), 
        join_by(pb == pb)
    )

all <- peptide_mapping %>% 
    filter(!GENCODE) %>% 
    filter(protein_classification_base == "pNNC") %>%
    distinct(pb) %>% 
    nrow()

detected <- peptide_mapping %>% 
    filter(!GENCODE) %>% 
    filter(protein_classification_base == "pNNC") %>%
    filter(detected == "True") %>%
    distinct(pb) %>% 
    nrow()

undetected <- all - detected

p1 <- tibble(
    is_detected = c(TRUE, FALSE),
    len = c(detected, undetected)
    ) %>% 
    mutate(
        percent = len / sum(len) * 100, # calculate percentage over all groups
        type = "known" # add literal column type
    ) %>% 
    arrange(desc(is_detected)) %>%
    mutate(ypos = cumsum(percent) - 0.5 * percent) %>% 
    ggplot(aes(x="", y=percent, fill=is_detected)) +
    geom_bar(width=1, stat="identity", color = "white") +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 9) +
    theme(
        legend.position = "none"
    )    

# Plot percentage of transcripts that are validated by novel peptides only

all <- peptide_mapping %>% 
    filter(GENCODE == FALSE) %>%
    filter(protein_classification_base == "pNIC") %>%
    distinct(pb) %>% 
    nrow()

detected <- peptide_mapping %>% 
    filter(GENCODE == FALSE) %>%
    filter(protein_classification_base == "pNIC") %>%
    filter(detected == "True") %>%
    distinct(pb) %>% 
    nrow()

undetected <- all - detected

p2 <- tibble(
    is_detected = c(TRUE, FALSE),
    len = c(detected, undetected)
    ) %>% 
    mutate(
        percent = len / sum(len) * 100, # calculate percentage over all groups
        type = "known" # add literal column type
    ) %>% 
    arrange(desc(is_detected)) %>%
    mutate(ypos = cumsum(percent) - 0.5 * percent) %>% 
    ggplot(aes(x="", y=percent, fill=is_detected)) +
    geom_bar(width=1, stat="identity", color = "white") +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 9)

p1 + p2

ggsave("figures/figure_2/validation_by_novel_peptides.pdf", width = 10, height = 5)