library(arrow)
library(dplyr)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(readr)
library(ggplot2)
library(patchwork)

#-----------------------------------Prepare data-----------------------------------#

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

#--------------------------------Second version-----------------------------------#

p1 <- peptide_mapping %>% 
    distinct(peptide, .keep_all = TRUE) %>% 
    filter(GENCODE == TRUE) %>% 
    group_by(detected) %>% 
    summarise(n = n()) %>% 
    mutate(percent = n / sum(n) * 100) %>% 
    arrange(desc(detected)) %>%
    mutate(ypos = cumsum(percent) - 0.5 * percent) %>% 
    ggplot(aes(x="", y=percent, fill=detected)) +
    geom_bar(width=1, stat="identity", color = "white") +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 9)

p2 <- peptide_mapping %>% 
    distinct(peptide, .keep_all = TRUE) %>% 
    filter(GENCODE == FALSE) %>% 
    group_by(detected) %>% 
    summarise(n = n()) %>% 
    mutate(percent = n / sum(n) * 100) %>% 
    arrange(desc(detected)) %>%
    mutate(ypos = cumsum(percent) - 0.5 * percent) %>% 
    ggplot(aes(x="", y=percent, fill=detected)) +
    geom_bar(width=1, stat="identity", color = "white") +
    coord_polar("y", start=0) +
    theme_void() +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 9)       

p1 + p2

ggsave("figures/figure_2/validation_by_novel_peptides.pdf", width = 10, height = 5)
