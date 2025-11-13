library(rtracklayer)
library(dplyr)
library(readr)
library(ggplot2)
library(arrow)

peptide_mapping <- read_parquet("nextflow_results/orfanage/peptide_mapping.parquet")
classification <- read_parquet("nextflow_results/V47/final_classification.parquet")
expression <- read_parquet("nextflow_results/V47/final_expression.parquet")

peptide_mapping <- peptide_mapping %>% 
    left_join(classification %>% select(isoform, structural_category), join_by(pb == isoform))

# Calculate mean expression in CPM from long-read data

lr_expression_cpm <- read_parquet("nextflow_results/V47/final_expression.parquet") %>% 
    dplyr::rename(
        NPC_1_2 = NPC_1_3,
        NPC_3_2 = NPC_3_3,
        CN_1_1 = CN_1_2,
        CN_1_2 = CN_1_3
    )%>%
    mutate(across( where(is.numeric), ~ .x / sum(.x) * 1e6 )) %>% 
    mutate(pb = isoform) %>% 
    # calculate mean across all numeric columns
    rowwise() %>%
    mutate(mean_cpm = mean(c_across(where(is.numeric)))) %>%
    ungroup() %>% 
    select(pb, mean_cpm)

# First definition

novel_detected_isoforms <- peptide_mapping %>% 
    filter(!GENCODE) %>% 
    distinct(pb) %>% 
    mutate(type="novel", detected=TRUE)

known_detected_isoforms <- peptide_mapping %>% 
    filter(structural_category=="full-splice_match") %>% 
    distinct(pb) %>% 
    mutate(type="known", detected=TRUE)

novel_undetected_isoforms <- classification %>% 
    filter(structural_category != "full-splice_match") %>% 
    filter(!(isoform %in% novel_detected_isoforms$pb)) %>% 
    distinct(isoform) %>% 
    mutate(pb=isoform, type="novel", detected=FALSE) %>%
    select(-isoform)

known_undetected_isoforms <- classification %>% 
    filter(structural_category == "full-splice_match") %>% 
    filter(!(isoform %in% known_detected_isoforms$pb)) %>% 
    distinct(isoform) %>% 
    mutate(pb=isoform, type="known", detected=FALSE) %>%
    select(-isoform)

combined <- bind_rows(
    novel_detected_isoforms,
    known_detected_isoforms,
    novel_undetected_isoforms,
    known_undetected_isoforms
)

# Second definition (less stringent)

novel_detected_isoforms <- peptide_mapping %>% 
    filter(structural_category!="full-splice_match") %>% 
    distinct(pb) %>% 
    mutate(type="novel", detected=TRUE)

known_detected_isoforms <- peptide_mapping %>% 
    filter(structural_category=="full-splice_match") %>% 
    distinct(pb) %>% 
    mutate(type="known", detected=TRUE)

novel_undetected_isoforms <- classification %>% 
    filter(structural_category != "full-splice_match") %>% 
    filter(!(isoform %in% novel_detected_isoforms$pb)) %>% 
    distinct(isoform) %>% 
    mutate(pb=isoform, type="novel", detected=FALSE) %>%
    select(-isoform)

known_undetected_isoforms <- classification %>% 
    filter(structural_category == "full-splice_match") %>% 
    filter(!(isoform %in% known_detected_isoforms$pb)) %>% 
    distinct(isoform) %>% 
    mutate(pb=isoform, type="known", detected=FALSE) %>%
    select(-isoform)

combined <- bind_rows(
    novel_detected_isoforms,
    known_detected_isoforms,
    novel_undetected_isoforms,
    known_undetected_isoforms
)

# Plotting

combined <- combined %>%
    left_join(lr_expression_cpm, join_by(pb == pb))

combined %>% 
    ggplot(aes(x=type, fill=detected)) +
    geom_bar(position="fill") +
    scale_y_continuous(labels=scales::percent_format()) +
    scale_fill_manual(values=c("TRUE"="#E69F00", "FALSE"="#56B4E9")) +
    labs(x="", y="Percentage of isoforms", fill="Detected by\npeptide mapping") +
    theme_minimal() +
    theme(text=element_text(size=16))
ggsave("figures/figure_2/transcript_pie_2.pdf", width=6, height=4)

combined %>% 
    ggplot(aes(x=type, fill=detected, y=mean_cpm)) +
    geom_boxplot(position=position_dodge(width=0.75), outlier.size=0.5) +
    scale_y_log10() +
    scale_fill_manual(values=c("TRUE"="#E69F00", "FALSE"="#56B4E9")) +
    labs(x="", y="Mean expression (CPM)", fill="Detected by\npeptide mapping") +
    theme_minimal() +
    theme(text=element_text(size=16))
ggsave("figures/figure_2/transcript_pie_expression_2.pdf", width=6, height=4)