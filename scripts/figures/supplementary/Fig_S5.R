library(patchwork)
library(ggplot2)
library(dplyr)
library(arrow)
library(edgeR)
library(tidyr)
library(readxl)
library(ggtext)
library(ggpubr)
library(rtracklayer)
library(readr)
library(purrr)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 14)        
    )
theme_set(my_theme)

colorVector <- c("TRUE" = "#00BFC4", "FALSE" = "#F8766D")

# SJ validation
LR_SJ_plot <- read_parquet("nextflow_results/V47/LR_SJ_2.parquet") %>% 
    mutate(SR = replace_na(SR, FALSE)) %>% 
    mutate(SR = factor(SR, levels=c(TRUE, FALSE))) %>%
    mutate(type = recode(type, "known" = "GENCODE", "novel" = "Novel", "shared" = "Shared")) %>%
    ggplot(aes(x=type, y=mean_log2_cpm, fill=SR)) +
    geom_violin() +
    scale_fill_manual("Supproted by\nshort reads", values=colorVector) +
    stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
    stat_compare_means(
        aes(label = after_stat(p.signif)),
        comparisons = list(c("GENCODE", "Novel"), c("Novel", "Shared"), c("GENCODE", "Shared")),
        method = "t.test",
        label.y = c(12, 13, 14)
    ) +
    labs(x = "Splice junctions", y = expression("Mean"~log[2]*"(CPM + 1) of isoforms containing the SJ across all samples"))
# Patowary validation
LR_patowary_plot <- read_parquet("nextflow_results/V47/LR_patowary.parquet") %>% 
    mutate(supported = factor(supported, levels=c(TRUE, FALSE))) %>%
    mutate(type = recode(type, "known" = "GENCODE", "novel" = "Novel")) %>%
    ggplot(aes(x=type, y=mean_log2_cpm, fill=supported)) +
    geom_boxplot() +
    scale_fill_manual("Supported by\nPatowary et al.", values=colorVector) +
    stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
    stat_compare_means(
        aes(label = after_stat(p.signif)),
        comparisons = list(c("GENCODE", "Novel")),
        method = "t.test",
        label.y = c(15)
    ) +
    labs(x = "Full-length isoforms", y = expression(atop("Mean"~log[2]*"(CPM + 1)", "across all samples")))

# Peptide validation
peptide_mapping <- read_parquet("nextflow_results/orfanage/peptide_mapping.parquet")
novel_peptides <- read_csv("nextflow_results/orfanage/novel_peptides.csv")
protein_classification <- read_csv("nextflow_results/V47/orfanage/isoforms_to_proteoforms.csv")
classification <- read_parquet("nextflow_results/V47/final_classification.parquet")
expression <- read_parquet("nextflow_results/V47/final_expression.parquet")
dge <- DGEList(counts=as.matrix(expression[, grep("_", colnames(expression))]), genes=expression[, "isoform"])
log2_cpm <- cpm(dge, log=TRUE)
mean_log2_cpm <- rowMeans(log2_cpm)
mean_log2_cpm <- as_tibble(data.frame(isoform=expression$isoform, mean_log2_cpm=mean_log2_cpm))

bin_size <- 10
detected_known_isoforms <- peptide_mapping %>% 
    filter(grepl("^ENST", transcript_id)) %>% 
    left_join(classification[, c("isoform", "associated_transcript")], join_by(transcript_id == associated_transcript)) %>% 
    distinct(isoform) %>% 
    drop_na() %>% 
    pull(isoform)

detected_novel_isoforms <- novel_peptides %>% 
    left_join(protein_classification, join_by(transcript_id == base_isoform)) %>%
    distinct(isoform) %>%
    pull(isoform)

known_only <- classification %>% 
    filter(structural_category == "full-splice_match") %>%
    select(c(isoform, structural_category)) %>%
    left_join(mean_log2_cpm, join_by(isoform == isoform)) %>% 
    mutate(
        detected = isoform %in% detected_known_isoforms
    ) %>% 
    mutate(type="known")

novel_only <- classification %>% 
    filter(structural_category != "full-splice_match") %>%
    left_join(mean_log2_cpm, join_by(isoform == isoform)) %>% 
    mutate(
        detected = isoform %in% detected_novel_isoforms
    ) %>% 
    mutate(type="novel")

combined <- bind_rows(
        known_only,
        novel_only
    ) %>% 
    mutate(detected = TRUE)

peptide_plot <- combined %>% 
    mutate(type = recode(type, "known" = "GENCODE", "novel" = "Novel")) %>%
    ggplot(aes(x=type, y=mean_log2_cpm, fill=detected)) +
    scale_fill_manual("Supported by\npeptides", values=colorVector) +
    geom_boxplot(position=position_dodge(width=0.75), outlier.size=0.5) +
    stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test", comparisons=list(c("GENCODE", "Novel"))) +
    labs(x="Full-length isoforms", y=expression(atop("Mean"~log[2]*"(CPM + 1)", "across all samples")))

# Ribo-seq data
file_path <- "data/katherine/Duffy_TE_filtered.xlsx" 

ribo_seq <- file_path %>% 
    excel_sheets() %>%               # 1. Get the names of the sheets
    set_names() %>%                  # 2. Name the vector (crucial for step 4)
    map(~ read_excel(path = file_path, sheet = .x)) %>% # 3. Read each sheet
    bind_rows(.id = "sheet_name") %>% 
    mutate(
    annotation = case_when(
        grepl("GENCODE", sheet_name) ~ "GENCODE",
        grepl("Novel", sheet_name) ~ "Novel"),
    regions = case_when(
        grepl("UTR", sheet_name) ~ "UTR",
        grepl("CDS", sheet_name) ~ "ORF")
    )

ribo_seq_plot <- ribo_seq %>% 
    ggplot(aes(x=annotation, y=`AvgTE (log2)`, fill=regions)) +
    geom_boxplot() +
    scale_fill_manual("Regions", values=c("ORF"="#6FA7D7", "UTR"="#DEDCDC")) +
    labs(y=expression("Relative tanslation efficiency ("*log[2]*")")) +
    stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test") +
    stat_compare_means(
        aes(label = after_stat(p.signif)),
        comparisons = list(c("GENCODE", "Novel")),
        method = "t.test",
        label.y = c(11)
    )

# Combining plots
(LR_patowary_plot + LR_SJ_plot) / (peptide_plot + ribo_seq_plot) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
ggsave("figures/supplementary/expression_whether_novel.pdf", width = 14, height = 8)