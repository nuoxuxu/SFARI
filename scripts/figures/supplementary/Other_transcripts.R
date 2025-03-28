library(dplyr)
library(ggplot2)
library(arrow)
library(scales)
library(patchwork)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

# Transcript classification for the transcripts in the Other category

transcript_classification_hist <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    filter(structural_category %in% c("fusion", "genic", "antisense", "intergenic", "moreJunctions")) %>% 
    group_by(structural_category) %>%
    summarise(len = n()) %>%
    mutate(
        percentage = (len / sum(len)) * 100,
        structural_category = factor(structural_category)
    ) %>%
    ggplot(
        aes(x = structural_category, y = len, fill = structural_category)
        ) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(percentage, 1), "%")), colour = "black", size = 4, vjust=-0.5) +
        scale_y_continuous(labels = function(x) x / 1000, limits=c(0, 5000)) +
        labs(x = NULL, y = expression("Transcripts (x" ~ 10^3 * ")"))
ggsave("figures/supplementary/other_transcript_classification_hist.pdf", width = 10, height = 7)

# Counts vs abundance

row_sum <- read_parquet("nextflow_results/V47/final_expression.parquet") %>%
    mutate(
        row_sum = rowSums(across(where(is.numeric)))
    ) %>%
    dplyr::select(row_sum)

count_vs_abundance <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    bind_cols(row_sum) %>%
    filter(
        structural_category %in% c("fusion", "genic", "antisense", "intergenic", "moreJunctions"),
        row_sum > 10,
        ) %>%
    ggplot(aes(x = row_sum, fill = structural_category)) +
    geom_histogram(position = position_fill(), alpha = .75) +
    scale_x_log10() +
    labs(x = "Total observed counts", y = "Transcript proportion")

# Number of exons vs abudance

nexon_vs_abundance <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    filter(structural_category %in% c("fusion", "genic", "antisense", "intergenic", "moreJunctions")) %>% 
    ggplot(aes(x = exons, fill = structural_category)) +
    geom_histogram(alpha = .75, binwidth = 1) +
    scale_x_continuous(limits = c(1, 40)) +
    scale_y_continuous(labels = function(x) x / 1000) +
    labs(
        x = "# Exons",
        y = expression("Transcripts (x" ~ 10^3 * ")")
    )

transcript_classification_hist + count_vs_abundance + nexon_vs_abundance + plot_layout(widths = c(1.5, 1, 1))
ggsave("figures/supplementary/Other_transcript_fig_combined.pdf", width = 7.5, height = 2.5)

# FSM transcript gene type percentages

read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    filter(
        structural_category == "full-splice_match"
    ) %>%
    distinct(associated_transcript) %>%
    left_join(
        GENCODE_gtf[, c("transcript_id", "gene_type")],
        join_by(associated_transcript == transcript_id)
    ) %>%
    distinct(associated_transcript, .keep_all=TRUE) %>%
    mutate(gene_type = if_else(gene_type %in% c("processed_pseudogene", "unprocessed_pseudogene"), "pseudogene", gene_type)) %>%
    filter(
        gene_type %in% c("lncRNA", "protein_coding", "pseudogene")
    ) %>%
    mutate(gene_type = factor(gene_type, levels = c("protein_coding", "lncRNA", "pseudogene"))) %>%
    ggplot(aes(x = gene_type)) +
        geom_bar(fill = "#009E73") +
        labs(x="", y="# FSM transcripts") +
        coord_flip()
ggsave("figures/supplementary/FSM_transcripts_gene_type.pdf", width = 8, height = 7)

# FSM transcript subcategory percentages

GENCODE_gtf <- rtracklayer::import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf")) %>%
    as_tibble()

read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    filter(
        structural_category == "full-splice_match"
    ) %>%
    distinct(associated_transcript) %>%
    left_join(
        GENCODE_gtf[, c("transcript_id", "gene_type")],
        join_by(associated_transcript == transcript_id)
    ) %>%
    distinct(associated_transcript, .keep_all=TRUE) %>%
    mutate(gene_type = if_else(gene_type %in% c("processed_pseudogene", "unprocessed_pseudogene"), "pseudogene", gene_type)) %>%
    filter(
        gene_type %in% c("lncRNA", "protein_coding", "pseudogene")
    ) %>%
    mutate(gene_type = factor(gene_type, levels = c("protein_coding", "lncRNA", "pseudogene"))) %>%
    ggplot(aes(x = gene_type)) +
        geom_bar(fill = "#009E73") +
        labs(x="", y="# FSM transcripts") +
        coord_flip()
ggsave("figures/supplementary/FSM_transcripts_gene_type.pdf", width = 8, height = 7)

