library(dplyr)
library(ggplot2)
library(arrow)
library(scales)
library(patchwork)
library(RColorBrewer)

colorVector <- brewer.pal(5, "Set2")

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
    filter(structural_category %in% c("fusion", "genic", "antisense", "intergenic")) %>% 
    mutate(
        structural_category = factor(structural_category, levels = c("fusion", "genic", "antisense", "intergenic"))
    ) %>% 
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
        labs(x = NULL, y = expression("Transcripts (x" ~ 10^3 * ")")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Counts vs abundance

row_sum <- read_parquet("nextflow_results/V47/final_expression.parquet") %>%
    mutate(
        row_sum = rowSums(across(where(is.numeric)))
    ) %>%
    dplyr::select(row_sum)

count_vs_abundance <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    bind_cols(row_sum) %>%
    filter(
        structural_category %in% c("fusion", "genic", "antisense", "intergenic"),
        row_sum > 10,
        ) %>%
    mutate(
        structural_category = factor(structural_category, levels = c("fusion", "genic", "antisense", "intergenic"))
    ) %>%         
    ggplot(aes(x = row_sum, fill = structural_category)) +
    geom_histogram(position = position_fill(), alpha = .75) +
    scale_x_log10() +
    labs(x = "Total observed counts", y = "Transcript proportion")

# Number of exons vs abudance

nexon_vs_abundance <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    filter(structural_category %in% c("fusion", "genic", "antisense", "intergenic")) %>% 
    mutate(
        structural_category = factor(structural_category, levels = c("fusion", "genic", "antisense", "intergenic"))
    ) %>%     
    ggplot(aes(x = exons, fill = structural_category)) +
    geom_histogram(alpha = .75, binwidth = 1) +
    scale_x_continuous(limits = c(1, 40)) +
    scale_y_continuous(labels = function(x) x / 1000) +
    labs(
        x = "# Exons",
        y = expression("Transcripts (x" ~ 10^3 * ")")
    )

transcript_classification_hist + count_vs_abundance + nexon_vs_abundance + plot_layout(widths = c(1, 1, 1))

ggsave("figures/supplementary/Other_transcript_fig_combined.pdf", width = 8, height = 3)
