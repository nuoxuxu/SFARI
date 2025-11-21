library(dplyr)
library(ggplot2)
library(arrow)
library(scales)
library(patchwork)

colorVector <- c(
    "FSM" = "#009E73",
    "ISM" = "#0072B2",
    "NIC" = "#D55E00",
    "NNC" = "#E69F00",
    "Other" = "#000000"
)

structural_category_labels <- c(
    "full-splice_match"        = "FSM",
    "incomplete-splice_match"  = "ISM",
    "novel_in_catalog"         = "NIC",
    "novel_not_in_catalog"     = "NNC",
    "Other"                    = "Other"
)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "none"
    )
update_geom_defaults("text", list(size = 5))    
theme_set(my_theme)

# Transcript Classification

transcript_classification_hist <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    mutate(
        # If category is in our known set, keep it; otherwise use "Other"
        structural_category2 = if_else(
            structural_category %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            structural_category,
            "Other"
        ),
        # Map to short labels (FSM, ISM, NIC, NNC, Other)
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%
    group_by(structural_category2) %>%
    summarise(len = n()) %>%
    mutate(
        percentage = (len / sum(len)) * 100,
        structural_category2 = factor(structural_category2, levels = c("FSM", "ISM", "NIC", "NNC", "Other"))
    ) %>%
    ggplot(aes(x = structural_category2, y = len, fill = structural_category2)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = 2, colour = "black", size = 3.5) +
        scale_y_continuous(labels = function(x) x / 1000) +
        scale_fill_manual("Structural Category", values = colorVector) +
        labs(x = "Structural Category", y = expression("Transcripts (x" ~ 10^3 * ")"))
# ggsave("figures/figure_1/transcript_classification_hist.pdf", width = 8, height = 7)

# Length vs abundance

read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    mutate(
        structural_category2 = if_else(
            structural_category %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            structural_category,
            "Other"
        ),
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%
    ggplot(aes(x = length, fill = structural_category2)) +
    geom_histogram(alpha = .75) +
    scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits = c(300, 10^4)) +
    scale_y_continuous(labels = function(x) x / 1000) +
    scale_fill_manual("Structural Category", values = colorVector) +
    labs(
        x = "Transcript Length (bp)",
        y = expression("Transcripts (x" ~ 10^3 * ")")
    ) +
    theme(
        legend.position = c(0.8, 0.8)
    )
ggsave("figures/figure_1/length_vs_abundance.pdf", width = 3.5, height = 3)

# Number of exons vs abudance

nexon_vs_abundance <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    mutate(
        structural_category2 = if_else(
            structural_category %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            structural_category,
            "Other"
        ),
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%
    ggplot(aes(x = exons, fill = structural_category2)) +
    geom_histogram(alpha = .75, binwidth = 1) +
    scale_x_continuous(limits = c(1, 40)) +
    scale_y_continuous(labels = function(x) x / 1000) +
    scale_fill_manual("Structural Category", values = colorVector) +
    labs(
        x = "# Exons",
        y = expression("Transcripts (x" ~ 10^3 * ")")
    )
# ggsave("figures/figure_1/nexon_vs_abundance.pdf", width = 8, height = 7)

# Counts vs abundance

row_sum <- read_parquet("nextflow_results/V47/final_expression.parquet") %>%
    mutate(
        row_sum = rowSums(across(where(is.numeric)))
    )
count_vs_abundance <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    mutate(
        structural_category2 = if_else(
            structural_category %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            structural_category,
            "Other"
        ),
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%
    left_join(row_sum[c("isoform", "row_sum")], by = "isoform") %>%
    filter(row_sum > 10) %>%
    ggplot(aes(x = row_sum, fill = structural_category2)) +
    geom_histogram(position = position_fill(), alpha = .75) +
    scale_x_log10() +
    scale_fill_manual("Structural Category", values = colorVector) +
    labs(x = "Total observed counts", y = "Transcript proportion")
# ggsave("figures/figure_1/count_vs_abundance.pdf", width = 8, height = 7)

transcript_classification_hist + count_vs_abundance + nexon_vs_abundance + plot_layout(widths = c(1, 1, 1))
ggsave("figures/figure_1/transcript_fig_combined_2.pdf", width = 7.5, height = 2.5)