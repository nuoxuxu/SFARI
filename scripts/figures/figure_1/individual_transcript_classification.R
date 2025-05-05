library(dplyr)
library(ggplot2)
library(arrow)
library(scales)
library(patchwork)
library(arrow)

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
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

# Transcript Classification

plot_transcript_classification_hist <- function(timepoint) {
    isoforms <- read_parquet("nextflow_results/V47/final_expression.parquet") %>% 
        select(isoform, starts_with({{timepoint}})) %>%
        filter(rowSums(select(., -isoform) >= 5) >= 2) %>%
        pull(isoform)

    read_parquet("nextflow_results/V47/final_classification.parquet") %>%
        filter(isoform %in% isoforms) %>%
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
            geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = -0.5, colour = "black", size = 3.5) +
            scale_y_continuous(labels = function(x) x / 1000) +
            scale_fill_manual("Structural Category", values = colorVector) +
            labs(x = "Structural Category", y = expression("Transcripts (x" ~ 10^3 * ")")) +
            ggtitle(timepoint)      
}

plot_transcript_classification_hist("iPSC") + plot_transcript_classification_hist("NPC") + plot_transcript_classification_hist("CN") + plot_layout(widths = c(1, 1, 1))
ggsave("figures/figure_1/transcript_classification_hist_by_timepoint.pdf", width = 7.5, height = 2.5)
