library(dplyr)
library(arrow)
library(ggplot2)
library(readr)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none"
    )

theme_set(my_theme)

classification <- read_parquet("nextflow_results/V47/final_classification.parquet")

classification %>%
    group_by(associated_gene) %>%
    summarise(len = n(), .groups = "drop") %>%
    group_by(len) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
        len = case_when(
            len == 1 ~ "1",
            len %in% c(2, 3) ~ "2-3",
            len %in% c(4, 5) ~ "4-5",
            len %in% c(6, 7) ~ "6-7",
            len >= 8 ~ ">=8"
        )
    ) %>%
    group_by(len) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(
        len = factor(len, levels = c("1", "2-3", "4-5", "6-7", ">=8"))
    ) %>%
    mutate(
        percentage = (count / sum(count)) * 100
    ) %>%
    ggplot(aes(x = len, y = count)) +
        geom_col() +
        geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = 2, colour = "white", size = 3) +
        scale_y_continuous(labels = function(x) x / 1000) +
        labs(
            x = "Isoforms per gene",
            y = expression("Transcripts (x" ~ 10^3 * ")")
        )
ggsave("figures/figure_1/isoforms_per_gene.pdf", width = 3, height = 2.5)

# classification <- read_tsv("nextflow_results/V47/merged_collapsed_classification.filtered_lite_classification.txt")

# classification %>%
#     filter(structural_category %in% c(
#         "full-splice_match",
#         "incomplete-splice_match",
#         "novel_in_catalog",
#         "novel_not_in_catalog"
#     )) %>%
#     group_by(associated_gene) %>%
#     summarise(len = n(), .groups = "drop") %>%
#     group_by(len) %>%
#     summarise(count = n(), .groups = "drop") %>%
#     mutate(
#         len = case_when(
#             len == 1 ~ "1",
#             len %in% c(2, 3) ~ "2-3",
#             len %in% c(4, 5) ~ "4-5",
#             len %in% c(6, 7) ~ "6-7",
#             len >= 8 ~ ">=8"
#         )
#     ) %>%
#     group_by(len) %>%
#     summarise(count = sum(count), .groups = "drop") %>%
#     mutate(
#         len = factor(len, levels = c("1", "2-3", "4-5", "6-7", ">=8"))
#     ) %>%
#     ggplot(aes(x = len, y = count)) +
#     geom_col()

# classification %>%
#     filter(structural_category %in% c(
#         "full-splice_match",
#         "incomplete-splice_match",
#         "novel_in_catalog",
#         "novel_not_in_catalog"
#     )) %>%
#     group_by(associated_gene) %>%
#     summarise(len = n(), .groups = "drop") %>%
#     group_by(len) %>%
#     summarise(count = n(), .groups = "drop") %>%
#     ggplot(aes(len, count)) +
#         geom_col()