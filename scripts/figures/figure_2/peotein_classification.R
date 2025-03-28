library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(rtracklayer)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

protein_classification <- read_tsv("nextflow_results/V47/orfanage/SFARI.protein_classification.tsv") %>%
    # distinct(base_isoform, .keep_all = TRUE) %>%
    mutate(structural_category2 = if_else(
        protein_classification_base %in% c("pNIC", "pFSM", "pISM", "pNNC"),
        protein_classification_base,
        "Other"
    )) %>% 
    distinct(base_isoform, .keep_all = TRUE)

summary_df <- protein_classification %>%
    group_by(structural_category2) %>%
    summarise(len = n(), .groups = "drop") %>%
    mutate(percentage = len / sum(len) * 100)

colorVector <- c(
    "pFSM" = "#009E73",
    "pISM" = "#0072B2",
    "pNIC" = "#D55E00",
    "pNNC" = "#E69F00",
    "Other" = "#000000"
)

summary_df$structural_category2 <- factor(summary_df$structural_category2, levels = c("pFSM", "pISM", "pNIC", "pNNC", "Other"))

summary_df %>%
    filter(structural_category2 %in% c("pFSM", "pNIC", "pNNC", "Other")) %>%
    ggplot(aes(x = structural_category2, y = len, fill = structural_category2)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = 2, colour = "white", size = 4) +
    scale_fill_manual("Structural Category", values = colorVector) +
    scale_y_continuous(labels = function(x) x / 1000) +
    labs(
        x = "",
        y = expression("Proteins (x" ~ 10^3 * ")")
        )

ggsave("figures/figure_2/protein_class.pdf", width = 3, height = 2.5)