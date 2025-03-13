library(arrow)
library(readr)
library(dplyr)

colorVector <- c(
    "FSM" = "#009E73",
    "ISM" = "#0072B2",
    "NIC" = "#D55E00",
    "NNC" = "#E69F00"
)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)
    )

theme_set(my_theme)

filtered_lite_classification <- read_tsv("nextflow_results/V47/merged_collapsed_classification.filtered_lite_classification.txt")

full_expression <- read_parquet("nextflow_results/V47/full_expression.parquet") %>% 
    filter(isoform %in% filtered_lite_classification$isoform) %>% 
    left_join(
        filtered_lite_classification[, c("isoform", "structural_category")],
        join_by(isoform == isoform)
    )

get_n_isoforms <- function(expression, percentage) {
  expression %>%
    dplyr::select(where(is.numeric)) %>%
    mutate(across(everything(), ~ . * percentage)) %>%
    filter(rowSums(across(everything(), ~ . > 5)) > 2) %>%
    nrow()
}

FSM <- tibble(
    x = seq(0, 1, 0.01),
    y = sapply(x, function(x) get_n_isoforms(filter(full_expression, structural_category=="full-splice_match"), x)),
    category = "FSM"
    )

ISM <- tibble(
    x = seq(0, 1, 0.01),
    y = sapply(x, function(x) get_n_isoforms(filter(full_expression, structural_category=="incomplete-splice_match"), x)),
    category = "ISM"
    )

NIC <- tibble(
    x = seq(0, 1, 0.01),
    y = sapply(x, function(x) get_n_isoforms(filter(full_expression, structural_category=="novel_in_catalog"), x)),
    category = "NIC"
    )

NNC <- tibble(
    x = seq(0, 1, 0.01),
    y = sapply(x, function(x) get_n_isoforms(filter(full_expression, structural_category=="novel_not_in_catalog"), x)),
    category = "NNC"
    )

df <- bind_rows(FSM, ISM, NIC, NNC)

df %>% 
    ggplot(aes(x, y, color=category)) +
    geom_line() +
    scale_color_manual("Structural\nCategory", values = colorVector) +
    labs(
        x = "Percentage of subsampled reads",
        y = expression("Transcripts (x" ~ 10^3 * ")")
        ) +
    scale_y_continuous(labels = function(x) x / 1000)

ggsave("figures/supplementary/long_read_saturation.pdf", width = 4, height = 3)
