library(arrow)
library(readr)
library(dplyr)

colorVector <- c(
    "FSM" = "#009E73",
    "ISM" = "#0072B2",
    "NIC" = "#D55E00",
    "NNC" = "#E69F00"
)

structural_category_labels <- c(
    "full-splice_match"        = "FSM",
    "incomplete-splice_match"  = "ISM",
    "novel_in_catalog"         = "NIC",
    "novel_not_in_catalog"     = "NNC"
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

get_one_df <- function(structural_category, sum_to_100 = TRUE) {
    df <- tibble(
        x = seq(0, 1, 0.01),
        y = sapply(x, function(x) get_n_isoforms(filter(full_expression, structural_category == {{structural_category}}), x)),
        category = structural_category
    )
    if (!sum_to_100) {
        return(df)
    } else {
        total_n_tx <- full_expression %>%
            filter(structural_category == {{structural_category}}) %>%
            filter(rowSums(across(everything(), ~ . > 5)) > 2) %>%
            nrow()
        df %>%
            mutate(
                y = y / filter(df, x == 1)$y * 100
            )
    }

}

plot_saturation <- function(sum_to_100 = TRUE) {
    if (!sum_to_100) {
        df <- c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog", "incomplete-splice_match") %>% 
            lapply(function(x) get_one_df(x, sum_to_100 = FALSE)) %>% 
            bind_rows()
        df %>%
            mutate(
                category = structural_category_labels[category]
            ) %>% 
            ggplot(aes(x, y, color = category)) +
            geom_line() +
            scale_color_manual("Structural\nCategory", values = colorVector, ) +
            scale_y_continuous(labels = function(x) x / 1000) +
            labs(
                x = "Percentage of subsampled reads",
                y = expression("Transcripts (x" ~ 10^3 * ")")
            )
    } else {
        df <- c("full-splice_match", "novel_in_catalog", "novel_not_in_catalog", "incomplete-splice_match") %>% 
            lapply(function(x) get_one_df(x)) %>% 
            bind_rows()
        df %>%
            mutate(
                category = structural_category_labels[category]
            ) %>% 
            ggplot(aes(x, y, color = category)) +
            geom_line() +
            scale_color_manual("Structural\nCategory", values = colorVector, ) +
            labs(
                x = "Percentage of subsampled reads",
                y = "Percentage of transcripts"
            )
    }
}

plot_saturation(sum_to_100 = FALSE)
ggsave("figures/supplementary/long_read_saturation_y_percentage.pdf", width = 4, height = 3)

plot_saturation(sum_to_100 = TRUE)
ggsave("figures/supplementary/long_read_saturation_y_number.pdf", width = 4, height = 3)