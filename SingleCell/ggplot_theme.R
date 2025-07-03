library(ggplot2)
library(scales)

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
    "Other"                    = "Other",
    "fusion"                   = "Fusion",
    "genic"                    = "Genic",
    "intergenic"               = "Intergenic",
    "antisense"                = "Antisense",
    "moreJunctions"            = "moreJunctions"
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

library(arrow)
library(dplyr)

final_classification <- read_parquet("nextflow_results/V47/final_classification.parquet")

final_classification %>%
    mutate(
        structural_category = structural_category_labels[structural_category]
    )

final_classification %>% distinct(structural_category)
