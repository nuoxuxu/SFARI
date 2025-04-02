library(ggalluvial)
library(readr)
library(dplyr)

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

protein_class <- read_tsv("nextflow_results/V47/orfanage/SFARI.protein_classification.tsv")

protein_class %>%
    mutate(
        structural_category2 = if_else(
            tx_cat %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            tx_cat,
            "Other"
        ),
        protein_classification_base = if_else(
            protein_classification_base %in% c("pNIC", "pFSM", "pISM", "pNNC"),
            protein_classification_base,
            "Other"
        ),
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%  
    group_by(structural_category2, protein_classification_base) %>% 
    summarise(Freq = n()) %>%
    mutate(
        protein_classification_base = factor(protein_classification_base, levels = c("pFSM", "pISM", "pNIC", "pNNC", "Other")),
        structural_category2 = factor(structural_category2, levels = c("FSM", "ISM", "NIC", "NNC", "Other"))
    ) %>%
    ggplot(
        aes(axis1 = structural_category2, axis2 = protein_classification_base, y = Freq)
        ) +
    geom_alluvium(aes(fill=structural_category2)) +
    geom_stratum(aes(fill=structural_category2)) +
    scale_x_discrete(limits = c("Transcript class", "Protein Class"), expand = c(.2, .05)) +
    scale_y_continuous(labels = function(x) x / 1000) +
    scale_fill_manual("Structural\nCategory", values = colorVector) + 
    labs(y=expression("Transcripts (x" ~ 10^3 * ")")) +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20)
    ) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    ggtitle("Transcript and Protein Classification")
ggsave("figures/supplementary/sankey.pdf", width = 6, height = 5)