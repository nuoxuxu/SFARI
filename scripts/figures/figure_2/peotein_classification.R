library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(rtracklayer)
library(reticulate)

use_condaenv("/scratch/s/shreejoy/nxu/SFARI/envs/r_env")

py_run_string("
from src.utils import collapse_isoforms_to_proteoforms, read_gtf
import polars as pl

tx_classification = pl.read_parquet('nextflow_results/V47/final_classification.parquet')

orfanage_gtf = read_gtf('nextflow_results/V47/orfanage/orfanage.gtf', attributes=['gene_id', 'transcript_id'])
protein_classification = pl.read_csv('nextflow_results/V47/orfanage/SFARI.protein_classification.tsv', separator='\t')
protein_classification = (
    protein_classification
    .join(
        collapse_isoforms_to_proteoforms(orfanage_gtf).rename({'isoform': 'pb'}),
        on = 'pb',
        how = 'left'        
    )
)
")

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

protein_classification <- as_tibble(py$protein_classification$to_pandas()) %>%
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
    geom_text(aes(label = paste0(round(percentage, 1), "%")), position = position_dodge(0.9), vjust = 0, colour = "black", size = 4) +
    scale_fill_manual("Structural Category", values = colorVector) +
    scale_y_continuous(labels = function(x) x / 1000) +
    labs(
        x = "",
        y = expression("Proteins (x" ~ 10^3 * ")")
        )

ggsave("figures/figure_2/protein_class.pdf", width = 3, height = 2.5)