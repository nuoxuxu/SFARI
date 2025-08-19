library(dplyr)
library(readr)
library(rtracklayer)
library(ggplot2)
library(patchwork)

# Set ggplot theme

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)
    )

theme_set(my_theme)

# Essential splice sites

known_splice_sites_cds <- read_csv("nextflow_results/V47/orfanage/splice_sites/known_splice_sites_cds_phyloP.csv")
novel_splice_sites_cds <- read_csv("nextflow_results/V47/orfanage/splice_sites/novel_splice_sites_cds_phyloP.csv")

splice_sites_df <- bind_rows(known_splice_sites_cds, novel_splice_sites_cds) %>% 
    mutate(
        start_or_end = ifelse(start_or_end == "start", "donor", "acceptor"),
    ) %>%
    sample_frac(0.1)

splice_sites_df <- splice_sites_df %>%
    dplyr::rename(column = start_or_end) %>%
    dplyr::select(column, phyloP, type)

# Exonic regions

CDS_regions <- read_csv("nextflow_results/V47/orfanage/exonic_regions/CDS_regions_with_phyloP.csv") %>% 
    mutate(
        region = "CDS"
    )

UTR_regions <- read_csv("nextflow_results/V47/orfanage/exonic_regions/UTR_regions_with_phyloP.csv") %>% 
    mutate(
        region = "UTR"
    )

combined_exonic_regions <- bind_rows(CDS_regions, UTR_regions) %>% 
    dplyr::rename(column = region) %>%
    dplyr::select(column, phyloP, type)

df <- bind_rows(splice_sites_df, combined_exonic_regions)

df %>%
    mutate(column = factor(column, levels = c("donor", "acceptor", "CDS", "UTR"))) %>%
    ggplot(aes(column, phyloP, fill = type)) +
    geom_boxplot(alpha = 0.5, outliers=FALSE) +
    geom_hline(yintercept = 2.27,
               colour = "red",
               linetype = "dashed",
               linewidth = 1) +
    labs(x = NULL) +
    ylim(c(-8, 10))

df %>% 
    filter(column=="donor", type=="known") %>% 
    mutate(phyloP = ifelse(phyloP < 2.27, TRUE, FALSE)) %>%
    summarise(percentage = sum(phyloP) / n() * 100)

df %>% 
    filter(column=="acceptor", type=="known") %>% 
    mutate(phyloP = ifelse(phyloP < 2.27, TRUE, FALSE)) %>%
    summarise(percentage = sum(phyloP) / n() * 100)    

ggsave("figures/figure_6/phyloP_boxplot.pdf", width = 8, height = 5.5)