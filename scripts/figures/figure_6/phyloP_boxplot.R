library(dplyr)
library(readr)
library(rtracklayer)
library(ggplot2)
library(tidyr)
library(patchwork)
library(stringr)

novel_maps_output <- read_tsv("/scratch/nxu/100KGP_splicing/nextflow_results/novel/maps_output.tsv") %>%
    mutate(
        type = "Novel"
    ) %>%
    filter(
        region %in% c("Acceptor", "Donor"),
        site >= -10
    )

known_maps_output <- read_tsv("/scratch/nxu/100KGP_splicing/nextflow_results/known/maps_output.tsv") %>%
    mutate(
        type = "Known"
    ) %>%
    filter(
        (site >= -10) | (region == "Coding")
    )

novel_coding_maps <- read_tsv("MAPS_output.csv") %>%
    filter(
        region_site == "CDS_novel"
    ) %>%
    mutate(
        region_site = c("CDS_novel" = "Coding")[region_site],
        type = "Novel",
        csq = c("synonymous_variant" = "Synonymous", "missense_variant" = "Missense", "stop_gained" = "Nonsense")[csq]
    ) %>%
    dplyr::rename(
        xlab = csq,
        region = region_site
    )

maps_output <- bind_rows(novel_maps_output, known_maps_output, novel_coding_maps) %>%
  mutate(
    zscore  = maps / se,
    p_value = 2 * pnorm(abs(zscore), lower.tail = FALSE)
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>% 
  mutate(
    p_value = ifelse(p_adj < 0.05, TRUE, FALSE)
  ) %>% 
  mutate(site = as.character(site)) %>% 
  mutate(site = coalesce(site, xlab)) %>% 
    mutate(
        site = factor(site, levels = c("-10", "-9", "-8", "-7", "-6", "-5", "-4", "-3", "-2", "-1", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
        "Synonymous", "Missense", "Nonsense"))
    )

my_theme <- theme_bw() +
    theme(
        axis.text.x = element_text(size = 20, vjust = 0.5, angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        strip.text = element_text(size = 20, color = "black"),
        strip.placement = "outside",
        legend.title = element_blank(),
        legend.position = c(0.08, 0.90),
        legend.text = element_text(size = 20),
        # panel.grid.major.x = element_line(size = 1),
        # panel.spacing.y = unit(1, "lines")
    )

theme_set(my_theme)

maps_plot <- maps_output %>%
    mutate(region = factor(region, levels = c("Acceptor", "Donor", "Coding"))) %>%
    ggplot(aes(site, maps, color = type)) +
    geom_pointrange(
        aes(
            ymin = maps - 1.96 * se,
            ymax = maps + 1.96 * se,
            shape = p_value
        ),
        size = 0.8
    ) +
    geom_hline(yintercept = 0.051, linetype = 3, color = "#00BFC4", size = 1.2) +
    ylab("Mutational Constraint (MAPS)") +
    scale_color_manual(values = c("Novel" = "#F8766D", "Known" = "#00BFC4")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1)) +
    scale_x_discrete(breaks = levels(maps_output$site)[seq(1, length(levels(maps_output$site)), by = 5)]) +
    facet_grid(cols=vars(region), rows=vars(type), scales = "free", space="free_x") +
    theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
    )

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

phyloP_boxplot <- df %>%
    mutate(column = factor(column, levels = c("donor", "acceptor", "CDS", "UTR"))) %>%
    ggplot(aes(column, phyloP, fill = type)) +
    geom_boxplot(alpha = 0.5, outliers=FALSE) +
    geom_hline(yintercept = 2.27,
               colour = "red",
               linetype = "dashed",
               linewidth = 1) +
    labs(x = NULL, y="Evolutionary Conservation (phyloP)") +
    ylim(c(-8, 10))

maps_plot + phyloP_boxplot + plot_layout(widths = c(1.5, 1))
ggsave("figures/figure_6/phyloP_boxplot.pdf", width = 14, height = 8)