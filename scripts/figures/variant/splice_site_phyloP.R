library(tidyverse)

my_theme <- theme_bw() +
    theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.text.x = element_text(size = 14, vjust = 0.5, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 18),
        strip.placement = "outside",
        legend.title = element_blank(),
        legend.position = c(0.08, 0.90),
        legend.background = element_blank(),
        legend.text = element_text(size = 11),
        panel.spacing.y = unit(1, "lines")
    )

theme_set(my_theme)

read_csv("export/variant/CDS_ss_phyloP.csv") %>%
    ggplot(aes(pos, phyloP)) +
    geom_pointrange(
        aes(ymin = phyloP - 1.96 * sem,
            ymax = phyloP + 1.96 * sem
        ),
        size=.1,
        position = position_dodge2(0.5)
    ) +
    scale_x_continuous(breaks = seq.int(-25, 10, 5)) +
    facet_grid(
            cols = vars(factor(region, levels = c("Acceptor", "Donor"))),
            rows = vars(factor(spl_type, levels = c("known", "novel_5prime", "novel_3prime", "novel_both"))),
            scale = "free", space = "free_x"
        ) +
    ggtitle("PhyloP scores for splice sites and near-splice regions in CDS") +
    labs(x=NULL, y=NULL)
ggsave(filename = "figures/variant/CDS_ss_phyloP.pdf", width = 12, height = 9, dpi = 300)

read_csv("export/variant/exon_ss_phyloP.csv") %>%
    ggplot(aes(pos, phyloP)) +
    geom_pointrange(
        aes(ymin = phyloP - 1.96 * sem,
            ymax = phyloP + 1.96 * sem
        ),
        size=.1,
        position = position_dodge2(0.5)
    ) +
    scale_x_continuous(breaks = seq.int(-25, 10, 5)) +
    facet_grid(
            cols = vars(factor(region, levels = c("Acceptor", "Donor"))),
            rows = vars(factor(spl_type, levels = c("known", "novel_5prime", "novel_3prime", "novel_both"))),
            scale = "free", space = "free_x"
        ) +
    ggtitle("PhyloP scores for splice sites and near-splice regions") +
    labs(x=NULL, y=NULL)        
ggsave(filename = "figures/variant/exon_ss_phyloP.pdf", width = 12, height = 9, dpi = 300)

read_csv("export/variant/exon_ss_phyloP_canonical.csv") %>%
    ggplot(aes(pos, phyloP)) +
    geom_pointrange(
        aes(ymin = phyloP - 1.96 * sem,
            ymax = phyloP + 1.96 * sem
        ),
        size=.1,
        position = position_dodge2(0.5)
    ) +
    scale_x_continuous(breaks = seq.int(-25, 10, 5)) +
    facet_grid(
            cols = vars(factor(region, levels = c("Acceptor", "Donor"))),
            rows = vars(factor(spl_type, levels = c("known", "novel_5prime", "novel_3prime", "novel_both"))),
            scale = "free", space = "free_x"
        ) +
    ggtitle("PhyloP scores for all canonical splice sites and their near-splice regions") +
    labs(x=NULL, y=NULL)
ggsave(filename = "figures/variant/exon_ss_phyloP_canonical.pdf", width = 12, height = 9, dpi = 300)

read_csv("export/variant/CDS_ss_phyloP_canonical.csv") %>%
    ggplot(aes(pos, phyloP)) +
    geom_pointrange(
        aes(ymin = phyloP - 1.96 * sem,
            ymax = phyloP + 1.96 * sem
        ),
        size=.1,
        position = position_dodge2(0.5)
    ) +
    scale_x_continuous(breaks = seq.int(-25, 10, 5)) +
    facet_grid(
            cols = vars(factor(region, levels = c("Acceptor", "Donor"))),
            rows = vars(factor(spl_type, levels = c("known", "novel_5prime", "novel_3prime", "novel_both"))),
            scale = "free", space = "free_x"
        ) +
    ggtitle("PhyloP scores for CDS canonical splice sites and their near-splice regions") +
    labs(x=NULL, y=NULL)
ggsave(filename = "figures/variant/CDS_ss_phyloP_canonical.pdf", width = 12, height = 9, dpi = 300)