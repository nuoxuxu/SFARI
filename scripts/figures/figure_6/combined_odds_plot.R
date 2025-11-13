library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(stringr)

my_theme <- theme_bw() +
    theme(
        axis.text.x = element_text(size = 20, vjust = 0.5, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        strip.text = element_text(size = 20, color = "black"),
        strip.placement = "outside",
        legend.title = element_blank(),
        legend.position = c(0.08, 0.90),
        legend.text = element_text(size = 20)
    )

theme_set(my_theme)

combined_odds_conserved_vars <- read_tsv("data/burdens_by_variant_type_known_vs_novel_WITHOUT_missense.tsv") %>% 
    pivot_longer(c(Affected_with, Unaffected_with), names_to="status", values_to="number") %>% 
    mutate(
      Category = recode(Category, "nonsense" = "Nonsense", "combined" = "Combined", "donor" = "Donor splice sites", "acceptor" = "Acceptor splice sites"),
      status = recode(status, "Affected_with" = "ASD proband", "Unaffected_with" = "Sibling"),
      Pvalue = ifelse(Pvalue < 0.05, format(Pvalue, scientific = TRUE, digits = 2), round(Pvalue, 2))
    ) %>% 
    mutate(
        Category = factor(Category, levels = c("Acceptor splice sites", "Donor splice sites", "Nonsense", "Combined"))
    ) %>% 
    mutate(
        status = factor(status, levels = c("Sibling", "ASD proband"))
    )

known_vars_plot <- combined_odds_conserved_vars %>% 
    filter(IsKnown=="Known") %>% 
    ggplot(aes(x=Category, y=number, fill=status)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = c("ASD proband" = "#fdc086", "Sibling" = "#beaed4")) +
    geom_text(aes(label = paste0("OR = ", round(OddsRatio, 2))), vjust = -2, size = 6) +
    geom_text(aes(label = paste0("p = ", Pvalue)), vjust = -0.5, size = 6) +
    labs(x = NULL, y = "Individuals with >= 1 de novo variants") +
    ylim(0, 3150)
    
novel_vars_plot <- combined_odds_conserved_vars %>% 
    filter(IsKnown=="Novel") %>% 
    ggplot(aes(x=Category, y=number, fill=status)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = c("ASD proband" = "#fdc086", "Sibling" = "#beaed4")) +
    geom_text(aes(label = paste0("OR = ", round(OddsRatio, 2))), vjust = -2, size=6) +
    geom_text(aes(label = paste0("p = ", Pvalue)), vjust = -0.5, size=6) +
    labs(x = NULL, y = NULL) +
    theme(
        legend.position = "none"
    ) +
    ylim(0, 210)

known_vars_plot + novel_vars_plot
ggsave("figures/figure_6/combined_odds_plot.pdf", width = 16, height = 5.5)
