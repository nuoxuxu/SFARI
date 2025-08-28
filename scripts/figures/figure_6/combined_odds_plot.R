library(dplyr)
library(readr)
library(rtracklayer)
library(ggplot2)
library(tidyr)
library(patchwork)
library(stringr)

my_theme <- theme_bw() +
    theme(
        axis.text.x = element_text(size = 20, vjust = 0.5, angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 25, color = "black"),
        strip.text = element_text(size = 20, color = "black"),
        strip.placement = "outside",
        legend.title = element_blank(),
        legend.position = c(0.08, 0.90),
        legend.text = element_text(size = 20)
    )

theme_set(my_theme)

combined_odds_conserved_vars <- read_tsv("data/combined_odds_conserved_vars.tsv") %>% 
    mutate(
        xlab = case_when(
            str_detect(Category, "_") ~ str_replace(Category, "^[^_]*_", ""),
            str_detect(Category, " ") ~ str_replace(Category, "^[^ ]* ", ""),
            TRUE ~ Category
        ),
        color = case_when(
            str_detect(Category, "Known") ~ "Known",
            str_detect(Category, "Novel") ~ "Novel",
            TRUE ~ "Other"
        )
    ) %>% 
    pivot_longer(c(Affected_with, Unaffected_with), names_to="status", values_to="number") %>% 
    filter(variant_type=="denovo") %>% 
    mutate(
        xlab = factor(xlab, levels = c("Acceptor", "Donor", "nonsense", "combined"))
    ) %>% 
    mutate(
      xlab = recode(xlab, "nonsense" = "Nonsense", "combined" = "Combined"),
      status = recode(status, "Affected_with" = "ASD proband", "Unaffected_with" = "Sibling"),
      Pvalue = ifelse(Pvalue < 0.05, format(Pvalue, scientific = TRUE, digits = 2), round(Pvalue, 2))
    ) %>% 
    mutate(
        status = factor(status, levels = c("Sibling", "ASD proband"))
    )

known_vars_plot <- combined_odds_conserved_vars %>% 
    filter(color=="Known") %>% 
    ggplot(aes(x=xlab, y=number, fill=status)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = c("ASD proband" = "#fdc086", "Sibling" = "#beaed4")) +
    geom_text(aes(label = paste0("OR = ", round(OddsRatio, 2))), vjust = -2, size = 6) +
    geom_text(aes(label = paste0("p = ", Pvalue)), vjust = -0.5, size = 6) +
    labs(x = NULL, y = "Individuals with >= 1 de novo variants") +
    ylim(c(0, 2400))
    
novel_vars_plot <- combined_odds_conserved_vars %>% 
    filter(color=="Novel") %>% 
    ggplot(aes(x=xlab, y=number, fill=status)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(values = c("ASD proband" = "#fdc086", "Sibling" = "#beaed4")) +
    geom_text(aes(label = paste0("OR = ", round(OddsRatio, 2))), vjust = -2, size=6) +
    geom_text(aes(label = paste0("p = ", Pvalue)), vjust = -0.5, size=6) +
    labs(x = NULL, y = NULL) +
    theme(
        legend.position = "none"
    ) +
    ylim(c(0, 60))

known_vars_plot + novel_vars_plot
ggsave("figures/figure_6/combined_odds_plot.pdf", width = 12, height = 5.5)
