library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(stringr)
library(glue)

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

burdens_by_variant <- read_tsv("export/variant/burdens_by_variant_type_ALL_COHORTS_MASTER.tsv") %>% 
    pivot_longer(c(Affected_with, Unaffected_with), names_to="status", values_to="number") %>% 
    mutate(
        Category = recode(Category, "donor"="Donor splice sites", "acceptor"="Acceptor splice sites"),
        status = recode(status, "Affected_with" = "Affected", "Unaffected_with" = "Unaffected"),
        Pvalue = ifelse(Pvalue < 0.05, format(Pvalue, scientific = TRUE, digits = 2), round(Pvalue, 2)),
    ) %>% 
    mutate(
        Category = factor(Category, levels = c("Acceptor splice sites", "Donor splice sites", "Nonsense", "Combined"))
    ) %>% 
    mutate(
        status = factor(status, levels = c("Unaffected", "Affected"))
    )    

plot_cohort <- function(cohort_name) {
    cohort_known <- burdens_by_variant %>% 
        filter(Cohort == cohort_name) %>% 
        filter(IsKnown == "Known") 
    ymax <- max(cohort_known$number) * 1.3
    p1 <- cohort_known %>% 
        ggplot(aes(x=Category, y=number, fill=status)) +
        geom_col(position = position_dodge()) +
        scale_fill_manual(values = c("Affected" = "#fdc086", "Unaffected" = "#beaed4")) +
        geom_text(aes(label = paste0("OR = ", round(OddsRatio, 2))), vjust = -2, size = 6) +
        geom_text(aes(label = paste0("p = ", Pvalue)), vjust = -0.5, size = 6) +
        labs(x = NULL, y = glue("Individuals with >= 1 \nde novo variants\nin {cohort_name}")) +
        ylim(0, ymax)
    
    cohort_novel <- burdens_by_variant %>% 
        filter(Cohort == cohort_name) %>% 
        filter(IsKnown == "Novel")
    ymax <- max(cohort_novel$number) * 1.3
    p2 <- cohort_novel %>% 
        ggplot(aes(x=Category, y=number, fill=status)) +
        geom_col(position = position_dodge()) +
        scale_fill_manual(values = c("Affected" = "#fdc086", "Unaffected" = "#beaed4")) +
        geom_text(aes(label = paste0("OR = ", round(OddsRatio, 2))), vjust = -2, size = 6) +
        geom_text(aes(label = paste0("p = ", Pvalue)), vjust = -0.5, size = 6) +
        labs(x = NULL, y = glue("Individuals with >= 1 \nde novo variants\nin {cohort_name}")) +
        ylim(0, ymax)
    p1 + p2
}

plot_cohort("MSSNG") / plot_cohort("SSC") / plot_cohort("SPARK_WGS") / plot_cohort("ALL_WES") + 
    plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("figures/supplementary/burden_test_by_cohort.pdf", width = 16, height = 22)
