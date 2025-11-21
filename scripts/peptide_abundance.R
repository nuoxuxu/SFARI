library(arrow)
library(dplyr)
library(ggplot2)
library(ggpubr)

peptide_mapping_for_plotting <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping_for_plotting.parquet")

peptide_mapping_for_plotting %>% 
    ggplot(
        aes(x=GENCODE, y=mean_expression)
    ) +
    geom_boxplot() +
    theme_minimal()

my_comparisons <- list( c("GENCODE", "Novel") )

peptide_mapping_for_plotting %>% 
    mutate(
        type = ifelse(GENCODE, "GENCODE", "Novel")
    ) %>% 
    ggboxplot(
        x = "type", y = "mean_expression",
        color = "type", fill = "type", palette =c("#00AFBB", "#E7B800"),
        outliers = FALSE, method = "t.test"
    ) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")
ggsave("figures/test.png")
