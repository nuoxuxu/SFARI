library(ggplot2)
library(dplyr)
library(arrow)
library(edgeR)

expression <- read_parquet("nextflow_results/V47/final_expression.parquet")
DGEList(counts=as.matrix(expression[, grep("_", colnames(expression))]), genes=expression[, "isoform"]) %>%
    cpm(log=TRUE, prior.count=1) %>% 
    as.data.frame(row.names = expression$isoform) %>% 
    as_tibble(rownames="NA")

LR_SJ <- read_parquet("nextflow_results/V47/LR_SJ_2.parquet")

LR_SJ %>% 
    ggplot(aes(x=type, y=mean_log2_cpm, fill=SR)) +
    geom_boxplot() +
    labs(x = "Splice junctions", y = "Long-read RNA-seq expression\n(log2(CPM + 1))")
ggsave("figures/figure_1/sj_expression_by_whether_validated_2.png", width=4, height=4)

LR_patowary <- read_parquet("nextflow_results/V47/LR_patowary.parquet")

LR_patowary %>% 
    ggplot(aes(x=type, y=mean_log2_cpm, fill=supported)) +
    geom_boxplot() +
    labs(x = "long-read transcript types", y = "Long-read RNA-seq expression\n(log2(CPM + 1))")
ggsave("figures/figure_1/lr_expression_by_patowary_validation.png", width=4, height=4)
