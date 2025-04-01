library(arrow)
library(stringr)
library(DESeq2)
library(dplyr)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
    )

theme_set(my_theme)

full_expression <- read_parquet("nextflow_results/V47/full_expression.parquet")

full_expression <- as.data.frame(full_expression)

row.names(full_expression) <- full_expression$isoform

full_expression <- full_expression %>% select(-isoform)

cts <- as.matrix(full_expression)

coldata <- data.frame(
    row.names = colnames(cts),
    time_point = str_extract(colnames(full_expression), "^[^_]*"),
    replicate = str_extract(colnames(full_expression), "(?<=_)[^_]+") 
)

coldata <- coldata %>% 
    mutate(
        time_point = factor(time_point),
        replicate = factor(replicate)
    )

dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ time_point + replicate
)

vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup = c("time_point"))

ggsave("figures/supplementary/lr_pca.pdf",  width = 5, height = 4)