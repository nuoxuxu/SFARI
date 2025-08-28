library(ggplot2)
library(tximport)
library(stringr)
library(purrr)
library(GenomicFeatures)
library(arrow)
library(RColorBrewer)
library(pheatmap)
library(reticulate)
library(dplyr)

classification <- read_parquet("nextflow_results/V47/final_classification.parquet")

# Get short-read count matrix

txdb <- makeTxDbFromGFF("nextflow_results/V47/final_transcripts.gtf")

tx2gene <- AnnotationDbi::select(txdb, keys = transcripts(txdb)$tx_name, columns = "GENEID", keytype = "TXNAME")

files <- Sys.glob("nextflow_results/salmon_results/*/quant.sf")

names(files) <- str_split(files, "/") %>%
    map_chr(3)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = TRUE, txIdCol = "isoform")

sr_expression_tpm <- txi$abundance %>% as_tibble(rownames = "isoform")

# Get long-read count matrix

lr_expression <- read_parquet("nextflow_results/V47/final_expression.parquet") %>% 
    dplyr::rename(
        NPC_1_2 = NPC_1_3,
        NPC_3_2 = NPC_3_3,
        CN_1_1 = CN_1_2,
        CN_1_2 = CN_1_3
    )

lr_expression_cpm <- lr_expression %>%
  mutate(across(
    where(is.numeric),
    ~ .x / sum(.x) * 1e6
  ))

# Get correlation matrix

cor_matrix <- lr_expression_cpm %>%
    left_join(sr_expression_tpm, join_by(isoform == isoform), suffix = c("", "_illumina")) %>%
    left_join(select(classification, c(structural_category, isoform)), join_by(isoform == isoform)) %>% 
    filter(structural_category=="full-splice_match") %>% 
    dplyr::select(where(is.numeric)) %>%
    cor(method="pearson", use="complete.obs")

sample_names <- c("iPSC_1", "iPSC_2", "iPSC_3", "CN_1_1", "CN_1_2", "CN_2_1", "CN_2_2", "CN_3_1", "CN_3_2", "NPC_1_1", "NPC_1_2", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_2")
sample_names <- c(sample_names, str_c(sample_names, "_illumina"))

cor_matrix <- cor_matrix[sample_names, sample_names]

# Plotting using pheatmap

annotation <- data.frame(
    row.names = rownames(cor_matrix),
    time_point = str_extract(str_remove(rownames(cor_matrix), "_illumina"), "^[^_]*"),
    technology = c(rep("pacbio", 15), rep("illumina", 15))
)

pheatmap(
    cor_matrix,
    color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(10),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = annotation,
    annotation_row = annotation,
    show_rownames = FALSE,
    show_colnames = FALSE,
    display_numbers = TRUE,
    number_color = "black",
    fontsize_number = 6,
    fontsize_row = 10,
    fontsize_col = 10,
    width = 10,
    height = 10,
    filename = "figures/supplementary/correlation_matrix_transcript.png"
    )