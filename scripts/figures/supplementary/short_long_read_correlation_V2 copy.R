library(ggplot2)
library(tximport)
library(stringr)
library(purrr)
library(GenomicFeatures)
library(arrow)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(reticulate)

# Input file paths

path_to_classification <- "nextflow_results/V47/final_classification.parquet"
path_to_gtf <- paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf")
path_to_salmon <- "nextflow_results/salmon_GENCODE47/*/quant.sf"
path_to_lr_exp_parquet <- "nextflow_results/V47/final_expression.parquet"

# Get classification

classification <- read_parquet(path_to_classification)

# Get short-read count matrix

txdb <- makeTxDbFromGFF(path_to_gtf)

tx2gene <- AnnotationDbi::select(txdb, keys = transcripts(txdb)$tx_name, columns = "GENEID", keytype = "TXNAME")

files <- Sys.glob(path_to_salmon)

names(files) <- str_split(files, "/") %>%
    map_chr(3)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = TRUE, txIdCol = "isoform")

sr_expression <- txi$counts %>% as_tibble(rownames = "isoform")

sr_expression$isoform <- str_extract(sr_expression$isoform, "^[^|]+")

# Get long-read count matrix

lr_expression <- read_parquet(path_to_lr_exp_parquet) %>% 
    rename(
        NPC_1_2 = NPC_1_3,
        NPC_3_2 = NPC_3_3,
        CN_1_1 = CN_1_2,
        CN_1_2 = CN_1_3
    ) %>% 
    left_join(select(classification, c(associated_transcript, isoform)), join_by(isoform == isoform)) %>% 
    filter(associated_transcript!="novel") %>% 
    select(-isoform) %>% 
    rename(isoform = associated_transcript) %>% 
    select(isoform, everything())

# Get correlation matrix

cor_matrix <- lr_expression %>%
    left_join(sr_expression, join_by(isoform == isoform), suffix = c("", "_illumina")) %>%
    left_join(select(classification, c(structural_category, isoform)), join_by(isoform == isoform)) %>% 
    filter(structural_category=="full-splice_match") %>% 
    dplyr::select(where(is.numeric)) %>%
    cor()

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

# What are these transcripts that are detected only in long-read but not in short-read?

library(rtracklayer)

gencode_gtf <- rtracklayer::import(path_to_gtf) %>% as_tibble() %>% filter(type == "transcript")

tx <- setdiff(sr_expression$isoform, lr_expression$isoform)

setdiff(sr_expression$isoform, lr_expression$isoform) %>% length()
setdiff(lr_expression$isoform, sr_expression$isoform) %>% length()
intersect(lr_expression$isoform, sr_expression$isoform) %>% length()
gencode_gtf %>% 
    filter(
        transcript_id %in% tx
    )
