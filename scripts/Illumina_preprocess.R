library(tximport)
library(GenomicFeatures)
library(jsonlite)
library(purrr)

txdb <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "Ensembl/Human/Release_103/Raw/Homo_sapiens.GRCh38.103.gtf.gz"))
tx2gene <- select(txdb, keys = transcripts(txdb)$tx_name, columns = "GENEID", keytype = "TXNAME")
files <- Sys.glob("data/short_read/salmon_outputs_second_ncRNA/salmon_results/*/quant.sf")
names(files) <- fromJSON("data/short_read/sample_name_mapping.json")[str_split(files, "/")%>%
    map_chr(5) %>%
    str_remove("_quant")] %>% 
    unname() %>%
    unlist()

# Write transcript TPM

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "scaledTPM", txOut = TRUE)
count_matrix <- txi$abundance
rownames(count_matrix) <- sapply(rownames(count_matrix), function(x) {sub("\\..*$", "", x)}) %>% unname()
count_matrix[, c("iPSC_1", "iPSC_2", "iPSC_3", "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3", "CN_1_2", "CN_1_3", "CN_2_1", "CN_2_2", "CN_3_1", "CN_3_2")] %>% 
    write.csv("data/short_read/transcript_tpm.csv")

# Write transcript count matrix

count_matrix <- txi$counts %>%
    round() %>%
    as.data.frame() %>% 
    .[, c("iPSC_1", "iPSC_2", "iPSC_3", "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3", "CN_1_2", "CN_1_3", "CN_2_1", "CN_2_2", "CN_3_1", "CN_3_2")]
write.csv(count_matrix, "data/short_read/transcript_count_matrix.csv")

# Write gene TPM

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, countsFromAbundance = "scaledTPM")
count_matrix <- txi$abundance
rownames(count_matrix) <- sapply(rownames(count_matrix), function(x) {sub("\\..*$", "", x)}) %>% unname()
count_matrix %>% 
    .[, c("iPSC_1", "iPSC_2", "iPSC_3", "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3", "CN_1_2", "CN_1_3", "CN_2_1", "CN_2_2", "CN_3_1", "CN_3_2")] %>% 
    write.csv("data/short_read/gene_tpm.csv")

# Write gene count matrix

count_matrix <- txi$counts %>%
    round() %>%
    as.data.frame() %>% 
    .[, c("iPSC_1", "iPSC_2", "iPSC_3", "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3", "CN_1_2", "CN_1_3", "CN_2_1", "CN_2_2", "CN_3_1", "CN_3_2")]
write.csv(count_matrix, "data/short_read/gene_count_matrix.csv")