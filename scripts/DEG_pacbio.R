library(dplyr)
library(biomaRt)
library(stringr)
library(edgeR)
library(arrow)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)

gene_counts <- read_parquet("proc/pacbio_count_matrix.parquet")

gene_info <- gene_counts %>% 
  dplyr::select(Gene_ID) %>% 
  left_join(ensembl, by = c("Gene_ID" = "ensembl_gene_id"))

gene_counts <- gene_counts[-1]

metadata <- data.frame(
    sample_name = colnames(gene_counts),
    CT = str_extract(colnames(gene_counts), "^[^_]*")
) %>%
    mutate(CT = ifelse(CT == "iNPC", "NPC", ifelse(CT == "CN", "neuron", CT)))

##### Set up Limma-voom model ####

cell_type <- factor(unlist(metadata["CT"]), levels = c("iPSC", "NPC", "neuron"))
log10_library_size <- log10(unlist((gene_counts %>% colSums)))
dge <- DGEList(gene_counts, genes = gene_info)

min_samples_expressing_gene <- metadata %>% nrow * 0.8
# The 0.8 here refers to the fraction of total samples that needs to express the gene

dge <- dge[rowSums(dge$counts >= 1) >= min_samples_expressing_gene, ]
# This step filters genes such that they need to be detected in at least 80% of samples

dge <- calcNormFactors(dge, method = "TMM")

##### Create model design, run voom, fit design and contrasts, run emperical Bayes ####

design <- model.matrix(~ 0 + cell_type + log10_library_size)
colnames(design) <- c("iPSC", "NPC", "Neuron", "log10_library_size")

# Perform voom transformation
vm <- voom(dge, design, plot = TRUE)

# Perform lmFit
fit <- lmFit(vm, design)

# Set up and fit contrasts
cont.matrix <- makeContrasts(NPC_vs_iPSC = NPC-iPSC, Neuron_vs_NPC = Neuron-NPC, Neuron_vs_iPSC = Neuron-iPSC,
                             levels=design)

fit <- contrasts.fit(fit, cont.matrix)

# Perform eBayes
fit <- eBayes(fit)

#### Run contrasts ####

npc_ipsc_table <- topTable(fit, coef = "NPC_vs_iPSC",  n = Inf, sort = "none", 
                         adjust.method = "BH")

neuron_npc_table <- topTable(fit, coef = "Neuron_vs_NPC",  n = Inf, sort = "none", 
                         adjust.method = "BH")

neuron_ipsc_table <- topTable(fit, coef = "Neuron_vs_iPSC",  n = Inf, sort = "none", 
                         adjust.method = "BH")

#### Combine contrast results ####

table_list <- list(npc_ipsc_table, neuron_ipsc_table, neuron_npc_table)

all_res <- table_list %>%
  purrr::reduce(left_join, by = "Gene_ID") %>%
  dplyr::select(Gene_ID, external_gene_name,
         NPC_iPSC.logFC = logFC.x, Neuron_iPSC.logFC = logFC.y, Neuron_NPC.logFC = logFC,
         NPC_iPSC.P.Value = P.Value.x, Neuron_iPSC.P.Value = P.Value.y, Neuron_NPC.P.Value = P.Value,
         NPC_iPSC.adj.P.Val = adj.P.Val.x, Neuron_iPSC.adj.P.Val = adj.P.Val.y, Neuron_NPC.adj.P.Val = adj.P.Val)

write.csv(all_res, "results/long_read/combined_batches_DE_results.csv")
