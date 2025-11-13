### CALCULATE CPMs, RUN DESeq2, AND MSsatsTMT AND CALCULATE CLUSTERS ####
# 
library(dplyr)
library(tidyr)
library(arrow)
library(edgeR)
library(DESeq2)
library(MSstatsTMT)
library(stringr)


setwd("")


######################

# CALCULATE TRANSCRIPT- AND GENE-LEVEL CPM

#####################

#Load in data
final_pb_ids <- read.table("./data/filtered_transcript_list.txt", 
                           header = T) #only using transcripts with an ORFanage-predicted ORF from FSM, ISM, NIC or NNC categories

#Subset the transcript file to these PB IDs. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)

filtered_tr_count <- tr_count[tr_count$isoform %in% final_pb_ids$x, ]

# CPM - TRANSCRIPTS -------------------------------------------------------
PB_counts <- as.matrix(filtered_tr_count[,c(2:16)])
rownames(PB_counts) <- filtered_tr_count$isoform


dge_PB <- DGEList(counts = PB_counts, genes = filtered_tr_count[,1])
dge_PB <- calcNormFactors(dge_PB)


PB_tr_cpm_norm <- cpm(dge_PB, normalized.lib.sizes=TRUE, log=FALSE)
PB_tr_cpm_norm <- as.data.frame(PB_tr_cpm_norm)
PB_tr_cpm_norm$pb_id <- row.names(PB_tr_cpm_norm)

#Export to save. 
write.csv(PB_tr_cpm_norm, "./code/expression/output/transcript_CPM.csv", row.names = F)


# CPM - GENES (SUMMED) ----------------------------------------------------
#Sum counts. 
classification <- read_parquet("./data/final_classification.parquet") #182371


#Subset the transcript file to these PB IDs. 
filtered_tr_count <- tr_count[tr_count$isoform %in% final_pb_ids$x, ]


filtered_tr_count <- merge(filtered_tr_count, classification[, c("isoform", "associated_gene")], 
                           by.x = "isoform", by.y = "isoform", all.x = TRUE)

gene_sum_count <- filtered_tr_count %>%
  group_by(associated_gene) %>%
  summarise(sum_iPSC_1 = sum(iPSC_1), 
            sum_iPSC_2 = sum(iPSC_2),
            sum_iPSC_3 = sum(iPSC_3), 
            sum_NPC_1_1 = sum(NPC_1_1), 
            sum_NPC_1_3 = sum(NPC_1_3), 
            sum_NPC_2_1 = sum(NPC_2_1), 
            sum_NPC_2_2 = sum(NPC_2_2), 
            sum_NPC_3_1 = sum(NPC_3_1), 
            sum_NPC_3_3 = sum(NPC_3_3), 
            sum_CN_1_2 = sum(CN_1_2), 
            sum_CN_1_3 = sum(CN_1_3), 
            sum_CN_2_1 = sum(CN_2_1), 
            sum_CN_2_2 = sum(CN_2_2), 
            sum_CN_3_1 = sum(CN_3_1), 
            sum_CN_3_2 = sum(CN_3_2)) %>%
  ungroup()



gene_sum_counts <- as.matrix(gene_sum_count[,c(2:16)])
rownames(gene_sum_counts) <- gene_sum_count$associated_gene

dge_gene_sum <- DGEList(counts = gene_sum_counts, genes = gene_sum_count[,1])
dge_gene_sum <- calcNormFactors(dge_gene_sum)

gene_sum_cpm_norm <- cpm(dge_gene_sum, normalized.lib.sizes=TRUE, log=FALSE)

gene_sum_cpm_norm <- as.data.frame(gene_sum_cpm_norm)
gene_sum_cpm_norm$gene_name <- row.names(gene_sum_cpm_norm)

#Export to save. 
write.csv(gene_sum_cpm_norm, "./code/expression/output/gene_CPM.csv", row.names = F)

# CPM - GENES (SUMMED - everything) ----------------------------------------------------
#Sum counts. 
classification <- read_parquet("./data/final_classification.parquet") #182371


#Subset the transcript file to these PB IDs.
tr_count <- tr_count[tr_count$isoform %in% classification$isoform, ]


tr_count <- merge(tr_count, classification[, c("isoform", "associated_gene")], 
                           by.x = "isoform", by.y = "isoform", all.x = TRUE)

gene_sum_count <- tr_count %>%
  group_by(associated_gene) %>%
  summarise(sum_iPSC_1 = sum(iPSC_1), 
            sum_iPSC_2 = sum(iPSC_2),
            sum_iPSC_3 = sum(iPSC_3), 
            sum_NPC_1_1 = sum(NPC_1_1), 
            sum_NPC_1_3 = sum(NPC_1_3), 
            sum_NPC_2_1 = sum(NPC_2_1), 
            sum_NPC_2_2 = sum(NPC_2_2), 
            sum_NPC_3_1 = sum(NPC_3_1), 
            sum_NPC_3_3 = sum(NPC_3_3), 
            sum_CN_1_2 = sum(CN_1_2), 
            sum_CN_1_3 = sum(CN_1_3), 
            sum_CN_2_1 = sum(CN_2_1), 
            sum_CN_2_2 = sum(CN_2_2), 
            sum_CN_3_1 = sum(CN_3_1), 
            sum_CN_3_2 = sum(CN_3_2)) %>%
  ungroup()



gene_sum_counts <- as.matrix(gene_sum_count[,c(2:16)])
rownames(gene_sum_counts) <- gene_sum_count$associated_gene

dge_gene_sum <- DGEList(counts = gene_sum_counts, genes = gene_sum_count[,1])
dge_gene_sum <- calcNormFactors(dge_gene_sum)

gene_sum_cpm_norm <- cpm(dge_gene_sum, normalized.lib.sizes=TRUE, log=FALSE)

gene_sum_cpm_norm <- as.data.frame(gene_sum_cpm_norm)
gene_sum_cpm_norm$gene_name <- row.names(gene_sum_cpm_norm)

#Export to save. 
write.csv(gene_sum_cpm_norm, "./code/expression/output/gene_CPM_all.csv", row.names = F)




######################

# RUN TRANSCRIPT- AND GENE-LEVEL DIFFERENTIAL EXPRESSION USING DESEQ2

#####################

#Load data. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)

final_pb_ids <- read.table("./data/filtered_transcript_list.txt", 
                           header = T) 

filtered_tr_count <- tr_count[tr_count$isoform %in% final_pb_ids$x, ]

# DESEQ2 TRANSCRIPTS ------------------------------------------------------
#Prep df.
rownames(filtered_tr_count) <- filtered_tr_count$isoform
filtered_tr_count <- filtered_tr_count[,-1]


#Get column data
colData <- data.frame(group = c(rep("t00", 3), rep("t04", 6), rep("t30", 6)))
rownames(colData) <- colnames(filtered_tr_count)
colData$group <- as.factor(colData$group)

#Construct DESEqDataSet object
transcript_dds <- DESeqDataSetFromMatrix(countData=filtered_tr_count, 
                                         colData=colData, 
                                         design=~ group)

transcript_dds <- DESeq(transcript_dds)

#Get each comparison.
res_t04_vs_t00 <- results(transcript_dds, contrast = c("group", "t04", "t00"))
summary(res_t04_vs_t00)

res_t30_vs_t00 <- results(transcript_dds, contrast = c("group", "t30", "t00"))
summary(res_t30_vs_t00)

res_t30_vs_t04 <- results(transcript_dds, contrast = c("group", "t30", "t04"))
summary(res_t30_vs_t04)

#Save as dfs.
res_t04_vs_t00_df <- as.data.frame(res_t04_vs_t00) 
res_t04_vs_t00_df$pb_id <- rownames(res_t04_vs_t00_df)
res_t30_vs_t00_df <- as.data.frame(res_t30_vs_t00)
res_t30_vs_t00_df$pb_id <- rownames(res_t30_vs_t00_df)
res_t30_vs_t04_df <- as.data.frame(res_t30_vs_t04)
res_t30_vs_t04_df$pb_id <- rownames(res_t30_vs_t04_df)


# Save the results as a CSV file
write.csv(res_t04_vs_t00_df, file = "./code/expression/output/DESeq2_tr_t04_vs_t00.csv",
          row.names = F)
write.csv(res_t30_vs_t00_df, file = "./code/expression/output/DESeq2_tr_t30_vs_t00.csv",
          row.names = F)
write.csv(res_t30_vs_t04_df, file = "./code/expression/output/DESeq2_tr_t30_vs_t04.csv",
          row.names = F)




# DESEQ2 GENE-LEVEL (SUMMED TR) ------------------------------------------------------
#Sum counts. 
classification <- read_parquet("./data/final_classification.parquet")

filtered_tr_count <- tr_count[tr_count$isoform %in% final_pb_ids$x, ]
filtered_tr_count <- merge(filtered_tr_count, classification[, c("isoform", "associated_gene")], 
                           by.x = "isoform", by.y = "isoform", all.x = TRUE)

gene_sum_count <- filtered_tr_count %>%
  group_by(associated_gene) %>%
  summarise(sum_iPSC_1 = sum(iPSC_1), 
            sum_iPSC_2 = sum(iPSC_2),
            sum_iPSC_3 = sum(iPSC_3), 
            sum_NPC_1_1 = sum(NPC_1_1), 
            sum_NPC_1_3 = sum(NPC_1_3), 
            sum_NPC_2_1 = sum(NPC_2_1), 
            sum_NPC_2_2 = sum(NPC_2_2), 
            sum_NPC_3_1 = sum(NPC_3_1), 
            sum_NPC_3_3 = sum(NPC_3_3), 
            sum_CN_1_2 = sum(CN_1_2), 
            sum_CN_1_3 = sum(CN_1_3), 
            sum_CN_2_1 = sum(CN_2_1), 
            sum_CN_2_2 = sum(CN_2_2), 
            sum_CN_3_1 = sum(CN_3_1), 
            sum_CN_3_2 = sum(CN_3_2)) %>%
  ungroup()


#Use summed counts for DESeq
gene_sum_count <- as.data.frame(gene_sum_count)
rownames(gene_sum_count) <- gene_sum_count$associated_gene
gene_sum_count <- gene_sum_count[,-1]

#Get col data
colData <- data.frame(group = c(rep("t00", 3), rep("t04", 6), rep("t30", 6)))
rownames(colData) <- colnames(gene_sum_count)
colData$group <- as.factor(colData$group)

#Construct DESEqDataSet object
dds_gene_sum <- DESeqDataSetFromMatrix(countData=gene_sum_count, 
                                       colData=colData, 
                                       design=~group)

dds_gene_sum <- DESeq(dds_gene_sum)

#Results: summed gene-level
gene_sum_t04_vs_t00 <- results(dds_gene_sum, contrast = c("group", "t04", "t00"))
summary(gene_sum_t04_vs_t00)

gene_sum_t30_vs_t00 <- results(dds_gene_sum, contrast = c("group", "t30", "t00"))
summary(gene_sum_t30_vs_t00)

gene_sum_t30_vs_t04 <- results(dds_gene_sum, contrast = c("group", "t30", "t04"))
summary(gene_sum_t30_vs_t04)

#Save as dfs.
gene_sum_t04_vs_t00_df <- as.data.frame(gene_sum_t04_vs_t00) 
gene_sum_t04_vs_t00_df$gene_name <- rownames(gene_sum_t04_vs_t00_df)

gene_sum_t30_vs_t00_df <- as.data.frame(gene_sum_t30_vs_t00)
gene_sum_t30_vs_t00_df$gene_name <- rownames(gene_sum_t30_vs_t00_df)

gene_sum_t30_vs_t04_df <- as.data.frame(gene_sum_t30_vs_t04)
gene_sum_t30_vs_t04_df$gene_name <- rownames(gene_sum_t30_vs_t04_df)


# Save the results as a CSV file
write.csv(gene_sum_t04_vs_t00_df, file = "./code/expression/output/DESeq2_gene_sum_t04_vs_t00.csv",
          row.names = F)
write.csv(gene_sum_t30_vs_t00_df, file = "./code/expression/output/DESeq2_gene_sum_t30_vs_t00.csv",
          row.names = F)
write.csv(gene_sum_t30_vs_t04_df, file = "./code/expression/output/DESeq2_gene_sum_t30_vs_t04.csv",
          row.names = F)



######################

###DIFFERENTIAL PROTEIN EXPRESSION ###

#####################

raw_peptides <- read.csv("./data/Raw_pq_3619_peptides_TMTMosaic.csv")


#1) Filter peptides
sn_cols <- grep("^rq_\\d{3}[NC]?_sn$", names(raw_peptides), value = TRUE)

filtered_raw <- raw_peptides[rowSums(raw_peptides[, sn_cols], na.rm = TRUE) > 110, ] #Filtering to total S/N >110 per peptide
filtered_raw <- filtered_raw %>%filter(!str_detect(ProteinId, "##")) #Remove
filtered_raw <- filtered_raw %>%filter(!str_detect(ProteinId, "contaminant")) #Removing contaminants

filtered_raw <- filtered_raw %>% filter(IsolationSpecificity >= 0.5) #TCMP recommends 

# Function to convert TCMP data to MSstatsTMT format
convert_tcmp_to_msstats <- function(raw_peptides) {
  
  # Step 1: Reshape the TMT intensity columns from wide to long format
  # Extract TMT channel columns (rq_*_sn)
  tmt_channels <- grep("^rq_.*_sn$", colnames(raw_peptides), value = TRUE)
  
  long_data <- raw_peptides %>%
    pivot_longer(cols = all_of(tmt_channels), 
                 names_to = "Channel_raw", 
                 values_to = "Intensity") %>%
    # Clean up channel names (remove rq_ prefix and _sn suffix)
    mutate(Channel = gsub("^rq_(.*)_sn$", "\\1", Channel_raw)) %>%
    # Convert channel names to standard format
    mutate(Channel = case_when(
      Channel == "126" ~ "126",
      Channel == "127N" ~ "127N", 
      Channel == "127C" ~ "127C",
      Channel == "128N" ~ "128N",
      Channel == "128C" ~ "128C", 
      Channel == "129N" ~ "129N",
      Channel == "129C" ~ "129C",
      Channel == "130N" ~ "130N",
      Channel == "130C" ~ "130C",
      Channel == "131N" ~ "131N",
      Channel == "131C" ~ "131C",
      TRUE ~ Channel
    )) %>%
    select(-Channel_raw)
  
  # Step 2: Create MSstatsTMT format
  msstats_data <- long_data %>%
    mutate(
      # Use ProteinId as ProteinName
      ProteinName = ProteinId,
      
      # Clean peptide sequence (remove flanking amino acids and format properly)
      PeptideSequence = gsub("^[A-Z]\\.", "", gsub("\\.[A-Z]$", "", PeptideSequence)),
      # Format as [K].SEQUENCE.[A] style
      PeptideSequence = paste0("[", substr(PeptideSequence, 1, 1), "].", 
                               substr(PeptideSequence, 2, nchar(PeptideSequence)-1), 
                               ".[", substr(PeptideSequence, nchar(PeptideSequence), nchar(PeptideSequence)), "]"),
      
      # Create unique PSM identifier 
      PSM = paste(PeptideSequence, Charge, sep = "_"),
      
      # Extract Run information from file path
      Run = gsub(".*/(.*)\\.raw$", "\\1", RunLoadPath),
      
      # Add required MSstatsTMT columns
      Mixture = paste0("Mixture1"), #We had one TMT mixture (which was split into 12 fractions - the 'runs' above)
      TechRepMixture = 1, 
      
      # Remove zero intensities 
      Intensity = ifelse(Intensity == 0, NA, Intensity)
    ) %>%
    # Select and rename columns for MSstatsTMT (in the correct order)
    select(
      ProteinName,
      PeptideSequence, 
      Charge,
      PSM,
      Mixture,
      TechRepMixture,
      Run,
      Channel,
      Intensity
    ) %>%
    # Remove rows with missing intensities
    filter(!is.na(Intensity), Intensity > 0)
  
  return(msstats_data)
}

msstats_formatted <- convert_tcmp_to_msstats(filtered_raw)

#2) Create annotation file using channel mapping
create_annotation_from_channel_map <- function(msstats_data, channel_map) {
  
  channel_map <- data.frame(
    Channel_raw = c("rq_126_sn", "rq_127N_sn", "rq_127C_sn", "rq_128N_sn", "rq_128C_sn",
                    "rq_129N_sn", "rq_129C_sn", "rq_130N_sn", "rq_130C_sn", "rq_131N_sn", "rq_131C_sn"),
    Condition = c("t00", "t4", "t30",
                  "t00", "t4", "t30", 
                  "t00", "t4", "t30",
                  "t00", "t30"),
    BioReplicate = c("R5", "R5", "R5",
                     "R6", "R6", "R6",
                     "R7", "R7", "R7", 
                     "R8", "R8")
  )
  
  # Convert channel names to match MSstatsTMT format
  channel_map <- channel_map %>%
    mutate(Channel = case_when(
      Channel_raw == "rq_126_sn" ~ "126",
      Channel_raw == "rq_127N_sn" ~ "127N",
      Channel_raw == "rq_127C_sn" ~ "127C", 
      Channel_raw == "rq_128N_sn" ~ "128N",
      Channel_raw == "rq_128C_sn" ~ "128C",
      Channel_raw == "rq_129N_sn" ~ "129N", 
      Channel_raw == "rq_129C_sn" ~ "129C",
      Channel_raw == "rq_130N_sn" ~ "130N",
      Channel_raw == "rq_130C_sn" ~ "130C",
      Channel_raw == "rq_131N_sn" ~ "131N",
      Channel_raw == "rq_131C_sn" ~ "131C"
    ))
  
  mixtures <- unique(msstats_data$Mixture)
  
  annotation <- expand.grid(Mixture = mixtures, Channel = channel_map$Channel, TechRepMixture = 1) %>%
    left_join(channel_map %>% select(Channel, Condition, BioReplicate), by = "Channel") %>%
    arrange(Mixture, Channel)
  
  return(annotation)
}

annotation <- create_annotation_from_channel_map(msstats_formatted, channel_map)

#3) Run MSstatsTMT analysis
msstats_data <- msstats_formatted
# Merge annotation info into the data
msstats_with_annotation <- msstats_data %>%
  left_join(annotation, by = c("Mixture", "Channel", "TechRepMixture")) %>%
  filter(!is.na(Condition), !is.na(BioReplicate))

# Check for and aggregate duplicates at PSM level
duplicates <- msstats_with_annotation %>%
  group_by(ProteinName, PeptideSequence, Charge, PSM, Mixture, TechRepMixture, Run, Channel, Condition, BioReplicate) %>%
  summarise(Intensity = sum(Intensity), .groups = "drop")

# Remove any rows with missing information
msstats_with_annotation <- duplicates %>%
  filter(!is.na(Intensity), 
         Intensity > 0,
         !is.na(ProteinName), 
         !is.na(PeptideSequence),
         !is.na(Charge))



# Protein summarization 
quant_data <- proteinSummarization(data = msstats_with_annotation, 
                                   method = "msstats",
                                   global_norm = TRUE,
                                   reference_norm = FALSE,
                                   remove_norm_channel = FALSE,
                                   remove_empty_channel = TRUE,
                                   MBimpute = TRUE)

# Save protein quant
save(quant_data, file = "./code/expression/output/MSstatsTMT_abundance.RData")

write.csv(quant_data[["ProteinLevelData"]], "./code/expression/output/MSstatsTMT_abundance.csv", row.names = F)


# Create contrast matrix for time course comparisons
conditions <- c("t00", "t30", "t4")  # Order matters for contrasts

contrast_matrix <- matrix(c(-1,  1,  0,   # t30 vs t00
                            -1,  0,  1,  # t4 vs t00
                            0,   1, -1), # t30 vs t4
                          nrow = 3, byrow = TRUE)
rownames(contrast_matrix) <- c("t30_vs_t00", "t4_vs_t00", "t30_vs_t4")
colnames(contrast_matrix) <- conditions

test.pairwise <- groupComparisonTMT(quant_data, contrast.matrix = contrast_matrix, moderated = TRUE)

write.csv(test.pairwise$ComparisonResult, "./code/expression/output/MSstatsTMT_diff_prot_expr.csv", row.names = F)


#####################

#CALCULATING GENE (mRNA AND PROTEIN) CLUSTERS

####################

# GENE CLUSTERS (SUMMED TR) -----------------------------------------------
gene_sum_t04_vs_t00_df <- read.csv("./code/expression/output/DESeq2_gene_sum_t04_vs_t00.csv")
gene_sum_t30_vs_t00_df <- read.csv("./code/expression/output/DESeq2_gene_sum_t30_vs_t00.csv")
gene_sum_t30_vs_t04_df <- read.csv("./code/expression/output/DESeq2_gene_sum_t30_vs_t04.csv")


#Make 9 clusters
clusters_gene_sum <- gene_sum_t04_vs_t00_df[,c(2,5:7)]
colnames(clusters_gene_sum)[c(1:3)] <- c("log2FC_t04_v_t00", "pvalue_t04_v_t00", "padj_t04_v_t00")

colnames(gene_sum_t30_vs_t00_df)[c(2,5,6)] <- c("log2FC_t30_v_t00", "pvalue_t30_v_t00", "padj_t30_v_t00")
colnames(gene_sum_t30_vs_t04_df)[c(2,5,6)] <- c("log2FC_t30_v_t04", "pvalue_t30_v_t04", "padj_t30_v_t04")


clusters_gene_sum <- merge(clusters_gene_sum, gene_sum_t30_vs_t00_df[, c("gene_name",
                                                                         "log2FC_t30_v_t00",
                                                                         "pvalue_t30_v_t00",
                                                                         "padj_t30_v_t00")],
                           by.x = "gene_name", by.y = "gene_name", all.x = T)

clusters_gene_sum <- merge(clusters_gene_sum, gene_sum_t30_vs_t04_df[, c("gene_name",
                                                                         "log2FC_t30_v_t04",
                                                                         "pvalue_t30_v_t04",
                                                                         "padj_t30_v_t04")],
                           by.x = "gene_name", by.y = "gene_name", all.x = T)


#Make clusters_gene_sum
clusters_gene_sum <- cbind(clusters_gene_sum, cluster = NA)

clusters_gene_sum <- clusters_gene_sum %>%
  mutate(
    dir1 = case_when(
      log2FC_t04_v_t00 >= 1 & pvalue_t04_v_t00 < 0.05 ~ "U",
      log2FC_t04_v_t00 <= -1 & pvalue_t04_v_t00 < 0.05 ~ "D",
      TRUE ~ "-"
    ),
    dir2 = case_when(
      log2FC_t30_v_t04 >= 1 & pvalue_t30_v_t04 < 0.05 ~ "U",
      log2FC_t30_v_t04 <= -1 & pvalue_t30_v_t04 < 0.05 ~ "D",
      TRUE ~ "-"
    ),
    cluster = paste0(dir1, dir2)
  )
clusters_gene_sum$cluster <- as.factor(clusters_gene_sum$cluster)


write.csv(clusters_gene_sum, "./code/expression/output/gene_clusters_v2.csv", 
          row.names = F)


# PROTEIN CLUSTERS --------------------------------------------------------
diff_prot <- read.csv("./code/expression/output/MSstatsTMT_diff_prot_expr.csv")
diff_prot$UniprotID <- sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", diff_prot$Protein)

#
raw_peptides <- read.csv("./data/Raw_pq_3619_peptides_TMTMosaic.csv")
sn_cols <- grep("^rq_\\d{3}[NC]?_sn$", names(raw_peptides), value = TRUE)

filtered_raw <- raw_peptides[rowSums(raw_peptides[, sn_cols], na.rm = TRUE) > 110, ]
filtered_raw <- filtered_raw %>%filter(!str_detect(ProteinId, "##"))
filtered_raw <- filtered_raw %>%filter(!str_detect(ProteinId, "contaminant"))

filtered_raw <- filtered_raw %>% filter(IsolationSpecificity >= 0.5) 


# Get Gene names:
protein <- read.csv("./data/proteinID_mapping.csv")

diff_prot <- merge(diff_prot, protein[,c("ProteinID", "GeneSymbol")], 
                   by.x = "UniprotID", by.y = "ProteinID", all.x = T)
test <- diff_prot %>% group_by(GeneSymbol) %>% distinct(UniprotID) %>% summarise(n = n()) #some proteins have >1 isoform

raw_peptides_summary <- filtered_raw %>% group_by(ProteinId, PeptideSequence) %>%
  summarise(unique_pep = n()) 
raw_peptides_summary <- raw_peptides_summary %>% summarise(num_unique_pep = sum(unique_pep))
diff_prot <- merge(diff_prot, raw_peptides_summary[,c("ProteinId", "num_unique_pep")], 
                   by.x = "Protein", by.y = "ProteinId", all.x = T)

#Since some proteins have multiple isoforms, 
#Remove the isoform with the lowest number of peptides. 
#If they have the same number, manually remove the non-canonical (according to Uniprot IDs)
prot_ranked <- diff_prot %>%
  group_by(GeneSymbol) %>%
  mutate(max_pep = max(num_unique_pep),
         is_max = num_unique_pep == max_pep) %>%
  filter(is_max) %>%
  mutate(num_ids_with_max = n_distinct(UniprotID)) %>%
  ungroup()

# A. Genes with a single best UniprotID
unique_winners <- prot_ranked %>%
  filter(num_ids_with_max == 1)

# B. Genes with ties — manual review needed
ties_for_manual_review <- prot_ranked %>%
  filter(num_ids_with_max > 1)

manually_selected_ids <- c("Q9BYI3", "O60333-3")  # Add more as needed

manual_picks <- ties_for_manual_review %>%
  filter(UniprotID %in% manually_selected_ids)

final_filtered_diff_prot <- bind_rows(unique_winners, manual_picks)
length(unique(final_filtered_diff_prot$GeneSymbol)) #8878

# write.csv(final_filtered_diff_prot, "./code/expression/output/MSstatsTMT_diff_prot_expr_filtered.csv", 
#           row.names = F)
final_filtered_diff_prot <- read.csv("./code/expression/output/MSstatsTMT_diff_prot_expr_filtered.csv")

wide_df <- final_filtered_diff_prot %>%
  filter(Label %in% c("t4_vs_t00", "t30_vs_t4")) %>%
  select(GeneSymbol, Label, log2FC, adj.pvalue) %>%
  pivot_wider(
    names_from = Label,
    values_from = c(log2FC, adj.pvalue),
    names_glue = "{.value}_{Label}"
  )

# Step 2: Classify direction and cluster
clustered_df <- wide_df %>%
  mutate(
    dir1 = case_when(
      log2FC_t4_vs_t00 >= 1 & adj.pvalue_t4_vs_t00 < 0.05 ~ "U",
      log2FC_t4_vs_t00 <= -1 & adj.pvalue_t4_vs_t00 < 0.05 ~ "D",
      TRUE ~ "-"
    ),
    dir2 = case_when(
      log2FC_t30_vs_t4 >= 1 & adj.pvalue_t30_vs_t4 < 0.05 ~ "U",
      log2FC_t30_vs_t4 <= -1 & adj.pvalue_t30_vs_t4 < 0.05 ~ "D",
      TRUE ~ "-"
    ),
    cluster = paste0(dir1, dir2)
  )


clustered_df$cluster <- as.factor(clustered_df$cluster)

# 
write.csv(clustered_df, "./code/expression/output/protein_clusters_v3.csv", 
          row.names = F)


