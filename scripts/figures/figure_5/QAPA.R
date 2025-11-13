##### QAPA ######
library(dplyr)
library(rtracklayer)
library(DESeq2)
library(stringr)
library(tidyr)
library(arrow)



setwd("")

classification <- read_parquet("./data/final_classification.parquet")


# PREPARE FOR QAPA BUILD --------------------------------------------------
ORFanage_replaced <- rtracklayer::import("./code/IsoformSwitchAnalyzeR/input/ORF_gene_id_replaced_gtf.gtf")
orfanage <- as.data.frame(ORFanage_replaced)
# 
ORFanage_test <- subset(ORFanage_replaced, type %in% c("transcript", "exon", "CDS"))
rtracklayer::export(ORFanage_test, './code/APA/input/orfanage_APA.gtf')

classification_filter <- classification %>% filter(isoform %in% unique(orfanage$transcript_id))

##Making Ensembl ID file. 
Ensembl <- classification_filter[,c(1,7)]
colnames(Ensembl) <- c("Transcript stable ID", "Gene stable ID")
Ensembl$`Gene name` <- Ensembl$`Gene stable ID`
Ensembl$`Gene type` <- "protein_coding"
Ensembl$`Transcript type` <- "protein_coding"

Ensembl <- Ensembl %>% relocate(`Gene stable ID`, .before = `Transcript stable ID`)
Ensembl <- Ensembl %>% relocate(`Gene name`, .after = `Transcript type`)

write.table(Ensembl, "./code/APA/input/Ensembl_custom.txt", 
            row.names = F, 
            quote = F, 
            col.names = T, 
            sep = "\t")


# QAPA BUILD RESULTS ------------------------------------------------------
#After running modified QAPA build function
BED <- read.table("./code/APA/output/QAPA/output_utrs.bed",
                 header = FALSE, sep="\t", stringsAsFactors = FALSE)
BED$cluster <- 1:nrow(BED) #34064
colnames(BED) <- c("Chr", "ExonStart", "ExonEnd", "FullName", "UTR3Length", "Strand", "GeneName", "Cluster")

BED <- BED %>% group_by(GeneName) %>%
  mutate(NumAPA = n()) %>%
  ungroup()
BED_filter <- BED %>% filter(NumAPA > 1)

#Expand BED file for each PB_ID
pattern <- "PB\\.\\d+\\.\\d+"


#Use str_extract_all to extract the matches and unnest them into separate rows
BED_seq <- BED_filter %>%
  mutate(matches = str_extract_all(FullName, pattern)) %>%  # Extract matches as lists
  unnest(matches) 



#Add CPM data for each sample.
cpm <- read.csv("./code/expression/output/transcript_CPM.csv")
BED_seq <- merge(BED_seq, cpm[, c("pb_id", "iPSC_1", "iPSC_2",  "iPSC_3",
                                  "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3",
                                  "CN_1_2" , "CN_1_3" ,"CN_2_1",  "CN_2_2" , "CN_3_1",  "CN_3_2")], 
                by.x = "matches", by.y = "pb_id", all.x = TRUE)

###NEED TO SET A THRESHOLD:
cutoff <- BED_seq %>% group_by(GeneName) %>%
  summarise(iPSC_1_sum = sum(iPSC_1),
            iPSC_2_sum = sum(iPSC_2),
            iPSC_3_sum = sum(iPSC_3),
            NPC_1_1_sum = sum(NPC_1_1), 
            NPC_1_3_sum = sum(NPC_1_3),
            NPC_2_1_sum = sum(NPC_2_1),
            NPC_2_2_sum = sum(NPC_2_2),
            NPC_3_1_sum = sum(NPC_3_1),
            NPC_3_3_sum = sum(NPC_3_3),
            CN_1_2_sum = sum(CN_1_2),
            CN_1_3_sum = sum(CN_1_3),
            CN_2_1_sum = sum(CN_2_1),
            CN_2_2_sum = sum(CN_2_2),
            CN_3_1_sum = sum(CN_3_1),
            CN_3_2_sum = sum(CN_3_2)) %>%
  ungroup()

cutoff$Avgt00 <- rowMeans(subset(cutoff, select = c(2:4)))
cutoff$Avgt04 <- rowMeans(subset(cutoff, select = c(5:10)))
cutoff$Avgt30 <- rowMeans(subset(cutoff, select = c(11:16)))


# T00 VS T30 --------------------------------------------------------------
#Calculate PAU for each cluster for each sample.
cutoff_filter <- cutoff %>%
  filter(Avgt00 >1 & Avgt30 >1) #6109 genes


#Remove the whole gene if it doesn't meet this cutoff. 
BED_seq_t00_v_t30 <- BED_seq[BED_seq$GeneName %in% cutoff_filter$GeneName,] #73,876
length(unique(BED_seq_t00_v_t30$GeneName)) #6109 genes

# 1. Sum counts per GeneName and Cluster
utr_sums <- BED_seq_t00_v_t30 %>%
  group_by(GeneName, Cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 2. Calculate total CPM per GeneName (summing across clusters per gene)
gene_totals <- utr_sums %>%
  group_by(GeneName) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 3. Join totals back to utr_sums and calculate PAU
PAU <- utr_sums %>%
  left_join(gene_totals, by = "GeneName", suffix = c("_cluster", "_gene")) %>%
  mutate(across(ends_with("_cluster"),
                ~ .x / get(gsub("_cluster", "_gene", cur_column())) * 100,
                .names = "{gsub('_cluster', '_PAU', .col)}")) %>%
  # Adjust this select line to include the desired identifiers and PAUs
  dplyr::select(GeneName, Cluster, ends_with("_PAU"))

#Calculate median PAUs.
PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t00 = median(c(iPSC_1_PAU, iPSC_2_PAU,iPSC_3_PAU), na.rm = TRUE))

PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t04 = median(c(NPC_1_1_PAU, NPC_1_3_PAU, NPC_2_1_PAU, NPC_2_2_PAU, NPC_3_1_PAU, NPC_3_3_PAU), na.rm = TRUE))

PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t30 = median(c(CN_1_2_PAU, CN_1_3_PAU, CN_2_1_PAU, CN_2_2_PAU, CN_3_1_PAU, CN_3_2_PAU), na.rm = TRUE))


###Per cluster, select the most proximal PAU (PPAU).
PAU_unique <- BED_seq_t00_v_t30 %>%
  distinct(GeneName, Cluster, .keep_all = T)

PAU <- merge(PAU, PAU_unique[,c('Cluster', 'UTR3Length')])

PPAU <- PAU %>% group_by(GeneName) %>%
  dplyr::slice(which.min(UTR3Length)) %>%
  ungroup()
PPAU$dPPAU <- PPAU$MedianPAU_t00 - PPAU$MedianPAU_t30

PPAU$change <- ifelse(PPAU$dPPAU>20, "Lengthening",
                      ifelse(PPAU$dPPAU < -20, "Shortening", "No change"))


#Export. 
write.csv(PPAU,"./code/APA/output/QAPA/QAPA_PPAU_t00_v_t30.csv", row.names = FALSE)


QAPA_PPAU <- read.csv("./code/APA/output/QAPA/QAPA_PPAU_t00_v_t30.csv")

##Re-run DESeq2, but sum counts per transcripts (per gene) in the QAPA analysis
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)


filtered_tr_count <- merge(tr_count, BED_seq_t00_v_t30[, c("matches", "Cluster", "GeneName")], 
                           by.x = "isoform", by.y = "matches", all.x = F)

cluster_count <- filtered_tr_count %>%
  group_by(GeneName) %>%
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
cluster_count <- as.data.frame(cluster_count)
rownames(cluster_count) <- cluster_count$GeneName
cluster_count <- cluster_count[,-1]

#Get col data
colData <- data.frame(group = c(rep("t00", 3), rep("t04", 6), rep("t30", 6)))
rownames(colData) <- colnames(cluster_count)
colData$group <- as.factor(colData$group)

#Construct DESEqDataSet object
dds_gene <- DESeqDataSetFromMatrix(countData= cluster_count, 
                                   colData=colData, 
                                   design=~ group)

dds_gene <- DESeq(dds_gene)

#Results: cluster-level
cluster_t30_vs_t00 <- results(dds_gene, contrast = c("group", "t30", "t00"))
summary(cluster_t30_vs_t00)

#Save as dfs.
cluster_t30_vs_t00_df <- as.data.frame(cluster_t30_vs_t00)

# Save the results as a CSV file
write.csv(cluster_t30_vs_t00_df, file = "./code/APA/output/QAPA/Deseq2_cluster_t00_v_t30.csv", row.names = T)



# T00 VS T04 --------------------------------------------------------------
#Calculate PAU for each cluster for each sample.
cutoff_filter <- cutoff %>%
  filter(Avgt00 >1 & Avgt04 >1) #5989 genes


#Remove the whole cluster if it doesn't meet this cutoff. 
BED_seq_t00_v_t04 <- BED_seq[BED_seq$GeneName %in% cutoff_filter$GeneName,] #73,233
length(unique(BED_seq_t00_v_t04$GeneName)) #5989 genes

# 1. Sum counts per GeneName and Cluster
utr_sums <- BED_seq_t00_v_t04 %>%
  group_by(GeneName, Cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 2. Calculate total CPM per GeneName (summing across clusters per gene)
gene_totals <- utr_sums %>%
  group_by(GeneName) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 3. Join totals back to utr_sums and calculate PAU
PAU <- utr_sums %>%
  left_join(gene_totals, by = "GeneName", suffix = c("_cluster", "_gene")) %>%
  mutate(across(ends_with("_cluster"),
                ~ .x / get(gsub("_cluster", "_gene", cur_column())) * 100,
                .names = "{gsub('_cluster', '_PAU', .col)}")) %>%
  # Adjust this select line to include the desired identifiers and PAUs
  dplyr::select(GeneName, Cluster, ends_with("_PAU"))

#Calculate median PAUs.
PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t00 = median(c(iPSC_1_PAU, iPSC_2_PAU,iPSC_3_PAU), na.rm = TRUE))

PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t04 = median(c(NPC_1_1_PAU, NPC_1_3_PAU, NPC_2_1_PAU, NPC_2_2_PAU, NPC_3_1_PAU, NPC_3_3_PAU), na.rm = TRUE))

PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t30 = median(c(CN_1_2_PAU, CN_1_3_PAU, CN_2_1_PAU, CN_2_2_PAU, CN_3_1_PAU, CN_3_2_PAU), na.rm = TRUE))


###Per cluster, select the most proximal PAU (PPAU).
PAU_unique <- BED_seq_t00_v_t04 %>%
  distinct(GeneName, Cluster, .keep_all = T)

PAU <- merge(PAU, PAU_unique[,c('Cluster', 'UTR3Length')])

PPAU <- PAU %>% group_by(GeneName) %>%
  dplyr::slice(which.min(UTR3Length)) %>%
  ungroup()
PPAU$dPPAU <- PPAU$MedianPAU_t00 - PPAU$MedianPAU_t04

PPAU$change <- ifelse(PPAU$dPPAU>20, "Lengthening",
                      ifelse(PPAU$dPPAU < -20, "Shortening", "No change"))



#Export. 
write.csv(PPAU,"./code/APA/output/QAPA/QAPA_PPAU_t00_v_t04.csv", row.names = FALSE)


QAPA_PPAU <- read.csv("./code/APA/output/QAPA/QAPA_PPAU_t00_v_t04.csv")


##Re-run DESeq2, but sum counts per cluster. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)


filtered_tr_count <- merge(tr_count, BED_seq_t00_v_t04[, c("matches", "Cluster", "GeneName")], 
                           by.x = "isoform", by.y = "matches", all.x = F)

cluster_count <- filtered_tr_count %>%
  group_by(GeneName) %>%
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
cluster_count <- as.data.frame(cluster_count)
rownames(cluster_count) <- cluster_count$GeneName
cluster_count <- cluster_count[,-1]

#Get col data
colData <- data.frame(group = c(rep("t00", 3), rep("t04", 6), rep("t04", 6)))
rownames(colData) <- colnames(cluster_count)
colData$group <- as.factor(colData$group)

#Construct DESEqDataSet object
library(DESeq2)
dds_gene <- DESeqDataSetFromMatrix(countData= cluster_count, 
                                   colData=colData, 
                                   design=~ group)

dds_gene <- DESeq(dds_gene)

#Results: cluster-level
cluster_t04_vs_t00 <- results(dds_gene, contrast = c("group", "t04", "t00"))
summary(cluster_t04_vs_t00)

#Save as dfs.
cluster_t04_vs_t00_df <- as.data.frame(cluster_t04_vs_t00)

# Save the results as a CSV file
write.csv(cluster_t04_vs_t00_df, file = "./code/APA/output/QAPA/Deseq2_cluster_t00_v_t04.csv", row.names = T)


# T04 VS T30 --------------------------------------------------------------
#Calculate PAU for each cluster for each sample.
cutoff_filter <- cutoff %>%
  filter(Avgt04 >1 & Avgt30 >1) #6262 genes


#Remove the whole cluster if it doesn't meet this cutoff. 
BED_seq_t04_v_t30 <- BED_seq[BED_seq$GeneName %in% cutoff_filter$GeneName,] #75,306
length(unique(BED_seq_t04_v_t30$GeneName)) #6262 genes

# 1. Sum counts per GeneName and Cluster
utr_sums <- BED_seq_t04_v_t30 %>%
  group_by(GeneName, Cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 2. Calculate total CPM per GeneName (summing across clusters per gene)
gene_totals <- utr_sums %>%
  group_by(GeneName) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 3. Join totals back to utr_sums and calculate PAU
PAU <- utr_sums %>%
  left_join(gene_totals, by = "GeneName", suffix = c("_cluster", "_gene")) %>%
  mutate(across(ends_with("_cluster"),
                ~ .x / get(gsub("_cluster", "_gene", cur_column())) * 100,
                .names = "{gsub('_cluster', '_PAU', .col)}")) %>%
  # Adjust this select line to include the desired identifiers and PAUs
  dplyr::select(GeneName, Cluster, ends_with("_PAU"))

#Calculate median PAUs.
PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t00 = median(c(iPSC_1_PAU, iPSC_2_PAU,iPSC_3_PAU), na.rm = TRUE))

PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t04 = median(c(NPC_1_1_PAU, NPC_1_3_PAU, NPC_2_1_PAU, NPC_2_2_PAU, NPC_3_1_PAU, NPC_3_3_PAU), na.rm = TRUE))

PAU <- PAU %>% 
  rowwise() %>% 
  mutate(MedianPAU_t30 = median(c(CN_1_2_PAU, CN_1_3_PAU, CN_2_1_PAU, CN_2_2_PAU, CN_3_1_PAU, CN_3_2_PAU), na.rm = TRUE))


###Per cluster, select the most proximal PAU (PPAU).
PAU_unique <- BED_seq_t04_v_t30 %>%
  distinct(GeneName, Cluster, .keep_all = T)

PAU <- merge(PAU, PAU_unique[,c('Cluster', 'UTR3Length')])

PPAU <- PAU %>% group_by(GeneName) %>%
  dplyr::slice(which.min(UTR3Length)) %>%
  ungroup()
PPAU$dPPAU <- PPAU$MedianPAU_t04 - PPAU$MedianPAU_t30

PPAU$change <- ifelse(PPAU$dPPAU>20, "Lengthening",
                      ifelse(PPAU$dPPAU < -20, "Shortening", "No change"))




#Export. 
write.csv(PPAU,"./code/APA/output/QAPA/QAPA_PPAU_t04_v_t30.csv", row.names = FALSE)



QAPA_PPAU <- read.csv("./code/APA/output/QAPA/QAPA_PPAU_t04_v_t30.csv")
# 

##Re-run DESeq2, but sum counts per cluster. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)


filtered_tr_count <- merge(tr_count, BED_seq_t04_v_t30[, c("matches", "Cluster", "GeneName")], 
                           by.x = "isoform", by.y = "matches", all.x = F)

cluster_count <- filtered_tr_count %>%
  group_by(GeneName) %>%
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
cluster_count <- as.data.frame(cluster_count)
rownames(cluster_count) <- cluster_count$GeneName
cluster_count <- cluster_count[,-1]

#Get col data
colData <- data.frame(group = c(rep("t04", 3), rep("t04", 6), rep("t30", 6)))
rownames(colData) <- colnames(cluster_count)
colData$group <- as.factor(colData$group)

#Construct DESEqDataSet object
library(DESeq2)
dds_gene <- DESeqDataSetFromMatrix(countData= cluster_count, 
                                   colData=colData, 
                                   design=~ group)

dds_gene <- DESeq(dds_gene)

#Results: cluster-level
cluster_t30_vs_t04 <- results(dds_gene, contrast = c("group", "t30", "t04"))
summary(cluster_t30_vs_t04)

#Save as dfs.
cluster_t30_vs_t04_df <- as.data.frame(cluster_t30_vs_t04)

# Save the results as a CSV file
write.csv(cluster_t30_vs_t04_df, file = "./code/APA/output/QAPA/Deseq2_cluster_t04_v_t30.csv", row.names = T)


