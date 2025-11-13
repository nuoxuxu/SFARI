### APA ###
library(rtracklayer)
library(seqinr)
library(Biostrings)
library(dplyr)
library(tidyr)
library(arrow)
library(edgeR)
library(stringr)
library(purrr)


setwd("")

classification <- read_parquet("./data/final_classification.parquet") #182371 

# PREP FILEs -------------------------------------------------
#Read in transcript regions
ORFanage_replaced <- rtracklayer::import("./code/IsoformSwitchAnalyzeR/input/ORF_gene_id_replaced_gtf.gtf") #only using 158,844 here
gtf_df <- as.data.frame(ORFanage_replaced)

# Keep only transcripts with a 3'UTR
transcripts_with_3utr <- gtf_df %>%
  group_by(transcript_id) %>%
  filter(any(type == "three_prime_utr")) %>%
  ungroup()

gtf_filtered <- ORFanage_replaced[ORFanage_replaced$transcript_id %in% transcripts_with_3utr$transcript_id] #158,841 remaining


# CHECKING FOR IDENTICAL SEQUENCES - using tr structure ----------------------
# 1. Keep only relevant features
features_keep <- gtf_filtered[gtf_filtered$type %in% c("five_prime_utr", "CDS", "stop_codon")]

# 2. Convert to tibble for easier manipulation
features_df <- as_tibble(features_keep)

# 3. For each transcript, sort features by start and create a "signature" string
transcript_signatures <- features_df %>%
  arrange(transcript_id, start) %>%
  group_by(gene_id, transcript_id, strand) %>%
  summarize(
    signature = paste0(type, ":", start, "-", end, collapse = "|"),
    .groups = "drop"
  )

# 4. Find identical transcripts within each gene
identical_transcripts <- transcript_signatures %>%
  group_by(gene_id, signature) %>%
  summarize(
    transcripts = list(transcript_id),
    n_transcripts = n(),
    .groups = "drop"
  )


# Step 1: Assign a sequence_cluster ID to each unique signature
identical_transcripts <- identical_transcripts %>%
  mutate(sequence_cluster = paste0("cluster_", row_number()))

# Step 2: Uncollapse all transcripts so each row has transcript_id + cluster
transcript_cluster_df <- identical_transcripts %>%
  dplyr::select(sequence_cluster, transcripts) %>%
  unnest(cols = transcripts)

# Add gene_id back
transcript_cluster_df <- transcript_cluster_df %>%
  left_join(
    transcript_signatures %>% select(transcript_id, gene_id),
    by = c("transcripts" = "transcript_id")
  ) %>%
  rename(transcript_id = transcripts)

utr3_gr <- gtf_filtered[gtf_filtered$type == "three_prime_utr"]

# Step 2: Compute 3'UTR length per transcript
utr3_lengths <- as_tibble(utr3_gr) %>%
  group_by(transcript_id) %>%
  summarize(
    utr3_length = sum(width),
    .groups = "drop"
  )

# Step 3: Merge with your transcript cluster table
transcript_cluster_df <- transcript_cluster_df %>%
  left_join(utr3_lengths, by = c("transcript_id"))

# Now each transcript has: transcript_id, gene_id, sequence_cluster, utr3_length
head(transcript_cluster_df)

transcript_cluster_df <- as.data.frame(transcript_cluster_df)
ident_seq <- transcript_cluster_df %>% 
  group_by(sequence_cluster) %>%
  mutate(UTR3_length_rel_to_proximal = utr3_length - min(utr3_length)) %>%
  ungroup()

# Function to assign dynamic bins per group (group pA sites by within 24 nt)
assign_dynamic_bins <- function(values, bin_size = 24) {
  sorted_values <- sort(values)
  groups <- numeric(length(values))
  current_bin_end <- -1
  group_number <- 1
  
  for (i in seq_along(sorted_values)) {
    val <- sorted_values[i]
    if (val > current_bin_end) {
      current_bin_end <- val + bin_size - 1
      group_number <- group_number + 1
    }
    bin_range <- val <= current_bin_end
    groups[values == val & groups == 0] <- group_number - 1
  }
  
  return(groups)
}

# Apply per group
ident_seq <- ident_seq %>%
  group_by(sequence_cluster) %>%
  arrange(UTR3_length_rel_to_proximal, .by_group = TRUE) %>%
  mutate(UTR3_length_cluster = assign_dynamic_bins(UTR3_length_rel_to_proximal)) %>%
  ungroup()

#need to have at least 2 UTR lengths per sequence cluster. 
ident_seq <- ident_seq %>% group_by(sequence_cluster) %>%
  filter(n_distinct(UTR3_length_cluster) > 1) %>%
  ungroup() #14,738 transcripts

###Sum each sequence cluster, and filter by expression. 
cpm <- read.csv("./code/expression/output/transcript_CPM.csv")

cpm$Avgt00 <- rowMeans(subset(cpm, select = c(1:3)))
cpm$Avgt04 <- rowMeans(subset(cpm, select = c(4:9)))
cpm$Avgt30 <- rowMeans(subset(cpm, select = c(10:15)))


#Add CPM data for each sample.
ident_seq <- merge(ident_seq, cpm[, c("pb_id", "iPSC_1", "iPSC_2",  "iPSC_3",
                                      "NPC_1_1", "NPC_1_3", "NPC_2_1", "NPC_2_2", "NPC_3_1", "NPC_3_3",
                                      "CN_1_2" , "CN_1_3" ,"CN_2_1",  "CN_2_2" , "CN_3_1",  "CN_3_2")], 
                   by.x = "transcript_id", by.y = "pb_id", all.x = TRUE)

#Summarise to the sum per cluster. 
#This will allow for some transcripts with 0s. 
cutoff <- ident_seq %>% group_by(sequence_cluster) %>%
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




##############
###T00 vs T30!
cutoff_filter <- cutoff %>%
  filter((Avgt00 >1 & Avgt30 >1)) #2879 clusters


#Remove the whole cluster if it doesn't meet this cutoff. 
ident_seq_t00_t30 <- ident_seq
ident_seq_t00_t30 <- ident_seq_t00_t30[ident_seq_t00_t30$sequence_cluster %in% cutoff_filter$sequence_cluster,] 
length(unique(ident_seq_t00_t30$gene_id)) #2044

#Calculate PAU for each cluster for each sample.
# 1. Sum CPMs per UTR3 per sequence_cluster
utr_sums <- ident_seq_t00_t30 %>%
  group_by(sequence_cluster, UTR3_length_cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 2. Total CPM sums per sequence_cluster
cluster_totals <- utr_sums %>%
  group_by(sequence_cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 3. Join totals back to UTRs and calculate PAUs
PAU <- utr_sums %>%
  left_join(cluster_totals, by = "sequence_cluster", suffix = c("_utr", "_total")) %>%
  mutate(across(ends_with("_utr"), 
                ~ .x / get(gsub("_utr", "_total", cur_column())) * 100, 
                .names = "{gsub('_utr', '_PAU', .col)}")) %>%
  dplyr::select(sequence_cluster, UTR3_length_cluster, ends_with("_PAU"))



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
PPAU <- PAU %>% group_by(sequence_cluster) %>%
  dplyr::slice(which.min(UTR3_length_cluster)) %>%
  ungroup()
PPAU$dPPAU <- PPAU$MedianPAU_t00 - PPAU$MedianPAU_t30

PPAU$change <- ifelse(PPAU$dPPAU>20, "Lengthening",
                      ifelse(PPAU$dPPAU < -20, "Shortening", "No change"))

#Bring back some info:
PAU_unique <- ident_seq_t00_t30 %>%
  distinct(sequence_cluster, UTR3_length_cluster, .keep_all = T)
PPAU <- merge(PPAU, PAU_unique[,c('sequence_cluster', "gene_id", "UTR3_length_cluster")],
              by.y = c("sequence_cluster", "UTR3_length_cluster"), by.x = c("sequence_cluster", "UTR3_length_cluster"), 
              all.y = F)

#Export. 
write.csv(PPAU,"./code/APA/output/APA/PPAU_t00_v_t30.csv", row.names = FALSE)



K_PPAU <- read.csv("./code/APA/output/APA/PPAU_t00_v_t30.csv")
# 

##Re-run DESeq2, but sum counts per cluster. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)


filtered_tr_count <- merge(tr_count, PAU_unique[, c("transcript_id", "sequence_cluster")], 
                           by.x = "isoform", by.y = "transcript_id", all.x = F) #keep only 6961

cluster_count <- filtered_tr_count %>%
  group_by(sequence_cluster) %>%
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
rownames(cluster_count) <- cluster_count$sequence_cluster
cluster_count <- cluster_count[,-1]

#Get col data
colData <- data.frame(group = c(rep("t00", 3), rep("t04", 6), rep("t30", 6)))
rownames(colData) <- colnames(cluster_count)
colData$group <- as.factor(colData$group)

#Construct DESEqDataSet object
library(DESeq2)
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
write.csv(cluster_t30_vs_t00_df, file = "./code/APA/output/APA/Deseq2_cluster_lrmethod_t00_v_t30.csv", row.names = T)




####T00 vs T04:
cutoff_filter <- cutoff %>%
  filter((Avgt00 >1 & Avgt04 >1)) #2914 clusters


#Remove the whole cluster if it doesn't meet this cutoff. 
ident_seq_t00_t04 <- ident_seq
ident_seq_t00_t04 <- ident_seq_t00_t04[ident_seq_t00_t04$sequence_cluster %in% cutoff_filter$sequence_cluster,] 
length(unique(ident_seq_t00_t04$gene_id)) #2041

#Calculate PAU for each cluster for each sample.
# 1. Sum CPMs per UTR3 per sequence_cluster
utr_sums <- ident_seq_t00_t04 %>%
  group_by(sequence_cluster, UTR3_length_cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 2. Total CPM sums per sequence_cluster
cluster_totals <- utr_sums %>%
  group_by(sequence_cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 3. Join totals back to UTRs and calculate PAUs
PAU <- utr_sums %>%
  left_join(cluster_totals, by = "sequence_cluster", suffix = c("_utr", "_total")) %>%
  mutate(across(ends_with("_utr"), 
                ~ .x / get(gsub("_utr", "_total", cur_column())) * 100, 
                .names = "{gsub('_utr', '_PAU', .col)}")) %>%
  dplyr::select(sequence_cluster, UTR3_length_cluster, ends_with("_PAU"))



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
PPAU <- PAU %>% group_by(sequence_cluster) %>%
  dplyr::slice(which.min(UTR3_length_cluster)) %>%
  ungroup()
PPAU$dPPAU <- PPAU$MedianPAU_t00 - PPAU$MedianPAU_t04

PPAU$change <- ifelse(PPAU$dPPAU>20, "Lengthening",
                      ifelse(PPAU$dPPAU < -20, "Shortening", "No change"))

#Bring back some info:
PAU_unique <- ident_seq_t00_t04 %>%
  distinct(sequence_cluster, UTR3_length_cluster, .keep_all = T)
PPAU <- merge(PPAU, PAU_unique[,c('sequence_cluster', "gene_id", "UTR3_length_cluster")],
              by.y = c("sequence_cluster", "UTR3_length_cluster"), by.x = c("sequence_cluster", "UTR3_length_cluster"), 
              all.y = F)

#Export. 
write.csv(PPAU,"./code/APA/output/APA/PPAU_t00_v_t04.csv", row.names = FALSE)



K_PPAU <- read.csv("./code/APA/output/APA/PPAU_t00_v_t04.csv")
# 

##Re-run DESeq2, but sum counts per cluster. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)


filtered_tr_count <- merge(tr_count, PAU_unique[, c("transcript_id", "sequence_cluster")], 
                           by.x = "isoform", by.y = "transcript_id", all.x = F) #keep only 7033

cluster_count <- filtered_tr_count %>%
  group_by(sequence_cluster) %>%
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
rownames(cluster_count) <- cluster_count$sequence_cluster
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
write.csv(cluster_t04_vs_t00_df, file = "./code/APA/output/APA/Deseq2_cluster_lrmethod_t00_v_t04.csv", row.names = T)




###T04 vs T30:
cutoff_filter <- cutoff %>%
  filter((Avgt04 >1 & Avgt30 >1)) #3325 clusters


#Remove the whole cluster if it doesn't meet this cutoff. 
ident_seq_t04_t30 <- ident_seq
ident_seq_t04_t30 <- ident_seq_t04_t30[ident_seq_t04_t30$sequence_cluster %in% cutoff_filter$sequence_cluster,] 
length(unique(ident_seq_t04_t30$gene_id)) #2259

#Calculate PAU for each cluster for each sample.
# 1. Sum CPMs per UTR3 per sequence_cluster
utr_sums <- ident_seq_t04_t30 %>%
  group_by(sequence_cluster, UTR3_length_cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 2. Total CPM sums per sequence_cluster
cluster_totals <- utr_sums %>%
  group_by(sequence_cluster) %>%
  summarise(across(iPSC_1:CN_3_2, ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 3. Join totals back to UTRs and calculate PAUs
PAU <- utr_sums %>%
  left_join(cluster_totals, by = "sequence_cluster", suffix = c("_utr", "_total")) %>%
  mutate(across(ends_with("_utr"), 
                ~ .x / get(gsub("_utr", "_total", cur_column())) * 100, 
                .names = "{gsub('_utr', '_PAU', .col)}")) %>%
  dplyr::select(sequence_cluster, UTR3_length_cluster, ends_with("_PAU"))



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
PPAU <- PAU %>% group_by(sequence_cluster) %>%
  dplyr::slice(which.min(UTR3_length_cluster)) %>%
  ungroup()
PPAU$dPPAU <- PPAU$MedianPAU_t04 - PPAU$MedianPAU_t30

PPAU$change <- ifelse(PPAU$dPPAU>20, "Lengthening",
                      ifelse(PPAU$dPPAU < -20, "Shortening", "No change"))

#Bring back some info:
PAU_unique <- ident_seq_t04_t30 %>%
  distinct(sequence_cluster, UTR3_length_cluster, .keep_all = T)
PPAU <- merge(PPAU, PAU_unique[,c('sequence_cluster', "gene_id", "UTR3_length_cluster")],
              by.y = c("sequence_cluster", "UTR3_length_cluster"), by.x = c("sequence_cluster", "UTR3_length_cluster"), 
              all.y = F)

#Export. 
write.csv(PPAU,"./code/APA/output/APA/PPAU_t04_v_t30.csv", row.names = FALSE)



K_PPAU <- read.csv("./code/APA/output/APA/PPAU_t04_v_t30.csv")
# 

##Re-run DESeq2, but sum counts per cluster. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)


filtered_tr_count <- merge(tr_count, PAU_unique[, c("transcript_id", "sequence_cluster")], 
                           by.x = "isoform", by.y = "transcript_id", all.x = F) #keep only 8008

cluster_count <- filtered_tr_count %>%
  group_by(sequence_cluster) %>%
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
rownames(cluster_count) <- cluster_count$sequence_cluster
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
write.csv(cluster_t30_vs_t04_df, file = "./code/APA/output/APA/Deseq2_cluster_lrmethod_t04_v_t30.csv", row.names = T)

