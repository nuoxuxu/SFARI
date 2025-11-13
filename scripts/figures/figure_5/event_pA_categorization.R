### AS CORRELATION EVENTS ####
library(arrow)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

setwd("")


###Load in data:
classification <- read_parquet("./data/final_classification.parquet")
# 

# Read in exon splicing events (from running SUPPA2)
ES_events <- read.table("./code/AS_APA/output/output_APA_AS_corr/ORFanage_events_SE_strict.ioe", 
                        header = T)

# Read in all transcripts, and group into pA sites
ORFanage_replaced <- rtracklayer::import("./code/IsoformSwitchAnalyzeR/input/ORF_gene_id_replaced_tr_exon.gtf") #only using 158,844 here

exons <- ORFanage_replaced[ORFanage_replaced$type == "exon"]
transcripts <- ORFanage_replaced[ORFanage_replaced$type == "transcript"]


#Get tr end positions
transcript_df <- as.data.frame(transcripts)

# Strand-aware transcript end site
transcript_df <- transcript_df %>%
  mutate(transcript_end = if_else(strand == "+", end, start)) %>%
  dplyr::select(transcript_id, gene_id, transcript_end, strand)

cluster_transcript_ends <- function(df, window_size = 24) {
  df %>%
    mutate(
      transcript_end = as.numeric(transcript_end),
      strand = as.character(strand)
    ) %>%
    arrange(gene_id, strand, transcript_end * ifelse(strand == "+", 1, -1)) %>%
    group_by(gene_id, strand) %>%
    mutate(cluster_id = {
      n_trans <- n()
      
      if (n_trans == 1) {
        "pA1"
      } else {
        cluster <- integer(n_trans)
        current_cluster <- 1L
        cluster[1] <- current_cluster
        start_pos <- transcript_end[1]
        
        for (i in 2:n_trans) {
          dist <- if (strand[1] == "+") {
            transcript_end[i] - start_pos
          } else {
            start_pos - transcript_end[i]
          }
          
          if (!is.na(dist) && dist <= window_size) {
            cluster[i] <- current_cluster
          } else {
            current_cluster <- current_cluster + 1L
            cluster[i] <- current_cluster
            start_pos <- transcript_end[i]
          }
        }
        paste0("pA", cluster)
      }
    }) %>%
    ungroup() %>%
    group_by(gene_id, strand, cluster_id) %>%
    mutate(
      cluster_proximal_coord = if (cur_group()$strand == "+") {
        min(transcript_end, na.rm = TRUE)
      } else {
        max(transcript_end, na.rm = TRUE)
      }
    ) %>%
    ungroup()
}


clustered_df <- cluster_transcript_ends(transcript_df)

#Remove transcripts with only 1 pA site:
filtered_df <- clustered_df %>%
  group_by(gene_id) %>%
  filter(n_distinct(cluster_id) > 1) %>%
  ungroup()
filtered_tr <- unique(filtered_df$transcript_id) #138,308



# Only keep relevant transcripts per splicing event (i.e. Remove all transcripts not present in the 'total_transcripts' col.)
all_transcripts <- unlist(strsplit(ES_events$total_transcripts, ","))
all_transcripts <- unique(all_transcripts) #127,451

#Now filter
filtered_df_ES <- filtered_df %>%
  filter(transcript_id %in% all_transcripts) #113,400

#Re-do the filtering to remove genes w/ only 1 pA site
filtered_df_ES <- filtered_df_ES %>%
  group_by(gene_id) %>%
  filter(n_distinct(cluster_id) > 1) %>%
  ungroup() #111,744


write.csv(filtered_df_ES, "./code/AS_APA/input/transcript_pA_sites.csv", 
          row.names = F)


# Generate counts matrix  --------------------------------------------------
filtered_df_ES <- read.csv("./code/AS_APA/input/transcript_pA_sites.csv")


#Load in data
final_pb_ids <- read.table("./data/filtered_transcript_list.txt", 
                           header = T) #only using transcripts with an ORFanage-predicted ORF from FSM, ISM, NIC or NNC categories

#Subset the transcript file to these PB IDs. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)

filtered_tr_count <- tr_count[tr_count$isoform %in% final_pb_ids$x, ]
filtered_tr_count$sumt00 <- rowSums(subset(filtered_tr_count, select = c(2:4)))
filtered_tr_count$sumt04 <- rowSums(subset(filtered_tr_count, select = c(5:10)))
filtered_tr_count$sumt30 <- rowSums(subset(filtered_tr_count, select = c(11:16)))

colnames(filtered_tr_count)[1] <- "transcript_id"
colnames(filtered_df_ES)[5] <- "pas_id"

#Filter counts:
counts_matrix <- filtered_tr_count %>% filter(transcript_id %in% filtered_df_ES$transcript_id)
counts_matrix <- counts_matrix[,c(1:16)]

# ---- 1. Reshape counts into long format ----
counts_long <- counts_matrix %>%
  as.data.frame() %>%
  pivot_longer(-transcript_id, names_to = "timepoint", values_to = "count")

# ---- 2. Add PAS info ----
counts_long <- counts_long %>%
  left_join(filtered_df_ES, by = "transcript_id")
# now we have: transcript_id | timepoint | count | pA_id

# Expand total_transcripts to long
tot_map <- ES_events %>%
  dplyr::select(event_id, total_transcripts, alternative_transcripts) %>%
  separate_rows(total_transcripts, sep = ",") %>%
  dplyr::rename(transcript_id = total_transcripts) %>%
  rowwise() %>%
  mutate(type = ifelse(transcript_id %in% strsplit(alternative_transcripts, ",")[[1]],
                       "Included", "Skipped")) %>%
  ungroup() %>%
  dplyr::select(event_id, transcript_id, type)

event_map <- tot_map


# ---- 4. Annotate counts with event info ----
counts_annotated <- counts_long %>%
  inner_join(event_map, by = "transcript_id")
# columns: transcript_id | timepoint | count | pA_id | event_id | type
write.csv(counts_annotated, "./code/AS_APA/input/tr_pA_counts.csv", row.names = F)


# 1. Sum counts per event_id, timepoint, and pA
pA_totals <- counts_annotated %>%
  group_by(event_id, type, pas_id, timepoint) %>%
  summarise(total = sum(count), .groups = "drop") #6619 genes

write.csv(pA_totals, "./code/AS_APA/input/pA_counts.csv", row.names = F)


