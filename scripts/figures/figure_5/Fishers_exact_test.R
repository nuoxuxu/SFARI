### AS CORRELATION EVENTS - FISHER'S EXACT TEST APPROACH ####
library(arrow)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(tibble)


setwd("")

#Read in data:
ES_events <- read.table("./code/AS_APA/output/output_APA_AS_corr/ORFanage_events_SE_strict.ioe", 
                        header = T) #From SUPPA2
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
counts_matrix <- counts_matrix[,c(1, 17:19)]

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

#Filter so each pA site is >= 100 counts (included + skipped) at EACH timepoint
# 1. Sum counts per event_id, timepoint, and pA
pA_totals <- counts_annotated %>%
  group_by(event_id, timepoint, pas_id) %>%
  summarise(total = sum(count), .groups = "drop")

# 2. Keep only pAs with total >=100 for ALL timepoints
pA_keep <- pA_totals %>%
  group_by(event_id, pas_id) %>%
  filter(all(total >= 100)) %>% ##setting this high!
  ungroup()

# 3. Remove events with <2 pAs left
pA_keep_valid <- pA_keep %>%
  group_by(event_id) %>%
  filter(n_distinct(pas_id) >= 2) %>%
  ungroup()

# 4. Filter original counts_annotated down to valid events/pAs
counts_filtered <- counts_annotated %>%
  semi_join(pA_keep_valid, by = c("event_id", "pas_id"))

# Summarise counts first
counts_grouped <- counts_filtered %>%
  group_by(timepoint, event_id, pas_id, type) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  ungroup()

counts_wide <- counts_grouped %>%
  pivot_wider(
    id_cols = c(timepoint, event_id, pas_id),
    names_from = type,
    values_from = total,
    values_fill = 0   # fill missing Included/Skipped with 0
  )

#Remove cases where included happens < 100 (total) or skipped happens < 100
counts_filtered <- counts_wide %>%
  group_by(timepoint, event_id) %>%
  filter(sum(Included) >= 100 & sum(Skipped) >= 100) %>%
  ungroup()

#Remove the entire event_id if it's not at all 3 timepoint
events_all_timepoints <- counts_filtered %>%
  distinct(timepoint, event_id) %>%
  group_by(event_id) %>%
  filter(n() == 3) %>%  
  pull(event_id)

# Filter the counts
counts_filtered_all3 <- counts_filtered %>%
  filter(event_id %in% events_all_timepoints)

#Add pseudocount of 1:
counts_wide_pseudo <- counts_filtered_all3 %>%
  mutate(
    Included = Included + 1,
    Skipped  = Skipped + 1
  )

#Pivot back to long: 
counts_long_ready <- counts_wide_pseudo %>%
  pivot_longer(
    cols = c(Included, Skipped),
    names_to = "type",
    values_to = "count"
  )


# Identify valid event_ids: must have both types at all timepoints
valid_events <- counts_long_ready %>%
  group_by(timepoint, event_id) %>%
  summarise(n_types = n_distinct(type), .groups = "drop") %>%
  group_by(event_id) %>%
  summarise(all_timepoints_ok = all(n_types == 2), .groups = "drop") %>%
  filter(all_timepoints_ok) %>%
  pull(event_id)

# Keep only valid events across all timepoints
counts_grouped_filtered <- counts_long_ready %>%
  filter(event_id %in% valid_events)


# Build named nested lists
tables_by_time <- counts_grouped_filtered %>%
  group_split(timepoint) %>%
  set_names(map_chr(., ~ unique(.x$timepoint))) %>%
  map(~ {
    .x %>%
      group_split(event_id) %>%
      set_names(map_chr(., ~ unique(.x$event_id))) %>%
      map(~ {
        tab <- pivot_wider(.x, names_from = type, values_from = count, values_fill = 0)
        mat <- as.data.frame(tab) %>%
          column_to_rownames("pas_id") %>%
          as.matrix()
        mat
      })
  })


###Apply Fisher's:
# Clean a character matrix: keep only Included/Skipped and convert to numeric
clean_matrix <- function(mat) {
  mat_num <- apply(mat[, c("Included", "Skipped")], 2, function(x) as.numeric(trimws(x)))
  rownames(mat_num) <- rownames(mat)
  mat_num
}

# Fisher test for one table
run_fisher <- function(mat_num, timepoint, event_id) {
  n <- nrow(mat_num)
  results <- vector("list", n)
  
  for (i in seq_len(n)) {
    this <- mat_num[i, ]
    others <- colSums(mat_num[-i, , drop = FALSE])
    table2x2 <- rbind(this, others)
    
    # Skip if table contains negative or NA
    if (any(is.na(table2x2)) || any(table2x2 < 0)) next
    
    ft <- fisher.test(table2x2)
    results[[i]] <- data.frame(
      timepoint = timepoint,
      event_id = event_id,
      pA_id = rownames(mat_num)[i],
      OR = unname(ft$estimate),
      pval = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
  
  bind_rows(results)
}

# Apply across nested list
results_list <- imap_dfr(tables_by_time, function(events, tp) {
  imap_dfr(events, function(mat, ev) {
    mat_num <- clean_matrix(mat)
    run_fisher(mat_num, timepoint = tp, event_id = ev)
  })
})

# Global p-value correction
results_list <- results_list %>%
  mutate(pval_adj = p.adjust(pval, method = "BH"))
results_list <- results_list %>%
  mutate(logOR = log(OR + 0.0000000000001)) #Add tiny value here 

write.csv(results_list, "./code/AS_APA/output/Fisher_result.csv", row.names = F)



## Identify significant events:
results_list <- read.csv("./code/AS_APA/output/Fisher_result.csv")

results_list <- results_list %>%
  mutate(
    gene_id = sub(";.*", "", event_id) )

sig_results <- results_list %>%
  group_by(timepoint, event_id) %>%
  # keep event_ids where at least one pA has p_adj < 0.05 AND |logOR| >=1
  filter(any(pval_adj < 0.05 & abs(logOR >= 1)) ) %>% 
  ungroup()


#OR = <1 = more skipped 
#OR = >1 = more included 
