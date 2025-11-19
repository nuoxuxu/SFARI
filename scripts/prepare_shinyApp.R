library(arrow)
library(dplyr)
library(rtracklayer)
library(tidyr)

structural_category_labels <- c(
    "full-splice_match"        = "FSM",
    "incomplete-splice_match"  = "ISM",
    "novel_in_catalog"         = "NIC",
    "novel_not_in_catalog"     = "NNC",
    "Other"                    = "Other"
)

time_point_labels <- c(
  "mean_iPSC" = "t00",
  "mean_NPC"  = "t04",
  "mean_CN"   = "t30"
)

classification <- read_parquet("nextflow_results/V47/final_classification.parquet") %>% 
  mutate(
    structural_category = if_else(structural_category %in% c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"), structural_category, "Other"),
    structural_category = structural_category_labels[structural_category],
    structural_category = factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC", "Other"))
  )

lr_log2_cpm <- read_parquet("nextflow_results/V47/final_expression.parquet") %>%
  dplyr::rename(
    NPC_1_2 = NPC_1_3,
    NPC_3_2 = NPC_3_3,
    CN_1_1 = CN_1_2,
    CN_1_2 = CN_1_3
  ) %>%
  mutate(
    across(
      where(is.numeric), 
      ~ log2((.x / sum(.x) * 1e6) + 1)
    )
  ) %>% 
  mutate(
    mean_iPSC = rowMeans(select(., starts_with("iPSC_")), na.rm = TRUE),
    mean_NPC  = rowMeans(select(., starts_with("NPC_")), na.rm = TRUE),
    mean_CN   = rowMeans(select(., starts_with("CN_")), na.rm = TRUE)
  ) %>% 
  select(isoform, mean_iPSC, mean_NPC, mean_CN) %>% 
  pivot_longer(
    cols = c("mean_iPSC", "mean_NPC", "mean_CN"),
    names_to = "time_point",
    values_to = "abundance"
  ) %>% 
  left_join(
    classification[, c("isoform", "structural_category", "associated_gene")],
    by = "isoform"
  ) %>%
  mutate(
    time_point = time_point_labels[time_point],
    time_point = factor(time_point, levels = c("t30", "t04", "t00"))
  )

orfanage_gtf <- rtracklayer::import("nextflow_results/V47/orfanage/orfanage.gtf") %>% as.data.frame()
full_gtf <- rtracklayer::import("nextflow_results/V47/final_transcripts.gtf") %>% as.data.frame()
tx_no_CDS <- setdiff(pull(distinct(full_gtf, transcript_id), transcript_id), pull(distinct(orfanage_gtf, transcript_id), transcript_id))
gtf <- bind_rows(orfanage_gtf, full_gtf %>% filter(transcript_id %in% tx_no_CDS))

gtf <- gtf %>%
  left_join(
    classification[, c("isoform", "structural_category", "associated_gene")],
    by = c("transcript_id" = "isoform")
  )

gtf <- gtf %>%
  select(
    -c(width, source, score, phase, gene_id, ID, orfanage_duplicity, orfanage_status, orfanage_template, orfanage_template_source, Parent, original_biotype)
  ) %>%
  filter(!(type %in% c("five_prime_utr", "three_prime_utr", "start_codon", "stop_codon", "gene")))

lr_log2_cpm %>% write.csv("transcript_vis_app/data/lr_log2_cpm.csv", row.names = FALSE)
gtf %>% write.csv("transcript_vis_app/data/gtf.csv", row.names = FALSE)