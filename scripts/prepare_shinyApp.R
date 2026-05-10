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

classification <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
  mutate(
    structural_category = if_else(structural_category %in% c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog"), structural_category, "Other"),
    structural_category = structural_category_labels[structural_category],
    structural_category = factor(structural_category, levels = c("FSM", "ISM", "NIC", "NNC", "Other"))
  )

# Keep per-replicate data (not just means) so the app can show variance
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
  pivot_longer(
    cols = -isoform,
    names_to = "replicate",
    values_to = "log2_cpm"
  ) %>%
  mutate(
    time_point = case_when(
      startsWith(replicate, "iPSC_") ~ "t00",
      startsWith(replicate, "NPC_")  ~ "t04",
      startsWith(replicate, "CN_")   ~ "t30"
    ),
    time_point = factor(time_point, levels = c("t30", "t04", "t00"))
  ) %>%
  left_join(
    classification[, c("isoform", "structural_category", "associated_gene")],
    by = "isoform"
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

# Flag transcripts that have CDS features (i.e. predicted protein-coding)
gtf <- gtf %>%
  mutate(
    has_CDS = case_when(!transcript_id %in% tx_no_CDS ~ TRUE, .default = FALSE)
  )

# Propagate has_CDS flag to expression data
cds_flag <- gtf %>%
  distinct(transcript_id, has_CDS) %>%
  dplyr::rename(isoform = transcript_id)

lr_log2_cpm <- lr_log2_cpm %>%
  left_join(cds_flag, by = "isoform")

dir.create("transcript_vis_app/data", showWarnings = FALSE)
lr_log2_cpm %>% write.csv("transcript_vis_app/data/lr_log2_cpm.csv", row.names = FALSE)
gtf %>% write.csv("transcript_vis_app/data/gtf.csv", row.names = FALSE)