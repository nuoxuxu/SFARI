library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(arrow)
library(ggvenn)

# Assuming you have a read_gtf function available in R

# Define functions
read_gtf <- function(file, attributes = c("transcript_id"), keep_attributes = TRUE) {
    library(readr)
    library(dplyr)
    library(stringr)

    # Read the GTF file
    df <- read_tsv(file,
        col_names = c(
            "seqname", "source", "feature", "start", "end",
            "score", "strand", "frame", "attributes"
        ),
        col_types = cols(seqname = col_character()),
        comment = "#",
        show_col_types = FALSE
    )

    # Extract attributes using regex
    for (attribute in attributes) {
        pattern <- paste0(attribute, ' "([^;]*)"')
        df <- df %>%
            mutate(!!attribute := str_extract(attributes, pattern) %>%
                str_replace(paste0(attribute, ' "'), "") %>%
                str_replace('"', ""))
    }

    if (!keep_attributes) {
        df <- df %>% select(-attributes)
    }

    return(df)
}

read_refmap <- function(file) {
    out <- read_tsv(file, col_types = cols()) %>%
        filter(class_code == "=") %>%
        mutate(
            qry_id_list = map(str_split(qry_id_list, ","), ~ map_chr(.x, ~ str_split(.x, "\\|")[[1]][2]))
        ) %>%
        unnest(qry_id_list) %>%
        select(qry_id_list, ref_id) %>%
        rename(transcript_id = qry_id_list)
    out
}

# Read input files
classification <- read_parquet("nextflow_results/V47/final_classification.parquet")
genomic_data_dir <- Sys.getenv("GENOMIC_DATA_DIR", "")
gencode_path <- file.path(genomic_data_dir, "GENCODE/gencode.v47.annotation.gtf")
gencode_gtf <- read_gtf(gencode_path, c("gene_name", "transcript_id", "gene_id")) %>%
    filter(feature == "exon")

TALON_gtf <- read_gtf(
        "nextflow_results/V47/compare/cp_vz_0.75_min_7_recovery_talon_hg38.gtf",
        c("gene_name", "transcript_id", "gene_id", "transcript_status")
    ) %>%
    filter(
        feature == "exon",
        !str_starts(gene_name, "TALON")
    )

classification <- arrow::read_parquet("nextflow_results/V47/final_classification.parquet")

final_transcripts_gtf <- read_gtf("nextflow_results/V47/final_transcripts.gtf", "transcript_id") %>%
    rename(isoform = transcript_id) %>%
    filter(feature == "exon") %>%
    left_join(classification %>% select(isoform, associated_transcript, associated_gene),
        by = "isoform"
    ) %>%
    filter(!str_starts(as.character(associated_gene), "novelGene")) %>%
    mutate(associated_transcript = as.character(associated_transcript))

TALON_ID_to_GENCODE_V47 <- read_refmap("nextflow_results/V47/compare/TALON_GENCODE_V47.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap")

TALON_SFARI <- read_refmap("nextflow_results/V47/compare/TALON_SFARI.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap")

# Start here

# First, identify transcripts that map to multiple ref_ids
multi_mapped_transcripts <- TALON_ID_to_GENCODE_V47 %>%
    group_by(transcript_id) %>%
    summarise(n_mappings = n(), .groups = "drop") %>%
    filter(n_mappings > 1) %>%
    pull(transcript_id)

n_transcripts_removed <- TALON_ID_to_GENCODE_V47 %>%
    filter(transcript_id %in% multi_mapped_transcripts) %>%
    mutate(
        transcript_id_clean = str_split(str_split(transcript_id, "_")[[1]][1], "\\.")[[1]][1],
        ref_id_clean = str_split(ref_id, "\\.")[[1]][1]
    ) %>%
    rowwise() %>%
    mutate(
        transcript_id_clean = str_split(str_split(transcript_id, "_")[[1]][1], "\\.")[[1]][1],
        ref_id_clean = str_split(ref_id, "\\.")[[1]][1]
    ) %>%
    ungroup() %>%
    mutate(is_same = as.integer(transcript_id_clean == ref_id_clean)) %>%
    group_by(transcript_id) %>%
    summarise(is_same_sum = sum(is_same), .groups = "drop") %>%
    filter(is_same_sum == 0) %>%
    nrow()

cat(paste0(n_transcripts_removed, " transcripts will be removed because they are not uniquely mapped to a ref_id and none of the matching ref_id matches the transcript_id\n"))

total_SFARI_tx <- final_transcripts_gtf %>%
    distinct(isoform) %>%
    nrow()

matched_SFARI_tx <- final_transcripts_gtf %>%
    distinct(isoform) %>%
    filter(isoform %in% unique(TALON_SFARI$ref_id)) %>%
    nrow()

SFARI_tx <- final_transcripts_gtf %>% distinct(isoform)
Patowary_tx <- unique(distinct(TALON_gtf, transcript_id))

SFARI_tx <- final_transcripts_gtf %>%
    left_join(select(classification, c(isoform, structural_category)), join_by(isoform == isoform)) %>%
    distinct(isoform, .keep_all = TRUE) %>%
    select(isoform, structural_category)

full_join_df <- TALON_gtf %>%
    distinct(transcript_id, .keep_all = TRUE) %>%
    select(transcript_id, transcript_status) %>%
    left_join(distinct(TALON_SFARI, transcript_id, .keep_all = TRUE), by = "transcript_id") %>%
    full_join(SFARI_tx, by = c("ref_id" = "isoform"), keep = TRUE)


full_join_df %>%
    mutate(
        transcript_id = ifelse(is.na(transcript_id), FALSE, TRUE),
        isoform = ifelse(is.na(isoform), FALSE, TRUE)
    ) %>%
    rename(Patowary = transcript_id, SFARI = isoform) %>%
    ggplot(aes(A = Patowary, B = SFARI)) +
    geom_venn(text_size = 4) +
    theme_void()
ggsave("figures/supplementary/venn_SFARI_Patowary.pdf", width = 4.5, height = 4.5)

full_join_df %>% 
    filter(structural_category!="full-splice_match") %>% 
    mutate(
        overlap = if_else(!is.na(transcript_id), "In Patowary", "Not in Patowary")
    ) %>% 
    group_by(overlap) %>%
    summarise(n = n())

full_join_df %>%
    filter(!is.na(isoform))

