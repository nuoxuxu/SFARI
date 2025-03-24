library(dplyr)
library(readr)
library(arrow)
library(rtracklayer)

transcript_classification <- read_parquet("nextflow_results/V47/final_classification.parquet")

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet")

peptides_gtf <- import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>%
    as_tibble() %>% 
    distinct(transcript_id, .keep_all = TRUE)

classification <- read_tsv("nextflow_results/V47/orfanage/SFARI.protein_classification.tsv")    

SFARI_gene <- read_csv("data/SFARI-Gene_genes_01-13-2025release_03-12-2025export.csv") %>% 
    filter(`gene-score` == 1)
    
peptide_mapping <- peptide_mapping %>% 
    left_join(
        peptides_gtf %>% 
            dplyr::select(transcript_id, detected),
        join_by(peptide == transcript_id)
    ) %>% 
    left_join(
        classification %>% dplyr::select(pb, protein_classification_base), 
        join_by(pb == pb)
    ) %>% 
    left_join(
        transcript_classification %>% dplyr::select(isoform, associated_gene),
        join_by(pb == isoform)
    )

peptide_mapping %>% 
    filter(detected == "True") %>%
    filter(!GENCODE) %>%
    group_by(peptide) %>%
    filter(all(protein_classification_base %in% c("pNIC", "pNNC"))) %>%
    ungroup() %>% 
    distinct(associated_gene, .keep_all = TRUE) %>% 
    filter(associated_gene %in% SFARI_gene$`gene-symbol`) %>% 
    View()
