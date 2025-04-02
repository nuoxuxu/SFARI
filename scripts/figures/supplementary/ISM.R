library(ggplot2)
library(dplyr)
library(arrow)
library(rtracklayer)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

gencode_gtf <- import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf")) %>%
    as_tibble() %>% 
    filter(type=="transcript")

# ISM transcript gene type percentages

read_parquet("nextflow_results/V47/final_classification.parquet") %>% 
    filter(structural_category == "incomplete-splice_match") %>% 
    left_join(
        gencode_gtf %>% select(gene_type, transcript_id),
        join_by(associated_transcript==transcript_id)
    ) %>% 
    mutate(gene_type = if_else(gene_type %in% c("unprocessed_pseudogene", "transcribed_unprocessed_pseudogene", "transcribed_processed_pseudogene"), "pseudogene", gene_type)) %>%
    filter(
        gene_type %in% c("lncRNA", "protein_coding", "pseudogene")
    ) %>%
    mutate(gene_type = factor(gene_type, levels = c("protein_coding", "lncRNA", "pseudogene"))) %>%
    ggplot(aes(x=gene_type, fill=gene_type)) +
        geom_bar(fill = "#0072B2") +
        scale_fill_brewer(palette = "Set2") +
        scale_y_continuous(labels = function(x) x / 1000) +
        labs(
            x = "",
            y = expression("Transcripts (x" ~ 10^3 * ")")
        ) +
        coord_flip()
ggsave("figures/supplementary/ISM_transcripts_gene_type.pdf", width = 3, height = 2.5)

# ISM subcategory percentages

read_parquet("nextflow_results/V47/final_classification.parquet") %>% 
    filter(structural_category == "incomplete-splice_match") %>% 
    mutate(
        subcategory = factor(subcategory, levels=rev(c("5prime_fragment", "3prime_fragment", "intron_retention", "internal_fragment")))
    ) %>% 
    ggplot(aes(x = subcategory, fill=subcategory)) +
        geom_bar() +
        scale_fill_brewer(palette = "Set1") +
        scale_y_continuous(labels = function(x) x / 1000) +
        labs(
            x = "",
            y = expression("Transcripts (x" ~ 10^3 * ")")
        ) +
        coord_flip()
ggsave("figures/supplementary/ISM_transcript_subtype.pdf", width = 4.5, height = 2.5)