library(arrow)
library(ggplot2)
library(dplyr)

# FSM transcript gene type percentages

read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    filter(
        structural_category == "full-splice_match"
    ) %>%
    distinct(associated_transcript) %>%
    left_join(
        GENCODE_gtf[, c("transcript_id", "gene_type")],
        join_by(associated_transcript == transcript_id)
    ) %>%
    distinct(associated_transcript, .keep_all=TRUE) %>%
    mutate(gene_type = if_else(gene_type %in% c("processed_pseudogene", "unprocessed_pseudogene"), "pseudogene", gene_type)) %>%
    filter(
        gene_type %in% c("lncRNA", "protein_coding", "pseudogene")
    ) %>%
    mutate(gene_type = factor(gene_type, levels = c("protein_coding", "lncRNA", "pseudogene"))) %>%
    ggplot(aes(x = gene_type)) +
        geom_bar(fill = "#009E73") +
        scale_y_continuous(labels = function(x) x / 1000) +
        labs(
            x = "",
            y = expression("Transcripts (x" ~ 10^3 * ")")
        ) +        
        coord_flip()
ggsave("figures/supplementary/FSM_transcripts_gene_type.pdf", width = 3, height = 2.5)

# FSM transcript subcategory percentages

read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    filter(
        structural_category == "full-splice_match"
    ) %>%
    mutate(
        subcategory = factor(subcategory, levels=rev(c("reference_match", "alternative_3end", "alternative_5end", "alternative_3end5end", "mono-exon")))
    ) %>% 
    ggplot(aes(x = subcategory, fill=subcategory)) +
        geom_bar() +
        scale_fill_manual("Structural category", values=colorVector) +
        scale_y_continuous(labels = function(x) x / 1000) +
        labs(
            x = "",
            y = expression("Transcripts (x" ~ 10^3 * ")")
        ) +
        coord_flip()
ggsave("figures/supplementary/FSM_transcript_subtype.pdf", width = 4.5, height = 2.5)