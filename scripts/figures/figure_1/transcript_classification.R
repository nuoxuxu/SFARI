library(dplyr)
library(ggplot2)
library(arrow)
library(scales)
library(patchwork)

colorVector <- c(
    "FSM" = "#009E73",
    "ISM" = "#0072B2",
    "NIC" = "#D55E00",
    "NNC" = "#E69F00",
    "Other" = "#000000"
)

structural_category_labels <- c(
    "full-splice_match"        = "FSM",
    "incomplete-splice_match"  = "ISM",
    "novel_in_catalog"         = "NIC",
    "novel_not_in_catalog"     = "NNC",
    "Other"                    = "Other"
)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

# Transcript Classification

transcript_classification_hist <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    mutate(
        # If category is in our known set, keep it; otherwise use "Other"
        structural_category2 = if_else(
            structural_category %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            structural_category,
            "Other"
        ),
        # Map to short labels (FSM, ISM, NIC, NNC, Other)
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%
    group_by(structural_category2) %>%
    summarise(len = n()) %>%
    mutate(
        percentage = (len / sum(len)) * 100,
        structural_category2 = factor(structural_category2, levels = c("FSM", "ISM", "NIC", "NNC", "Other"))
    ) %>%
    ggplot(aes(x = structural_category2, y = len, fill = structural_category2)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(percentage, 1), "%")), vjust = 2, colour = "white", size = 2.5) +
        scale_y_continuous(labels = function(x) x / 1000) +
        scale_fill_manual("Structural Category", values = colorVector) +
        labs(x = "Structural Category", y = expression("Transcripts (x" ~ 10^3 * ")")) +
        xlab(NULL)
ggsave("figures/figure_1/transcript_classification_hist.pdf", width = 8, height = 7)

# Length vs abundance

read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    mutate(
        structural_category2 = if_else(
            structural_category %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            structural_category,
            "Other"
        ),
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%
    ggplot(aes(x = length, fill = structural_category2)) +
    geom_histogram(alpha = .75) +
    scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)), limits = c(300, 10^4)) +
    scale_y_continuous(labels = function(x) x / 1000) +
    scale_fill_manual("Structural Category", values = colorVector) +
    labs(
        x = "Transcript Length (bp)",
        y = expression("Transcripts (x" ~ 10^3 * ")")
    ) +
    ggtitle("Transcript length distribution")
ggsave("figures/figure_1/length_vs_abundance.pdf", width = 8, height = 7)

# Number of exons vs abudance

nexon_vs_abundance <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    mutate(
        structural_category2 = if_else(
            structural_category %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            structural_category,
            "Other"
        ),
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%
    ggplot(aes(x = exons, fill = structural_category2)) +
    geom_histogram(alpha = .75, binwidth = 1) +
    scale_x_continuous(limits = c(1, 40)) +
    scale_y_continuous(labels = function(x) x / 1000) +
    scale_fill_manual("Structural Category", values = colorVector) +
    labs(
        x = "# Exons",
        y = expression("Transcripts (x" ~ 10^3 * ")")
    )
ggsave("figures/figure_1/nexon_vs_abundance.pdf", width = 8, height = 7)

# Counts vs abundance

row_sum <- read_parquet("nextflow_results/V47/final_expression.parquet") %>%
    mutate(
        row_sum = rowSums(across(where(is.numeric)))
    ) %>%
    dplyr::select(row_sum)

count_vs_abundance <- read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    mutate(
        structural_category2 = if_else(
            structural_category %in% c(
                "full-splice_match",
                "incomplete-splice_match",
                "novel_in_catalog",
                "novel_not_in_catalog"
            ),
            structural_category,
            "Other"
        ),
        structural_category2 = structural_category_labels[structural_category2]
    ) %>%
    bind_cols(row_sum) %>%
    filter(row_sum > 10) %>%
    ggplot(aes(x = row_sum, fill = structural_category2)) +
    geom_histogram(position = position_fill(), alpha = .75) +
    scale_x_log10() +
    scale_fill_manual("Structural Category", values = colorVector) +
    labs(x = "Total observed counts", y = "Transcript proportion")
ggsave("figures/figure_1/count_vs_abundance.pdf", width = 8, height = 7)

# Transcript classification for the transcripts in the Other category

read_parquet("nextflow_results/V47/final_classification.parquet") %>%
    filter(structural_category %in% c("fusion", "genic", "antisense", "intergenic", "moreJunctions")) %>% 
    group_by(structural_category) %>%
    summarise(len = n()) %>%
    mutate(
        percentage = (len / sum(len)) * 100,
        structural_category = factor(structural_category)
    ) %>%
    ggplot(aes(x = structural_category, y = len, fill = structural_category)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(percentage, 1), "%")), colour = "black", size = 7, vjust=-0.5) +
        scale_y_continuous(labels = function(x) x / 1000, limits=c(0, 5000)) +
        labs(x = NULL, y = expression("Transcripts (x" ~ 10^3 * ")"))
ggsave("figures/figure_1/other_transcript_classification_hist.pdf", width = 10, height = 7)

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
        labs(x="", y="# FSM transcripts") +
        coord_flip()
ggsave("figures/figure_1/FSM_transcripts_gene_type.pdf", width = 8, height = 7)

# FSM transcript subcategory percentages

GENCODE_gtf <- rtracklayer::import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf")) %>%
    as_tibble()

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
        labs(x="", y="# FSM transcripts") +
        coord_flip()
ggsave("figures/figure_1/FSM_transcripts_gene_type.pdf", width = 8, height = 7)

transcript_classification_hist + count_vs_abundance + nexon_vs_abundance + plot_layout(widths = c(1.5, 1, 1))
ggsave("figures/figure_1/transcript_fig_combined.pdf", width = 7, height = 2.5)
