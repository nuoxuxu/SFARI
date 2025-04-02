library(arrow)
library(readr)
library(ggplot2)
library(rtracklayer)
library(ggpubr)
library(tidyr)
library(dplyr)

#---------------------------Read in datasets---------------------------------#

protein_class <- read_tsv("nextflow_results/V47/orfanage/SFARI.protein_classification.tsv")

expression <- read_parquet(("nextflow_results/V47/final_expression.parquet"))

annot_peptides_hybrid <- rtracklayer::import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>% 
    as_tibble() %>% 
    distinct(transcript_id, .keep_all = TRUE)

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet") %>% 
    left_join(
        annot_peptides_hybrid,
        join_by(peptide == transcript_id)
    ) %>% 
    dplyr::filter(detected=="True")

mean_expression <- expression %>% 
    left_join(
        protein_class[, c("pb", "base_isoform")],
        by = c("isoform" = "pb")
    ) %>% 
    filter(!is.na(base_isoform)) %>% 
    mutate(across(-c(isoform, base_isoform), ~ log2((.x / sum(.x) * 1e6) + 1))) %>% 
    group_by(base_isoform) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))%>%
    mutate(row_mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE))

# sum_expression <- expression %>% 
#     left_join(
#         protein_class[, c("pb", "base_isoform")],
#         by = c("isoform" = "pb")
#     ) %>% 
#     filter(!is.na(base_isoform)) %>% 
#     mutate(across(-c(isoform, base_isoform), ~ log2((.x / sum(.x) * 1e6) + 1))) %>% 
#     group_by(base_isoform) %>%
#     summarise(across(where(is.numeric), sum, na.rm = TRUE))%>%
#     mutate(row_mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE))

sum_expression <- expression %>% 
    left_join(
        protein_class[, c("pb", "base_isoform")],
        by = c("isoform" = "pb")
    ) %>% 
    filter(!is.na(base_isoform)) %>% 
    group_by(base_isoform) %>%
    summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 
    mutate(across(-c(base_isoform), ~ log2((.x / sum(.x) * 1e6) + 1))) %>% 
    mutate(row_mean = rowMeans(dplyr::select(., where(is.numeric)), na.rm = TRUE))

sum_expression %>% write_parquet("export/sum_expression.parquet")

#---------------------------Set themes---------------------------------#

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

#---------------------------Start Here---------------------------------#

# Read in datasets

sum_expression <- read_parquet("export/sum_expression.parquet")

annot_peptides_hybrid <- rtracklayer::import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>% 
    as_tibble() %>% 
    distinct(transcript_id, .keep_all = TRUE)

#  At least 1 peptides

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet") %>% 
    left_join(
        annot_peptides_hybrid,
        join_by(peptide == transcript_id)
    ) %>%  
    dplyr::filter(detected=="True") %>% 
    group_by(pb) %>%
    summarise(n_peptides = n()) %>% 
    filter(n_peptides >= 1)

p1 <- sum_expression %>% 
    mutate(
        validated = case_when(
            base_isoform %in% peptide_mapping$pb ~ TRUE,
            .default = FALSE
        )
    ) %>% 
    ggplot(
        aes(validated, row_mean, fill=validated)
        ) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "Protein isoforms", y = "Expression\n(log2(CPM + 1))") +
        ylim(0, 2) +
        scale_x_discrete(labels = c("FALSE" = "Not validated", "TRUE" = "Validated"))

# At least 5 peptides

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet") %>% 
    left_join(
        annot_peptides_hybrid,
        join_by(peptide == transcript_id)
    ) %>%  
    dplyr::filter(detected=="True") %>% 
    group_by(pb) %>%
    summarise(n_peptides = n()) %>% 
    filter(n_peptides >= 5)

p2 <- sum_expression %>% 
    mutate(
        validated = case_when(
            base_isoform %in% peptide_mapping$pb ~ TRUE,
            .default = FALSE
        )
    ) %>% 
    ggplot(
        aes(validated, row_mean, fill=validated)
        ) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "Protein isoforms", y = "") +
        ylim(0, 2) +
        scale_x_discrete(labels = c("FALSE" = "Not validated", "TRUE" = "Validated"))

# At least 10 peptides

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet") %>% 
    left_join(
        annot_peptides_hybrid,
        join_by(peptide == transcript_id)
    ) %>%  
    dplyr::filter(detected=="True") %>% 
    group_by(pb) %>%
    summarise(n_peptides = n()) %>% 
    filter(n_peptides >= 10)

p3 <- sum_expression %>% 
    mutate(
        validated = case_when(
            base_isoform %in% peptide_mapping$pb ~ TRUE,
            .default = FALSE
        )
    ) %>% 
    ggplot(
        aes(validated, row_mean, fill=validated)
        ) +
        geom_boxplot(outlier.shape = NA) +
        labs(x = "Protein isoforms", y = "") +
        ylim(0, 2) +
        scale_x_discrete(labels = c("FALSE" = "Not validated", "TRUE" = "Validated"))

#---------------------------Export---------------------------------#        
library(patchwork)
p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave("figures/figure_2/protein_isoform_abundance_by_whether_validated_by_n_peptides.pdf", width = 9, height = 2.5)

df <- sum_expression %>% 
    mutate(
        validated = case_when(
            base_isoform %in% peptide_mapping$pb ~ TRUE,
            .default = FALSE
        )
    )

t.test(df %>% filter(validated) %>% pull(row_mean), df %>% filter(!validated) %>% pull(row_mean))

# "Validated" means that the protein isoform is validated by at least 5 peptide

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet") %>% 
    left_join(
        annot_peptides_hybrid,
        join_by(peptide == transcript_id)
    ) %>% 
    dplyr::filter(detected=="True") %>% 
    group_by(pb) %>%
    summarise(n_peptides = n())

sum_expression %>% 
    mutate(
        validated = case_when(
            base_isoform %in% peptide_mapping$pb ~ TRUE,
            .default = FALSE
        )
    ) %>% 
    ggplot(aes(validated, row_mean)) +
        geom_boxplot() +
        labs(x = "Protein isoforms validated by peptides", y = "Mean expression across all samples all time points\n(log2(CPM + 1))") +
        ggtitle('"Validated" means that the protein isoform is validated by at least 5 peptide')        

# novel peptides only

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet") %>% 
    left_join(
        annot_peptides_hybrid,
        join_by(peptide == transcript_id)
    ) %>% 
    dplyr::filter(detected=="True") %>% 
    dplyr::filter(!GENCODE)

mean_expression %>% 
    left_join(
        protein_class[, c("base_isoform", "protein_classification_base")] %>% distinct(base_isoform, .keep_all = TRUE),
        join_by(base_isoform == base_isoform)
    ) %>% 
    filter(protein_classification_base %in% c("pNNC")) %>% 
    mutate(
        validated = case_when(
            base_isoform %in% peptide_mapping$pb ~ TRUE,
            .default = FALSE
        )
    ) %>% 
    ggplot(aes(validated, row_mean)) +
        geom_boxplot() +
        labs(x = "Protein isoforms validated by peptides", y = "Mean expression across all samples all time points\n(log2(CPM + 1))")
ggsave("figures/figure_2/protein_isoform_abundance_by_whether_validated.pdf", width = 5, height = 5)
# histogram as in Figure 5 BCD in Sinitcyn et al.

peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet") %>% 
    left_join(
        annot_peptides_hybrid,
        join_by(peptide == transcript_id)
    ) %>% 
    dplyr::filter(detected=="True")

mean_expression %>% 
    left_join(
        peptide_mapping
    )