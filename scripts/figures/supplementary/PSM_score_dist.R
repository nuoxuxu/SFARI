library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(rtracklayer)

percolator_res <- read_tsv("nextflow_results/V47/orfanage/hybrid_percolator.tsv") %>%
    # Replace proteinIds: take the first ID (splitting on comma)
    mutate(proteinIds = sapply(str_split(proteinIds, ","), `[`, 1)) %>%
    # Replace all occurrences of "M[15.9949]" with "M" in the peptide column
    mutate(peptide = str_replace_all(peptide, "M\\[15\\.9949\\]", "M")) %>%
    # Split the peptide string (by period) and create new columns:
    mutate(
        prev_aa = sapply(str_split(peptide, "\\."), `[`, 1),
        pep = sapply(str_split(peptide, "\\."), `[`, 2),
        next_aa = sapply(str_split(peptide, "\\."), `[`, 3)
    ) %>%
    # Keep only unique rows based on the 'pep' column
    distinct(pep, .keep_all = TRUE)

annot_peptides_hybrid <- import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>%
    as_tibble() %>%
    distinct(transcript_id, .keep_all = TRUE)

my_theme <- theme_bw() +
    theme(
        plot.title = element_text(size = 25),
        axis.text.x = element_text(size = 18),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 20)
    )

theme_set(my_theme)

annot_peptides_hybrid %>%
    left_join(
        percolator_res[, c("PSMId", "score")],
        join_by(gene_id == PSMId)
    ) %>%
    mutate(
        novelty = factor(novelty, levels = c("novel", "known"))
    ) %>% 
    ggplot(aes(x = score.y, fill = novelty)) +
    geom_density(alpha=0.4) +
    labs(y="", x="PSM score")
ggsave("figures/figure_2/PSM_score_dist.pdf")