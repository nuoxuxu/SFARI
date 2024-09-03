library(tidyverse)
library(ggvenn)
library(ggplot2)

illumina <- read_csv("results/short_read/combined_batches_DE_results.csv")
pacbio <- read_csv("results/long_read/combined_batches_DE_results.csv")

contrast <- "NPC_iPSC.adj.P.Val"

plot_ggvenn <- function(contrast) {
    illumina_genes <- illumina[illumina[, contrast] < 0.05, ] %>% pull(Gene_ID)
    pacbio_genes <- pacbio[pacbio[, contrast] < 0.05, ] %>% pull(Gene_ID)
    ggvenn(list(illumina_DEGs = illumina_genes, pacbio_DEGs = pacbio_genes), fill_color = c("lightblue", "lightgreen"))
}

plot_ggvenn("Neuron_iPSC.adj.P.Val")
ggsave("figures/venn_Neuron_iPSC.png")
plot_ggvenn("NPC_iPSC.adj.P.Val")
ggsave("figures/venn_NPC_iPSC.png")
plot_ggvenn("Neuron_NPC.adj.P.Val")
ggsave("figures/venn_Neuron_NPC.png")