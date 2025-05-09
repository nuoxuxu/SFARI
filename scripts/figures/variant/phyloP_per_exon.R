library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(ggplot2)
library(patchwork)

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "none"
    )

theme_set(my_theme)

novel_exonic_regions <- readRDS("export/variant/novel_exonic_regions.rds")

gencode_exons <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
    exonsBy(by = "tx", use.names = TRUE) %>% 
    unlist() %>% 
    reduce()

get_mean_phylop <- function(exon_gr) {
    bw_file <- BigWigFile("data/hg38.phyloP100way.bw")
    phylop <- import(bw_file, which = exon_gr, as = "NumericList")
    mean_phylop <- sapply(phylop, function(scores) {
        ifelse(length(scores) == 0, NA, mean(scores, na.rm = TRUE))
    })
    mcols(exon_gr)$avg_phylop <- mean_phylop
    exon_gr
}    

novel_exonic_regions <- get_mean_phylop(novel_exonic_regions)
gencode_exons <- get_mean_phylop(gencode_exons)

novel_exonic_regions_g <- mcols(novel_exonic_regions) %>% 
    as_tibble() %>% 
        ggplot(aes(x = avg_phylop)) +
        geom_histogram() +
        xlim(-1, 6)

gencode_exons_g <- mcols(gencode_exons) %>% 
    as_tibble() %>% 
        ggplot(aes(x = avg_phylop)) +
        geom_histogram() +
        xlim(-1, 6)

novel_exonic_regions_g / gencode_exons_g
ggsave("figures/test/phyloP_per_exon.pdf", width = 7.5, height = 7)