library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(VariantAnnotation)
library(tidyr)
library(parallel)

gencode_exons <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
    exonsBy(by = "tx", use.names = TRUE) %>% 
    unlist() %>% 
    unique()

SFARI_exons <- rtracklayer::import("nextflow_results/V47/orfanage/orfanage.gtf") %>%
    subset(type == "exon") %>% 
    unique()

novel_exonic_regions <- GenomicRanges::subtract(SFARI_exons, gencode_exons, ignore.strand=TRUE) %>% unlist() %>% unique() %>% reduce(ignore.strand=TRUE)

novel_exonic_regions %>% saveRDS("export/variant/novel_exonic_regions.rds")

# Is it possible for all the coordinates in a gene to be 1 bp apart?
# Use GENCODE to test this

SFARI_exons <- rtracklayer::import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf"))

sample(pull(distinct(as_tibble(mcols(SFARI_exons)), gene_id), gene_id), size=10) %>% 
    mclapply(get_dcoord_per_gene_id) %>% 
    bind_rows() %>% 
    filter(dcoord < 200)  %>% 
    ggplot(aes(dcoord)) +
    geom_bar()





get_dcoord_per_gene_id <- function(gene_id) {
    subset(SFARI_exons, (gene_id=={{gene_id}})&(type=="exon")) %>% 
        as_tibble() %>% 
        dplyr::select(start, end) %>% 
        pivot_longer(cols=c(start, end), names_to="start_end", values_to="coord") %>% 
        arrange(coord) %>% 
        distinct(coord) %>% 
        mutate(dcoord = coord - lag(coord))
}

data.frame(width = width(novel_exonic_regions)) %>%
    filter(width < 200)  %>% 
    ggplot(aes(width)) +
    geom_bar()