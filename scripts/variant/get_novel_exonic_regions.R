library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(tidyr)

gencode_exons <- makeTxDbFromGFF(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
    exonsBy(by = "tx", use.names = TRUE) %>% 
    unlist() %>% 
    unique()

SFARI_exons <- rtracklayer::import("nextflow_results/V47/orfanage/orfanage.gtf") %>%
    subset(type == "exon") %>% 
    unique()

novel_exonic_regions <- GenomicRanges::subtract(SFARI_exons, gencode_exons, ignore.strand=TRUE) %>% unlist() %>% unique() %>% reduce(ignore.strand=TRUE)

novel_exonic_regions %>% saveRDS("export/variant/novel_exonic_regions.rds")