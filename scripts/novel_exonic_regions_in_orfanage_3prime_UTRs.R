library(GenomicRanges)
library(dplyr)
library(readr)
library(rtracklayer)

# Jimmy is asking if we have a gtf with novel regions that are in the 3' UTRs 
# (this will act as his control for his Ribo-seq analyses). 
# I think it would make sense to use orfanage 3'UTRs here since he's working with the orfanage ORF?

novel_exonic_regions <- import("nextflow_results/V47/orfanage/UCSC_tracks/novel_exonic_regions.gtf") %>% 
    subset(width(.) >= 10)
    
orfanage_gtf <- import("nextflow_results/V47/orfanage/orfanage.gtf")

orfanage_three_prime_utr <- subset(orfanage_gtf, mcols(orfanage_gtf)$original_biotype=="three_prime_UTR")

findOverlaps(novel_exonic_regions, orfanage_three_prime_utr, type="within") %>% 
    queryHits() %>% 
    unique() %>% 
    novel_exonic_regions[.] %>% 
    export("export/novel_exonic_regions_in_orfanage_3prime_UTRs.gtf")
