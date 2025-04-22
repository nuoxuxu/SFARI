library(dplyr)
library(rtracklayer)
library(GenomicFeatures)
library(GenomicAlignments)

annotation_gtf <- paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf")

SFARI_exon <- makeTxDbFromGFF(predicted_cds_gtf) %>% 
    cdsBy(by = "tx", use.names=TRUE)

# Read in the final_transcript GTF file
gtf <- rtracklayer::import("nextflow_results/V47/final_transcripts.gtf")

# Filter for exons only
exons <- gtf[gtf$type == "exon"]

# Subset the exons to only unique
exons.unique <- unique(exons)    

GENCODE_exon <- makeTxDbFromGFF(annotation_gtf) %>%
    exonsBy(by = "tx", use.names=TRUE) %>% 
    GRangesList(compress=FALSE)
    
GRangesList(GENCODE_exon) %>% class()

queryHits(findSpliceOverlaps(exons.unique, GENCODE_exon))

pbid_containing_novel_exons <- names(SFARI_CDS[idx])



matches <- findOverlaps(exons.unique, GENCODE_exon, type = "equal")
novelty_labels <- rep("novel", length(exons.unique))
novelty_labels[unique(queryHits(matches))] <- "known"
mcols(exons.unique)$novelty <- novelty_labels

length(which(exons.unique$novelty == "known"))
length(which(exons.unique$novelty =="novel"))
