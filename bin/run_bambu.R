#!/usr/bin/env Rscript
library(bambu)

args = commandArgs(trailingOnly=TRUE)
annotation_gtf <- args[1]
aligned_bam <- args[2]
ref_genome_fasta <- args[3]
ncore <- args[4]
output_rds <- args[5]

bambuAnnotations <- prepareAnnotations(annotation_gtf)
bam_files <- Sys.glob(aligned_bam)

bambu_result <- bambu::bambu(
    reads = bam_files, 
    annotations = bambuAnnotations,
    genome = ref_genome_fasta,
    ncore = ncore, 
    opt.discovery = list(remove.subsetTx = FALSE, min.readCount=1),
    quant = FALSE, 
    discovery = TRUE,
    lowMemory = TRUE,
)

saveRDS(bambu_result, output_rds)