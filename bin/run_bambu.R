#!/usr/bin/env Rscript
library(bambu)
library(rtracklayer)
library(arrow)
library(readr)

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
    processByChromosome = TRUE,
)

saveRDS(bambu_result, output_rds)

do.call(cbind, lapply(bambu_result, function(df) assays(df)$counts)) %>% 
    as.matrix() %>% 
    write.csv("nextflow_results/comapre_other_LRS_tools/bambu/bambu_expression.csv")

bambu_result <- bambu_result[[1]]
supportedTx <- !is.na(rowData(bambu_result)$readCount)
writeToGTF(rowRanges(bambu_result)[supportedTx, ], file = "nextflow_results/comapre_other_LRS_tools/bambu/supportedTranscriptModels.gtf")
data.frame(rowData(bambu_result)[supportedTx, ]) %>% 
    write_tsv(file="nextflow_results/comapre_other_LRS_tools/bambu/supportedTxClassification.txt")