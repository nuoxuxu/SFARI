#!/usr/bin/env Rscript
library(Gviz)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

options(ucscChromosomeNames=FALSE)

args <- commandArgs(trailingOnly=TRUE)
bam_file    <- args[1]
gtf_file    <- args[2]
gene_name   <- args[3]
gencode_gtf <- args[4]
out_pdf     <- args[5]

gencode <- import(gencode_gtf) %>%
    as_tibble() %>%
    filter(type == "gene")

txdb <- makeTxDbFromGFF(gtf_file)

chr    <- gencode %>% filter(gene_name == !!gene_name) %>% pull(seqnames) %>% as.vector()
aStart <- gencode %>% filter(gene_name == !!gene_name) %>% pull(start) %>% as.vector()
aEnd   <- gencode %>% filter(gene_name == !!gene_name) %>% pull(end) %>% as.vector()

alTrack <- AlignmentsTrack(bam_file)
txTr    <- GeneRegionTrack(txdb, chromosome = chr, start = aStart - 2000, end = aEnd)

pdf(out_pdf, width = 10, height = 6)
plotTracks(
    list(txTr, alTrack),
    chromosome = chr, from = aStart - 2000,
    to = aEnd,
    transcriptAnnotation = "transcript"
)
dev.off()
