#!/usr/bin/env Rscript
library(bambu)
library(dplyr)
library(rtracklayer)
library(arrow)
library(readr)

args = commandArgs(trailingOnly=TRUE)
input_rds <- args[1]

bambu_result <- readRDS(input_rds)

do.call(cbind, lapply(bambu_result, function(df) assays(df)$counts)) %>%
    as.matrix() %>%
    write.csv("bambu_expression.csv")

bambu_result <- bambu_result[[1]]
supportedTx <- !is.na(rowData(bambu_result)$readCount)
writeToGTF(rowRanges(bambu_result)[supportedTx, ], file = "supportedTranscriptModels.gtf")
data.frame(rowData(bambu_result)[supportedTx, ]) %>%
    write_tsv(file = "supportedTxClassification.txt")
