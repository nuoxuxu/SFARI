#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
import polars as pl
from src.single_cell import SingleCell
from src.utils import read_gtf
from src.ryp import r, to_r
import argparse

r(
"""
library(IsoformSwitchAnalyzeR)
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(scales)
library(dplyr)
"""
)

def main():
    parser = argparse.ArgumentParser(description="Prepare files for running IsoformSwitchAnalyzeR")
    parser.add_argument("--h5ad_file", action="store", type=str, required=True)
    parser.add_argument("--filtered_gff", action="store", type=str, required=True)
    parser.add_argument("--fasta_file", action="store", type=str, required=True)
    parser.add_argument("--rds_file", action="store", type=str, required=True)

    params = parser.parse_args()
    
    lr_bulk = SingleCell(params.h5ad_file)

    Isoseq_Expression = lr_bulk.to_frame().rename({"isoform": "isoform_id"})
    to_r(Isoseq_Expression, "Isoseq_Expression")

    read_gtf(params.filtered_gff)\
        .filter(pl.col("transcript_id").is_in(Isoseq_Expression["isoform_id"]) )\
        .drop("transcript_id")\
        .write_csv("isoformExonAnnoation.gtf", separator="\t", include_header=False, quote_style="never")

    r(
f"""
sampleID <- colnames(Isoseq_Expression)[c(-1)]
time_point <- str_split(sampleID, "_", 2) %>% map_chr(~ .x[1])
time_point <- factor(time_point, levels = c("iPSC", "NPC", "CN"))
myDesign <- data.frame(sampleID = sampleID, condition = time_point)

IsoseqsSwitchList <- importRdata(
    isoformCountMatrix = Isoseq_Expression,
    designMatrix = myDesign,
    isoformExonAnnoation = "isoformExonAnnoation.gtf",
    isoformNtFasta = "{params.fasta_file}",
    addAnnotatedORFs = FALSE,
    fixStringTieAnnotationProblem = FALSE
)
saveRDS(IsoseqsSwitchList, "{params.rds_file}")
"""
    )

if __name__ == "__main__":
    main()