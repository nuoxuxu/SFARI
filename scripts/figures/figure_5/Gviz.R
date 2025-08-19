library(Gviz)
library(biomaRt)
library(rtracklayer)
library(GenomicFeatures)
library(readr)
library(dplyr)
library(arrow)

mart <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl")
df <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = "SMARCA4",
  mart       = mart
)
gr <- GRanges(df$chromosome_name, IRanges(df$start_position, df$end_position),
  strand = ifelse(df$strand == 1, "+", "-")
)

txdb <- makeTxDbFromGFF("/project/s/shreejoy/Genomic_references/GENCODE/gencode.v47.annotation.gtf")
txTr <- GeneRegionTrack(txdb, chromosome = chr, start = from, end = to)

chr <- "chr19"
from <- 10960932
to <- 11079426

iPSC <- DataTrack(
  range      = "data/katherine/iPSC_merged.bw", # Gviz will import just the requested window
  genome     = "hg38",
  chromosome = chr,
  name       = "iPSC",
  type       = "h" # "h" = histogram; try "l" for a line
)

NPC <- DataTrack(
  range      = "data/katherine/NPC_merged.bw", # Gviz will import just the requested window
  genome     = "hg38",
  chromosome = chr,
  name       = "NPC",
  type       = "h" # "h" = histogram; try "l" for a line
)

CN <- DataTrack(
  range      = "data/katherine/CN_merged.bw", # Gviz will import just the requested window
  genome     = "hg38",
  chromosome = chr,
  name       = "CN",
  type       = "h" # "h" = histogram; try "l" for a line
)

plotTracks(
  list(GenomeAxisTrack(), iPSC, NPC, CN, txTr),
  chromosome = chr, from = from, to = to
)