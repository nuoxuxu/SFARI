library(Gviz)
library(GenomicFeatures)
library(GenomicRanges)
library(glue)
library(tidyr)
library(rtracklayer)
library(dplyr)
library(arrow)

# Settings

options(stringsAsFactors = FALSE)
options(Gviz.scheme = "myScheme")
options(ucscChromosomeNames = FALSE)

scheme <- getScheme()
scheme$GeneRegionTrack$col <- NULL
addScheme(scheme, "myScheme")

# Current gene

# current_gene <- "EXOC1"
# novel_pb_id <- "PB.28078.226"
# tx_to_plot <- "ENST00000381295.7"

current_gene <- "ZMYND8"
novel_pb_id <- "PB.104062.23"
tx_to_plot <- "ENST00000360911.7"

# GENCODE track

GENCODE_GRList <- paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf") %>%
  makeTxDbFromGFF(format = "gtf") %>%
  exonsBy(by = "tx", use.names = TRUE)

gencode_gtf <- rtracklayer::import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
  as_tibble()

# tx_to_plot <- gencode_gtf %>%
#   filter( gene_name == current_gene, tag == "CCDS") %>%
#   distinct(transcript_id, .keep_all=TRUE) %>%
#   pull(transcript_id)

gencode_track <- GENCODE_GRList[tx_to_plot, ] %>%
  unlist() %>%
  reduce() %>%
  GeneRegionTrack(name = "Known")

ranges(gencode_track)$transcript <- rep("transcript", length(ranges(gencode_track)))

displayPars(gencode_track) <- list(
  stacking = "squish",
  background.panel = "transparent",
  fill = "#009E73",
  col = "#009E73",
  lwd = 0.3,
  col.line = "black",
  fontcolor.title = "black",
  background.title = "#d1861d"
)

# Peptide track

peptide_annot <- rtracklayer::import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>% as_tibble()

peptides <- peptide_annot %>%
  dplyr::filter(
    gene_.me == current_gene,
    detected == "True",
    novelty == "novel"
  ) %>%
  pull(transcript_id) %>%
  unique()

peptide_gr <- makeTxDbFromGFF("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf", format = "gtf") %>%
  exonsBy(by = "tx", use.names = TRUE) %>%
  .[peptides] %>%
  unlist()

elementMetadata(peptide_gr)$transcript <- names(peptide_gr)
peptide_track <- GeneRegionTrack(peptide_gr, group = "transcript", name = "Novel\nPeptides")

displayPars(peptide_track) <- list(
  stacking = "squish",
  background.panel = "transparent",
  fill = "black",
  col = "black",
  lwd = 0.3,
  col.line = "black",
  showId = TRUE,
  transcriptAnnotation = "transcript",
  fontcolor.title = "black",
  background.title = "#ef6548"
)

# Isoseq track

isoseq_gr <- "nextflow_results/V47/orfanage/orfanage.gtf" %>%
  makeTxDbFromGFF(format = "gtf") %>%
  exonsBy(by = "tx", use.names = TRUE) %>%
  .[novel_pb_id] %>%
  unlist()

elementMetadata(isoseq_gr)$transcript <- names(isoseq_gr)
isoseq_track <- GeneRegionTrack(isoseq_gr, group = "transcript", name = "Novel")

displayPars(isoseq_track) <- list(
  stacking = "squish",
  fill = "#E69F00",
  col = "#E69F00",
  lwd = 0.3,
  col.line = "black",
  showId = TRUE,
  background.panel = "transparent",
  transcriptAnnotation = "transcript",
  fontcolor.title = "black",
  background.title = "#045a8d"
)

# Combine all tracks

chr <- as.character(seqnames(isoseq_gr))[1]

leftmost <- isoseq_gr@ranges@start %>% min()
rightmost <- end(isoseq_gr)  %>% max()
extra <- (rightmost - leftmost) * 0.05

pdf(glue("figures/figure_2/genome_track_{current_gene}.pdf"), width = 9, height = 3)

plotTracks(
  list(
    gencode_track, isoseq_track, peptide_track
  ),
  chromosome = chr,
  from = leftmost - extra,
  to = rightmost + extra,
  sizes = c(
    1, 1, 1
  )
)
dev.off()