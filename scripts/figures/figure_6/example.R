library(Gviz)
library(txdbmaker)
library(rtracklayer)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggtranscript)
library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg38)

# Known track

GENCODE_GRList <- paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf") %>%
  makeTxDbFromGFF(format = "gtf") %>%
  exonsBy(by = "tx", use.names = TRUE)

gencode_gtf <- import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
  as_tibble()

gencode_tx <- gencode_gtf %>% 
  as_tibble() %>% 
  filter(gene_name=="ZNF292") %>% 
  filter(type=="transcript") %>% 
  filter(tag=="GENCODE_Primary") %>%  
  pull(transcript_id)

gencode_track <- GENCODE_GRList[gencode_tx, ] %>% 
  unlist() %>% 
  reduce() %>% 
  GeneRegionTrack(name = "Known")

ranges(gencode_track)$transcript <- rep("transcript", length(ranges(gencode_track)))

# Novel track

novel_pb_id <- "PB.44758.107"
SFARI_txdb <- "nextflow_results/V47/orfanage/orfanage.gtf" %>% 
  makeTxDbFromGFF(format = "gtf")
SFARI_track <- GeneRegionTrack(SFARI_txdb, group = "transcript", name = "Novel")
SFARI_track_subset <- SFARI_track[transcript(SFARI_track) == novel_pb_id]

# Plot gencode track

bm <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome = "hg38", name = "ENSEMBL", transcript = c("ENST00000699925"), biomart = bm)
grtrack_subset <- biomTrack[transcript(biomTrack) == "ENST00000699925"]

vline <- AnnotationTrack(
  start = 87256249,
  end   = 87256250,   # same start and end → vertical line
  chromosome = "chr6",
)

# Plot sequence track

sTrack <- SequenceTrack(Hsapiens)

# Combine all tracks

chr <- as.character(seqnames(SFARI_track_subset))[1]
axisTrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = chr)

pdf("figures/figure_6/sequence.pdf", width = 10, height = 2)
plotTracks(
  list(
    sTrack,
    vline),
  chromosome = chr,
  from = 87256214,
  to = 87256280
)
dev.off()

# ggtranscript

my_theme <- theme_bw() +
    theme(
        axis.text.x = element_text(size = 20, angle = 45, color = "black", vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 20, color = "black"),
        strip.placement = "outside",
        legend.title = element_blank(),
        legend.position = c(0.08, 0.90),
        legend.text = element_text(size = 20)
    )

theme_set(my_theme)

ZNF292_known_annotation <- import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
  as_tibble() %>% 
  filter(transcript_id == "ENST00000369577.8")

ZNF292_novel_annotation <- import("nextflow_results/V47/orfanage/orfanage.gtf") %>% 
  as_tibble() %>% 
  filter(transcript_id == novel_pb_id)

# Combine

ZNF292_known_annotation <- ZNF292_known_annotation %>% select(c(seqnames, start, end, strand, transcript_id, type))
ZNF292_novel_annotation <- ZNF292_novel_annotation %>% select(c(seqnames, start, end, strand, transcript_id, type))

combined_ZNF292_annotation <- bind_rows(ZNF292_known_annotation, ZNF292_novel_annotation) %>% 
  mutate(
    transcript_id = factor(transcript_id, levels = c(novel_pb_id, "ENST00000369577.8"))
  )

ZNF292_exons <- combined_ZNF292_annotation %>%
  filter(type == "exon")

ZNF292_CDS <- combined_ZNF292_annotation %>%
  filter(type == "CDS")

ZNF292_exons %>% 
  ggplot(
    aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )
  ) +
  geom_range(
    fill = "white",
    height = 0.25
  ) +
  geom_intron(
      data = to_intron(ZNF292_exons, "transcript_id"),
      aes(strand = strand)
  ) + 
  geom_range(
      data = ZNF292_CDS
  ) +
  geom_intron(
      data = to_intron(ZNF292_CDS, "transcript_id"),
      aes(strand = strand),
      arrow.min.intron.length = 500,
  ) +
  geom_vline(xintercept = 87256249, color = "red", linetype = "dashed")

ggsave("figures/figure_6/ZNF292_full_transcript.pdf", width = 9, height = 5)
