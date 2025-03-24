library(Gviz)
library(GenomicFeatures)
library(GenomicRanges)
library(GenomicInteractions)
library(glue)
library(rtracklayer)
library(dplyr)
library(data.table)

options(stringsAsFactors = F)
options(Gviz.scheme = "myScheme")
options(ucscChromosomeNames = F)

args <- commandArgs(trailingOnly = TRUE)

scheme <- getScheme()
scheme$GeneRegionTrack$col <- NULL
addScheme(scheme, "myScheme")

gencode= paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")
gencode_txdb=makeTxDbFromGFF(gencode, format="gtf")
gencode_transcript=exonsBy(gencode_txdb,by="tx",use.names=T)

isoseq="nextflow_results/V47/orfanage/orfanage.gtf"
isoseq_txdb=makeTxDbFromGFF(isoseq, format="gtf")
isoseq_transcript=exonsBy(isoseq_txdb,by="tx", use.names=T)

peptide_txdb <- makeTxDbFromGFF("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf", format="gtf")
peptide_transcript=exonsBy(peptide_txdb,by="tx",use.names=T)
peptide_annot <- rtracklayer::import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>% as_tibble()

current_gene=args[1]
genes_to_plot=read.delim("export/plotgene_transcripts_codingOnly.txt",header = F)
chr=unique(genes_to_plot$V1[genes_to_plot$V2 == current_gene])
strd=unique(genes_to_plot$V4[genes_to_plot$V2 == current_gene])

options(ucscChromosomeNames=FALSE)
gencode_only=base::setdiff(with(genes_to_plot,V3[V2 == current_gene & V5=="gencode"]) , with(genes_to_plot,V3[V2 == current_gene & V5=="isoseq"]))
shared=base::intersect(with(genes_to_plot,V3[V2 == current_gene & V5=="gencode"]) , with(genes_to_plot,V3[V2 == current_gene & V5=="isoseq"]))

peptides <- peptide_annot %>% 
  dplyr::filter(
    gene_name == current_gene,
    detected == "True",
    novelty == "novel"
    ) %>%
  pull(transcript_id) %>% 
  unique()

gencode_only_transcript_onegene=gencode_transcript[gencode_only,]
gencode_only_transcript_onegene=unlist(gencode_only_transcript_onegene)
gencode_only_transcript_onegene <- reduce(gencode_only_transcript_onegene)
# elementMetadata(gencode_only_transcript_onegene)$transcript <- "collapsed"
# gencode_only_track=GeneRegionTrack(gencode_only_transcript_onegene,group = "transcirpt",name = "Gencode undetected")
gencode_only_track <- GeneRegionTrack(gencode_only_transcript_onegene, name="Known")
# ranges(geneTrack)$feature <- as.character(ranges(geneTrack)$feature)
ranges(gencode_only_track)$transcript <- rep("transcript", length(ranges(gencode_only_track)))

novel <- args[2]
isoseq_transcript_onegene=isoseq_transcript[novel,]
isoseq_transcript_onegene=unlist(isoseq_transcript_onegene)
elementMetadata(isoseq_transcript_onegene)$transcript=names(isoseq_transcript_onegene)
isoseq_track=GeneRegionTrack(isoseq_transcript_onegene,group = "transcript",name = "Novel")

peptide_onegene <- peptide_transcript[peptides,]
peptide_onegene <- unlist(peptide_onegene)
elementMetadata(peptide_onegene)$transcript=names(peptide_onegene)
peptide_track=GeneRegionTrack(peptide_onegene,group = "transcript",name = "Novel\nPeptides")

displayPars(gencode_only_track)=list(
  stacking="squish",
  background.panel="transparent",
  fill="#009E73",
  col="#009E73",
  lwd=0.3,
  col.line="black",
  fontcolor.title="black",
  background.title="#d1861d")

displayPars(isoseq_track)=list(
  stacking="squish",
  fill="#E69F00",
  col="#E69F00",
  lwd=0.3,
  col.line="black",
  showId = TRUE,
  background.panel="transparent",
  transcriptAnnotation = "transcript",
  fontcolor.title="black",
  background.title="#045a8d")
displayPars(peptide_track)=list(
  stacking="squish",
  background.panel="transparent",
  fill="black",
  col="black",
  lwd=0.3,
  col.line="black",
  showId = TRUE,
  transcriptAnnotation = "transcript",
  fontcolor.title="black",
  background.title="#ef6548"
  )                               

axisTrack <- GenomeAxisTrack()
ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = chr)

leftmost=min(c(gencode_only_transcript_onegene@ranges@start,isoseq_transcript_onegene@ranges@start))
rightmost=max(c(gencode_only_transcript_onegene@ranges@start,isoseq_transcript_onegene@ranges@start))

extra=(rightmost - leftmost)*0.05


pdf(glue("figures/figure_2/genome_track_{current_gene}.pdf"), width = 9, height = 3)
plotTracks(
  list(ideoTrack,axisTrack,gencode_only_track,isoseq_track, peptide_track),
  chromosome = chr,
  from = leftmost - extra,
  to = rightmost + extra,
  sizes = c(
    1,1,1,1,1
    )
  )
dev.off()