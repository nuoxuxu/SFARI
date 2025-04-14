library(dplyr)
library(ggplot2)
library(patchwork)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(arrow)
library(GenomicRanges)

LR_SJ <- read_parquet("data/riboseq/riboseq_SJ.parquet")

#-------------------------------------Get annotated long-read exons----------------------------

ribo <- BamFile("/scratch/s/shreejoy/ehogan/sfari/riboseq/analysis/peak-calling/merged.sorted.bam") %>% 
    scanBam(param = ScanBamParam(what = c("pos", "rname", "strand", "qwidth")))

ribo <- GRanges(
    seqnames = Rle(ribo[[1]]$rname), 
    ranges = IRanges(start = ribo[[1]]$pos, 
    width = ribo[[1]]$qwidth), 
    strand = ribo[[1]]$strand
    ) %>% 
    unique()

saveRDS(ribo, "data/riboseq/ribo.rds")

#--------------------------------------------END----------------------------------------------

#-------------------------------------Annotate long-read exons----------------------------

# Load in pre-processed data

ribo <- readRDS("data/riboseq/ribo.rds")

gencode <- rtracklayer::import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf"))
gencode.cds <- gencode[gencode$type == "CDS"]    
gencode.cds.unique <- unique(gencode.cds)

orfanage  <- rtracklayer::import("nextflow_results/V47/orfanage/orfanage.gtf")
orfanage.cds <- orfanage[orfanage$type == "CDS"]
orfanage.cds.unique <- unique(orfanage.cds)
        
# Annotate the Orfanage CDS as known/novel
matches <- findOverlaps(orfanage.cds.unique, gencode.cds.unique, type = "equal")
novelty_labels <- rep("novel", length(orfanage.cds.unique))
novelty_labels[queryHits(matches)] <- "known"
mcols(orfanage.cds.unique)$novelty <- novelty_labels

# Split into known/novel
known <- orfanage.cds.unique[orfanage.cds.unique$novelty == "known"]
novel <- orfanage.cds.unique[orfanage.cds.unique$novelty == "novel"]

# Trim the novel CDS to remove overlap with gencode CDS
overlap <- findOverlaps(novel, gencode.cds.unique)
novel.overlap <- extractList(gencode.cds.unique, as(overlap, "List"))
novel.subset <- psetdiff(novel, novel.overlap) %>% unlist()

# For each object, remove rows fully contained within other rows to avoid double-counting

# Known: find all of the hits within the object that aren't from matching itself
hits <- findOverlaps(known, known, type = "within")
hits <- hits[queryHits(hits) != subjectHits(hits)]

# Known: remove the hits from the GRanges object
contained.index <- unique(queryHits(hits))
known.subset <- known[-contained.index]

# Novel: find all of the hits within the object that aren't from matching itself
hits <- findOverlaps(novel.subset, novel.subset, type = "within")
hits <- hits[queryHits(hits) != subjectHits(hits)]

# Novel: remove the hits from the GRanges object
contained.index <- unique(queryHits(hits))
novel.subset <- novel.subset[-contained.index]

# Add riboseq information to the Orfanage object
known.subset$n_riboseq <- countOverlaps(known.subset, ribo)
known.subset$avg_riboseq <- known.subset$n_riboseq / width(known.subset)

novel.subset$n_riboseq <- countOverlaps(novel.subset, ribo)
novel.subset$avg_riboseq <- novel.subset$n_riboseq / width(novel.subset)

# Add a boolean validation column using 0.04 as the cutoff
known.subset$validated <- ifelse(known.subset$avg_riboseq > 0.04, T, F)
novel.subset$validated <- ifelse(novel.subset$avg_riboseq > 0.04, T, F)

# Create summary tables for graphing
exon.known.summary_df <- as.data.frame(mcols(known.subset)) %>%
    count(validated, name = "len")

exon.novel.summary_df <- as.data.frame(mcols(novel.subset)) %>%
    count(validated, name = "len")

SJ.known.summary_df <- LR_SJ %>%
    filter(LR, GENCODE) %>%
    count(SR, name = "len") %>%
    rename(
        "SR" = "validated"
    )

SJ.novel.summary_df <- LR_SJ %>%
    filter(LR, !GENCODE) %>%
    count(SR, name = "len") %>%
    rename(
        "SR" = "validated"
    )

known.summary_df <- rbind(exon.known.summary_df, SJ.known.summary_df) %>% 
    group_by(validated) %>%
    summarise(
        len = sum(len)
    ) %>% 
    mutate(
        percent = len / sum(len) * 100
    ) %>%     
    arrange(desc(validated)) %>%
    mutate(ypos = cumsum(percent) -0.5 * percent)

novel.summary_df <- rbind(exon.novel.summary_df, SJ.novel.summary_df) %>% 
    group_by(validated) %>%
    summarise(
        len = sum(len)
    ) %>% 
    mutate(
        percent = len / sum(len) * 100
    ) %>%     
    arrange(desc(validated)) %>%
    mutate(ypos = cumsum(percent) -0.5 * percent)        

#--------------Plotting---------------------------

my_theme <- theme_void() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)
    )

theme_set(my_theme)

p1 <- known.summary_df %>%
    ggplot(aes(x = "", y = percent, fill = validated)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    labs(fill = "Validation by\n Ribo-seq\n splice junctions") +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 6)

p2 <- novel.summary_df %>%
    ggplot(aes(x = "", y = percent, fill = validated)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    labs(fill = "Validation by\n Ribo-seq\n splice junctions") +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 6)

p1 + p2
ggsave("figures/supplementary/Ribo-seq_pie.pdf", width = 200, height = 100, units ="mm")