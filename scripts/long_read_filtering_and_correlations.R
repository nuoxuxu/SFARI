library(tidyverse)
library(stringr)
library(plyr)
library(tidyr)
library(wesanderson)
library(cowplot)
library(beepr)
library(data.table)
library(reshape2)
library(gridExtra)
library(purrr)
library(gplots)
library(RColorBrewer)
library(GenomicRanges)

setDTthreads(threads = 192)
#Gets first and last exons
firs_last_fun <- function(reads){
  reads_pos <- reads[reads$strand=='+',]
  reads_neg <- reads[reads$strand=='-',]
  
  reads_all_pos_first <- reads_pos[reads_pos[, .I[start == min(start)], by=read_name]$V1] %>% as.data.frame()
  reads_all_pos_last <- reads_pos[reads_pos[, .I[start == max(start)], by=read_name]$V1] %>% as.data.frame()
  
  reads_all_neg_first <- reads_neg[reads_neg[, .I[end == max(end)], by=read_name]$V1] %>% as.data.frame()
  reads_all_neg_last <- reads_neg[reads_neg[, .I[end == min(end)], by=read_name]$V1] %>% as.data.frame()
  
  reads_all_first <- rbind(reads_all_pos_first,reads_all_neg_first)
  reads_all_last <- rbind(reads_all_pos_last,reads_all_neg_last)
  
  return(list(reads_all_first=reads_all_first,reads_all_last=reads_all_last))
}

#Runs intersect
intersect_function <- function(reads_all_list){
  #Split starts and ends from the classification output
  read_starts <- reads_all_list[["reads_all_first"]][,c('chr','start','end','read_name','gene_id','strand')]
  read_ends <- reads_all_list[["reads_all_last"]][,c('chr','start','end','read_name','gene_id','strand')]
  gr_starts <- GRanges(
    seqnames = read_starts$chr,
    ranges = IRanges(start=read_starts$start+1, end=read_starts$end),
    strand = read_starts$strand,
    name = read_starts$read_name,
    gene_name = read_starts$gene_id
  )
  gr_ends <- GRanges(
    seqnames = read_ends$chr,
    ranges = IRanges(start=read_ends$start+1, end=read_ends$end),
    strand = read_ends$strand,
    name = read_ends$read_name,
    gene_name = read_ends$gene_id
  )  
  polyA_sites <- fread("data/atlas.clusters.2.0.GRCh38.96.bed",fill = T,sep = '\t') %>% 
    as.data.frame()
  gr_polyA <- GRanges(
    seqnames = polyA_sites$V1,
    ranges = IRanges(start=polyA_sites$V2, end=polyA_sites$V3),
    strand = polyA_sites$V6
  )
  # Find all overlaps between full last exon ranges and polyA clusters
  hits <- findOverlaps(gr_ends, gr_polyA, ignore.strand = FALSE)

  # Compute the width of each pairwise intersection
  overlap_widths <- width(pintersect(
      gr_ends[queryHits(hits)],
      gr_polyA[subjectHits(hits)]
  ))
  
  # Build a data frame of all hits with their overlap widths
  hits_df <- data.frame(
      exon_idx    = queryHits(hits),
      cluster_idx = subjectHits(hits),
      overlap     = overlap_widths
  )

  # For each last exon, keep only the polyA cluster with the largest overlap
  best_hits <- hits_df |>
      dplyr::group_by(exon_idx) |>
      dplyr::slice_max(overlap, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
  
  # Run bedtools to extract the sequence
  system("bedtools intersect -s -wao -a read_start_tmp.bed -b first_exonstmp.bed > starts_intersect.bed")
  system("bedtools intersect -s -wao -a read_end_tmp.bed -b polyAtmp.bed > ends_intersect.bed")
  
  starts_intersect <- fread('starts_intersect.bed',sep = '\t') %>% as.data.frame()
  ends_intersect <- fread('ends_intersect.bed',sep = '\t') %>% as.data.frame()
  
  #Now classify reads
  reads_no_start <- subset(starts_intersect,V10=='.') #Don't match first
  reads_no_end <- subset(ends_intersect,V10=='.') #Don't match polyA
  reads_matching_start <- starts_intersect[starts_intersect$V5==starts_intersect$V10,] #Match first
  reads_matching_end <- ends_intersect[ends_intersect$V5==ends_intersect$V10,] #Match polyA
  reads_matching_end <- setDT(reads_matching_end)[reads_matching_end[, .I[V11 == max(V11)], by=V4]$V1] %>% as.data.frame() #Keep the most downstream polyA
  
  #Merge reads using their names if they passed CAGE and polyA filters
  reads_passed_filters <- merge(reads_matching_start[,c('V1','V2','V3','V4','V5','V6','V11')],reads_matching_end[,c('V2','V3','V4','V10','V11')],by='V4')
  reads_passed_filters <- reads_passed_filters[reads_passed_filters$V5==reads_passed_filters$V10,]#Using gene ids for overlapping genes. I'm not removing reads mapping to multiple genes because only the most up and down features PER READ will be used when splitting, thus, there is no need to do that
  colnames(reads_passed_filters) <- c('read_name','chr','start_FE','end_FE','gene_id_FE','strand','FE_index','start_LE','end_LE','gene_id_LE','polyA_index')
  
  
  return(list(reads_no_start=reads_no_start,reads_no_end=reads_no_end,reads_passed_filters=reads_passed_filters))
  
}

#####
#Initiate variables
spearman_per_sample <- list()
gene_list_no_FE <- list()
gene_list_no_polyA <- list()
reads_across_datasets <- list()
reads_across_datasets_3utr <- list()
truncation_metrics_list <- list()
file_counter <- 1

bed_files <- list.files(path = 'nextflow_results/alternative_start_codon',pattern = '*.bed')
for (i in 1:length(bed_files)){
  reads_all <- fread("nextflow_results/alternative_start_codon/CN_1_3.aligned.bed",fill = T,sep = '\t') %>% 
    as.data.frame() %>% 
    select(-V5)
  colnames(reads_all) <- c('chr','start','end','read_name','strand','gene_id')
  
  reads_all$end <- as.numeric(reads_all$end)
  reads_all$start <- as.numeric(reads_all$start)
  
  #Split read starts and ends to make bed file 
  reads_all_list <- firs_last_fun(reads = setDT(reads_all))
  
  #Run intersect and store outputs
  intersect_function_out <- intersect_function(reads_all_list)
  gene_list_no_FE[[file_name]] <- intersect_function_out[["reads_no_start"]] #Store reads not matching peaks
  gene_list_no_polyA[[file_name]] <- intersect_function_out[["reads_no_end"]]
  reads_passed_filters <- intersect_function_out[["reads_passed_filters"]]
  
  #Get filter metrics
  counts_per_gene_raw <- table(reads_all$V10) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "freq_no_filter"))
  terminal_truncation <- table(reads_across_datasets_3utr[[file_name]]$gene_id) %>% as.data.frame()   %>% `colnames<-`(c("gene_id", "freq_terminal_reads"))
  truncation_5 <- table(reads_no_start$V5) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "freq_5end_filter"))
  truncation_3 <- table(reads_no_end$V5) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "freq_3end_filter"))
  counts_per_gene_after_filters <- table(reads_passed_filters$gene_id_FE) %>% as.data.frame() %>% `colnames<-`(c("gene_id", "final_freq"))
  
  truncation_metrics <- merge(merge(merge(merge(counts_per_gene_raw,terminal_truncation,by='gene_id',all.x=T),truncation_5,by='gene_id',all.x=T),truncation_3,by='gene_id',all.x=T),counts_per_gene_after_filters,,by='gene_id',all.x=T)
  
  truncation_metrics_list[[file_name]] <- truncation_metrics
  
  #Make vectors of genes using AFE & ALE or Unique F | L exons
  #Keep reads with TSS and TES that are used at least three times for each gene
  reads_passed_filters <- reads_passed_filters %>% dplyr::group_by(gene_id_FE,FE_index) %>% dplyr::mutate(TSS_counts=n())%>% filter(TSS_counts>2)
  reads_passed_filters <- reads_passed_filters %>% dplyr::group_by(gene_id_FE,polyA_index) %>% dplyr::mutate(polyA_counts=n()) %>% filter(polyA_counts>2) %>% as.data.frame()
  #Remove reads that match more than one gene
  reads_to_use_uniquelymapped <- table(reads_passed_filters$read_name) %>% as.data.frame() %>% filter(Freq==1) %>% pull(Var1)
  reads_passed_filters <- reads_passed_filters[reads_passed_filters$read_name%in%reads_to_use_uniquelymapped,]
  #Store on a list
  reads_across_datasets[[file_name]] <- reads_passed_filters
  
  data_count_1 <- aggregate(data = reads_passed_filters,FE_index ~ gene_id_FE,function(x) length(unique(x)))
  data_count_1 <- merge(data_count_1,aggregate(data = reads_passed_filters,polyA_index ~ gene_id_FE,function(x) length(unique(x))),by='gene_id_FE')
  
  alt_genes <- subset(data_count_1,FE_index>1 & polyA_index>1) %>% pull(gene_id_FE)
  unique_genes <- subset(data_count_1,FE_index==1 | polyA_index==1) %>% pull(gene_id_FE)
  
  #Calculate spearman Rs
  spearman_all <- by(reads_passed_filters, reads_passed_filters$gene_id_FE, FUN = function(X) cor(X$start_FE, X$end_LE, method = "spearman")) #From the start of the first exon to the end of the last
  spearman_all <- data.frame(gene_id = dimnames(spearman_all)[[1]],corr = as.vector(spearman_all))

  #Assign alternative/unique exons
  spearman_all$exon_type[spearman_all$gene_id %in% alt_genes] <- 'Alternative FE & polyA'
  spearman_all$exon_type[spearman_all$gene_id %in% unique_genes] <- 'Unique FE or polyA'
  
  #Get counts used for Spearman
  reads_all_counts <- counts_per_gene_after_filters %>% dplyr::mutate(cpm=10^6*final_freq/sum(final_freq))
  spearman_all <- merge(spearman_all,reads_all_counts,by = 'gene_id') #Add counts (after filtering) information
  
  #Add useful information
  spearman_all <- merge(spearman_all,data_count_1,by.x='gene_id',by.y='gene_id_FE')
  
  colnames(spearman_all)[c(6,7)] <- c('number_of_used_FE','number_of_used_polyA')
  spearman_per_sample[[file_name]] <- data.frame(spearman_all,sample=file_name)
}

#save.image('~/Dropbox (UMass Medical School)/PaiLab/members/Ezequiel/3rd year/TSS-TES/ENCODE_full/ENCODE_HITindex_intersect.RData')
#names(spearman_per_sample) <- gsub("\\..*","",names(spearman_per_sample))

spearman_per_sample_df[is.na(spearman_per_sample_df$corr),'corr'] <- 0 #NAs are things without any correlation so I reassign as 0

save(spearman_per_sample,spearman_per_sample_df,truncation_metrics_list,reads_across_datasets,gene_list_no_FE,gene_list_no_polyA,reads_across_datasets_3utr,file = '~/Dropbox (UMass Medical School)/PaiLab/members/Ezequiel/3rd year/TSS-TES/ENCODE_full/ENCODE_HITindex_intersect.RData')

load('~/Dropbox (UMass Medical School)/PaiLab/members/Ezequiel/3rd year/TSS-TES/ENCODE_full/ENCODE_HITindex_intersect.RData')

write.table(spearman_per_sample_df,'~/Dropbox (UMass Medical School)/PaiLab/members/Ezequiel/3rd year/TSS-TES/ENCODE_full/spearman_ENCODE_HITindex_intersect.tsv',sep = '\t',col.names = T,row.names = F,quote = F)