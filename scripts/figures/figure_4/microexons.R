### MICROEXONS ###
library(dplyr)
library(tidyr)
library(arrow)
library(purrr)
library(edgeR)
library(IsoformSwitchAnalyzeR)


setwd("")



classification <- read_parquet("./data/final_classification.parquet")

#Subset the transcript file to these PB IDs. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)

classif_final_tr <- read_parquet("./data/final_classification.parquet")
classification <- classif_final_tr

filtered_tr_count <- tr_count[tr_count$isoform %in% classification$isoform, ] 




# ISOFORM SWITCHING MICs --------------------------------------------------
###IsoformSwitches
aSwitchList_part2 <- readRDS("./code/IsoformSwitchAnalyzeR/output/isoformswitch_part2_v4.rds")
aSwitchList_part2[["isoformFeatures"]]$gene_name <- aSwitchList_part2[["isoformFeatures"]]$gene_id


MIC <- aSwitchList_part2$AlternativeSplicingAnalysis
MIC <- MIC %>% dplyr::select(isoform_id, 
                             ES, ES_genomic_start, ES_genomic_end, 
                             MES, MES_genomic_start, MES_genomic_end)
MIC$MES_genomic_start <- gsub(",", ";", MIC$MES_genomic_start)
MIC$MES_genomic_end <- gsub(",", ";", MIC$MES_genomic_end)


# Function to pair starts and ends into rows
extract_events <- function(df, start_col, end_col, event_type) {
  df %>%
    mutate(
      start_list = strsplit(as.character(.data[[start_col]]), ";"),
      end_list = strsplit(as.character(.data[[end_col]]), ";")
    ) %>%
    # Keep only rows where both start and end are not NA
    filter(!is.na(.data[[start_col]]), !is.na(.data[[end_col]])) %>%
    # Pair them using map2 and turn into data frame
    mutate(paired = map2(start_list, end_list, ~ tibble(genomic_start = .x, genomic_end = .y))) %>%
    dplyr::select(isoform_id, paired) %>%
    unnest(paired) %>%
    group_by(isoform_id) %>%
    mutate(
      event_id = row_number(),
      event_type = event_type
    ) %>%
    ungroup()
}

# Apply the function for MES and ES
mes_long <- extract_events(MIC, "MES_genomic_start", "MES_genomic_end", "MES")
es_long  <- extract_events(MIC, "ES_genomic_start", "ES_genomic_end", "ES")

# Combine both
all_events <- bind_rows(mes_long, es_long)


all_events$exon_length <- (as.numeric(all_events$genomic_end) - as.numeric(all_events$genomic_start)) + 1 

all_events$MIC <- ifelse(all_events$exon_length <= 27, "MIC", "not_MIC") 

MIC_transcripts <- all_events %>% filter(MIC == "MIC") #1236 transcripts 
MIC_events <- MIC_transcripts %>% distinct(genomic_start, genomic_end, .keep_all = T) #281 unique events


###NOTE: I AM NOT REMOVING ANYTHING FROM MES OR ES. I AM ADDING A NEW 'SUBCATEGORY'
AS_full <- aSwitchList_part2$AlternativeSplicingAnalysis

#Collapse MICs
# Step 1: Aggregate by isoform_id and event_type
collapsed <- MIC_transcripts %>%
  group_by(isoform_id) %>%
  summarise(
    MIC = n(),
    MIC_genomic_start = paste(genomic_start, collapse = ";"),
    MIC_genomic_end = paste(genomic_end, collapse = ";"),
    .groups = "drop"
  )

AS_full <- merge(AS_full, collapsed[, c("isoform_id", "MIC", "MIC_genomic_start", "MIC_genomic_end")], 
                        by.x = "isoform_id", by.y = "isoform_id", 
                        all.x = TRUE)


# Replace NAs with 0 
AS_full <- AS_full %>%
  mutate(across(c(MIC), ~ replace_na(., 0)))

# Add this back to isoformswitch_part2_v2
aSwitchList_part2$AlternativeSplicingAnalysis <- AS_full


saveRDS(aSwitchList_part2, file="./code/microexons/output/isoformswitch_MIC.rds") 




# COMPLETE ISOFORM SWITCH LIST --------------------------------------------
#Re-run all IsoformSwitchAnalysis, but only for AS

PB_counts <- as.matrix(filtered_tr_count[,c(2:16)])
rownames(PB_counts) <- filtered_tr_count$isoform


dge_PB <- DGEList(counts = PB_counts, genes = filtered_tr_count[,1])
dge_PB <- calcNormFactors(dge_PB)


PB_tr_cpm_norm <- cpm(dge_PB, normalized.lib.sizes=TRUE, log=FALSE)
PB_tr_cpm_norm <- as.data.frame(PB_tr_cpm_norm)
PB_tr_cpm_norm$pb_id <- row.names(PB_tr_cpm_norm)

myDesign <- data.frame(sampleID = colnames(tr_count)[-1])
myDesign$condition <- c(rep("t00", 3), rep("t04", 6), rep("t30", 6))

#Prep expression matrix
colnames(tr_count)[1] <- "isoform_id"
colnames(PB_tr_cpm_norm)[16] <- "isoform_id"
tr_CPM <- PB_tr_cpm_norm %>% relocate(isoform_id, .before = iPSC_1)


#Replace gene_id with gene_name:
transcripts <- rtracklayer::import("./data/final_transcripts.gtf")
transcripts <- subset(transcripts, transcript_id %in% classification$isoform) #182,371
transcript_ids <- data.frame(pb_id = transcripts@elementMetadata@listData[["transcript_id"]])
transcript_ids <- merge(transcript_ids, classification[, c("isoform", "associated_gene")], 
                        by.x = "pb_id", by.y = "isoform", 
                        all.x = TRUE, sort = F)

# Direct assignment to the listData
transcripts@elementMetadata@listData[["gene_id"]] <- as.character(transcript_ids$associated_gene)
rtracklayer::export(transcripts, './code/microexons/input/final_transcripts_gene_id_replaced_gtf.gtf')

colnames(filtered_tr_count)[1] <- "isoform_id"
#Make a switch list
aSwitchList_all_tr <- importRdata(
  isoformCountMatrix   = filtered_tr_count,
  isoformRepExpression = tr_CPM,
  designMatrix         = myDesign,
  isoformExonAnnoation = ('./code/microexons/input/final_transcripts_gene_id_replaced_gtf.gtf'),             
 # isoformNtFasta       = ('./data/PacBio_Nuo/IsoformSwitch/full_sequences_fixed_head.fasta.gz'),
  addAnnotatedORFs     = FALSE,
  detectUnwantedEffects = F, 
  fixStringTieAnnotationProblem = F,
  showProgress = T)


aSwitchList_all_tr <- analyzeAlternativeSplicing(aSwitchList_all_tr, 
                                             onlySwitchingGenes=F, #only analyze isoform switches (as indicated by alpha and dIF below)
                                             alpha=1, #FDR corrected pval
                                             dIFcutoff = 0, #min change in abs. isoform usage
                                             showProgress=TRUE)
saveRDS(aSwitchList_all_tr, file="./code/microexons/output/isoformswitch_all_tr.rds") 


##Read in all
aSwitchList_all_tr <- readRDS("./code/microexons/output/isoformswitch_all_tr.rds")


all_AS <- read.table("./code/microexons/input/vastDB/EVENT_INFO-hg38.tab", 
                     header = T, 
                     sep = "\t") #downloaded from vastDB
all_AS_EX <- all_AS %>% filter(grepl("HsaEX[0-9]+", EVENT))
remove(all_AS)
vast_MIC_EX <- all_AS_EX %>% filter(LE_n <= 28) #4191


MIC <- aSwitchList_all_tr$AlternativeSplicingAnalysis
MIC <- MIC %>% dplyr::select(isoform_id, 
                             ES, ES_genomic_start, ES_genomic_end, 
                             MES, MES_genomic_start, MES_genomic_end)
MIC$MES_genomic_start <- gsub(",", ";", MIC$MES_genomic_start)
MIC$MES_genomic_end <- gsub(",", ";", MIC$MES_genomic_end)




# Function to pair starts and ends into rows
extract_events <- function(df, start_col, end_col, event_type) {
  df %>%
    mutate(
      start_list = strsplit(as.character(.data[[start_col]]), ";"),
      end_list = strsplit(as.character(.data[[end_col]]), ";")
    ) %>%
    # Keep only rows where both start and end are not NA
    filter(!is.na(.data[[start_col]]), !is.na(.data[[end_col]])) %>%
    # Pair them using map2 and turn into data frame
    mutate(paired = map2(start_list, end_list, ~ tibble(genomic_start = .x, genomic_end = .y))) %>%
    dplyr::select(isoform_id, paired) %>%
    unnest(paired) %>%
    group_by(isoform_id) %>%
    mutate(
      event_id = row_number(),
      event_type = event_type
    ) %>%
    ungroup()
}

# Apply the function for MES and ES
mes_long <- extract_events(MIC, "MES_genomic_start", "MES_genomic_end", "MES")
es_long  <- extract_events(MIC, "ES_genomic_start", "ES_genomic_end", "ES")

# Combine both
all_events <- bind_rows(mes_long, es_long)


all_events$exon_length <- (as.numeric(all_events$genomic_end) - as.numeric(all_events$genomic_start)) + 1 


all_events$MIC <- ifelse(all_events$exon_length <= 27, "MIC", "not_MIC")

MIC_transcripts <- all_events %>% filter(MIC == "MIC") #5795 transcripts
MIC_events <- MIC_transcripts %>% distinct(genomic_start, genomic_end, .keep_all = T) #483 unique events


#Find how many of these 483 are in vastDB??
MIC_events <- merge(MIC_events, classification[, c("isoform", "chrom", "associated_gene")], 
                        by.x = "isoform_id", by.y = "isoform", 
                        all.x = TRUE)
MIC_events$event <- paste0(MIC_events$chrom, ":", MIC_events$genomic_start, "-", MIC_events$genomic_end)

intersect(vast_MIC_EX$CO_A, MIC_events$event) #293!!!!!



# MICROEXONS IN GENCODE - BY LENGTH ---------------------------------------
gencode <- rtracklayer::import("./data/Gencode_v47/gencode.v47.annotation.gtf")

#To make it a fair comparison, only using genes we have in our data:
genes <- unique(classification$associated_gene)
gencode <- subset(gencode, gene_name %in% genes)

gencode_exon <- subset(gencode, type == "exon")
gencode_MIC <- subset(gencode_exon, width <28) #Extract all microexons (by length)
gencode_MIC_df <- as.data.frame(gencode_MIC)

gencode_unique_MIC <- gencode_MIC_df %>% distinct(seqnames, start, end,
                                                  .keep_all = T) #Count unique ones (by coordinates)

#9,774 unique MICs 
gencode_unique_MIC$MIC_coord <- paste0(gencode_unique_MIC$seqnames, ":",
                                       gencode_unique_MIC$start, "-", 
                                       gencode_unique_MIC$end)



# LONG-READ MICROEXONS - BY LENGTH ----------------------------------------
transcripts <- rtracklayer::import("./data/final_transcripts.gtf")
transcripts <- subset(transcripts, transcript_id %in% classification$isoform) #182,371

PB_exon <- subset(transcripts, type == "exon")
PB_MIC <- subset(PB_exon, width <28) #Extract all microexons (by length)
PB_MIC_df <- as.data.frame(PB_MIC)

PB_unique_MIC <- PB_MIC_df %>% distinct(seqnames, start, end,
                                                  .keep_all = T) #Count unique ones (by coordinates)

#1930 unique PB MICs
PB_unique_MIC$MIC_coord <- paste0(PB_unique_MIC$seqnames, ":",
                                  PB_unique_MIC$start, "-", 
                                  PB_unique_MIC$end)



# MICROEXONS INTERSECTION -------------------------------------------------
all_AS <- read.table("./code/microexons/input/vastDB/EVENT_INFO-hg38.tab", 
                     header = T, 
                     sep = "\t")
all_AS_EX <- all_AS %>% filter(grepl("HsaEX[0-9]+", EVENT))
remove(all_AS)
vast_MIC_EX <- all_AS_EX %>% filter(LE_n <= 28) #4191
vast_MIC_EX <- vast_MIC_EX %>% filter(GENE %in% genes) #3238

#Gencode MICs and vastDB MICs
length(intersect(gencode_unique_MIC$MIC_coord, vast_MIC_EX$CO_A)) #1189 are in both
length(setdiff(gencode_unique_MIC$MIC_coord, vast_MIC_EX$CO_A)) #8585 in gencode, not PB
length(setdiff(vast_MIC_EX$CO_A, gencode_unique_MIC$MIC_coord)) #2049 in vast, not gencode

#Our MICs and Gencode MICs
length(intersect(gencode_unique_MIC$MIC_coord, PB_unique_MIC$MIC_coord)) #846 are in both
length(setdiff(gencode_unique_MIC$MIC_coord, PB_unique_MIC$MIC_coord)) #8928 in gencode, not PB
length(setdiff(PB_unique_MIC$MIC_coord, gencode_unique_MIC$MIC_coord)) #1084 in PB, not gencode

##Combine gencode and vast
combined_MIC <- union(gencode_unique_MIC$MIC_coord, vast_MIC_EX$CO_A)
combined_MIC <- unique(combined_MIC) #11823

#Our MICs vs all other MICs
length(intersect(combined_MIC, PB_unique_MIC$MIC_coord)) #884 are in all (40% of all our exons)
length(setdiff(combined_MIC, PB_unique_MIC$MIC_coord)) #10939 in gencode or vast, not PB
length(setdiff(PB_unique_MIC$MIC_coord, combined_MIC)) #1046 in PB, not gencode or vast

#884 'known' vs 1046 'novel'
known_MIC <- intersect(combined_MIC, PB_unique_MIC$MIC_coord) 
novel_MIC <- setdiff(PB_unique_MIC$MIC_coord, combined_MIC) 

#get exon coords back to tr mapping
PB_exon_df <- as.data.frame(PB_exon)
PB_exon_df$coord <- paste0(PB_exon_df$seqnames, ":",
                        PB_exon_df$start, "-", 
                        PB_exon_df$end)

all_exons <- PB_exon_df %>%
  mutate(MIC_type = case_when(
    coord %in% known_MIC ~ "known_MIC",
    coord %in% novel_MIC ~ "novel_MIC"))

MIC_exons_mapping <- all_exons %>%
  group_by(transcript_id) %>%
  drop_na(MIC_type) %>%
  ungroup()

