### IsoformSwitchAnalyzeR ###

library(arrow)
library(IsoformSwitchAnalyzeR)

setwd("")



final_pb_ids <- read.table("./data/filtered_transcript_list.txt", header = T)

#Subset the transcript file to these PB IDs. 
tr_count <- read_parquet("./data/final_expression.parquet")
tr_count <- as.data.frame(tr_count)

filtered_tr_count <- tr_count[tr_count$isoform %in% final_pb_ids$x, ]

classification <- read_parquet("./data/final_classification.parquet")

tr_CPM <- read.csv("./code/expression/output/transcript_CPM.csv")
gene_CPM <- read.csv("./code/expression/output/gene_CPM.csv")



# BUILD SWITCH LIST ----------------------------------------------------------
### Make design
myDesign <- data.frame(sampleID = colnames(filtered_tr_count)[-1])
myDesign$condition <- c(rep("t00", 3), rep("t04", 6), rep("t30", 6))

#Prep expression matrix
colnames(filtered_tr_count)[1] <- "isoform_id"
colnames(tr_CPM)[16] <- "isoform_id"
tr_CPM <- tr_CPM %>% relocate(isoform_id, .before = iPSC_1)


#Make a switch list
aSwitchList <- importRdata(
  isoformCountMatrix   = filtered_tr_count,
  isoformRepExpression = tr_CPM,
  designMatrix         = myDesign,
  isoformExonAnnoation = ('./code/IsoformSwitchAnalyzeR/input/ORF_gene_id_replaced_tr_exon.gtf.gz'), #gtf with only transcript and exon types            
  isoformNtFasta       = ('./code/IsoformSwitchAnalyzeR/input/final_transcripts.fasta.gz'), 
  addAnnotatedORFs     = FALSE, 
  detectUnwantedEffects = F, #FALSE because I'm importing my pre-calculated CPMs
  fixStringTieAnnotationProblem = F,
  showProgress = T)

# Step 1 of 10: Checking data...
# Step 2 of 10: Obtaining annotation...
# importing GTF (this may take a while)...
# Step 3 of 10: Fixing StringTie gene annoation problems...
# Was skipped as instructed via the "fixStringTieAnnotationProblem" argument...
# Step 4 of 10: Calculating expression estimates from count data...
# Skipped as user supplied expression via the "isoformRepExpression" argument...
# Step 5 of 10: Testing for unwanted effects...
# Skipped due to "detectUnwantedEffects=FALSE". 
# No unwanted effects added
# Step 6 of 10: Batch correcting expression estimates...
# Skipped as no batch effects were found or annoated...
# Step 7 of 10: Extracting data from each condition...
# =======================================================================| 100%
# Step 8 of 10: Making comparisons...
# =======================================================================| 100%
# Step 9 of 10: Making switchAnalyzeRlist object...
# Step 10 of 10: Guestimating differential usage...
# The GUESSTIMATED number of genes with differential isoform usage are:
#  comparison estimated_genes_with_dtu
# 1 t00 vs t04              2378 - 3964
# 2 t00 vs t30              2922 - 4871
# 3 t04 vs t30              1451 - 2418
# Done
# 
                 
aSwitchList$orfAnalysis #NULL
aSwitchList$isoformFeatures$PTC #NULL
# 

##Need to manually add gene_overall_value, gene_value_1, and gene_value_2
gene_CPM$t00 <- rowMeans(subset(gene_CPM, select = c(1:3)))
gene_CPM$t04 <- rowMeans(subset(gene_CPM, select = c(4:9)))
gene_CPM$t30 <- rowMeans(subset(gene_CPM, select = c(10:15)))

# Pivot longer to get condition-value pairs
df_long <- gene_CPM %>%
  pivot_longer(cols = c(t00, t04, t30), names_to = "condition", values_to = "value")

# Do pairwise combinations by joining the long data frame with itself
df_pairs <- df_long %>%
  inner_join(df_long, by = "gene_name") %>%
  filter(condition.x < condition.y) %>%
  dplyr::select(
    gene_name,
    condition_1 = condition.x,
    gene_value_1 = value.x,
    condition_2 = condition.y,
    gene_value_2 = value.y
  )


gene_CPM$gene_overall_mean <- rowMeans(subset(gene_CPM, select = c(1:15)))
aSwitchList$isoformFeatures$gene_overall_mean <- NULL

aSwitchList$isoformFeatures <- merge(aSwitchList$isoformFeatures, gene_CPM[,c('gene_name', 'gene_overall_mean')],
                                           by.x = c("gene_id"),
                                           by.y = c("gene_name"))

aSwitchList$isoformFeatures$gene_value_1 <- NULL
aSwitchList$isoformFeatures$gene_value_2 <- NULL
aSwitchList$isoformFeatures <- merge(aSwitchList$isoformFeatures, df_pairs[,c('gene_name', 'condition_1', 'condition_2', 'gene_value_1', 'gene_value_2')],
                                     by.x = c("gene_id", 'condition_1', 'condition_2'),
                                     by.y = c("gene_name", 'condition_1', 'condition_2'))


#Adding ORFanage annotated ORFs
aSwitchList <- addORFfromGTF(
  switchAnalyzeRlist     = aSwitchList,
  pathToGTF              = './code/IsoformSwitchAnalyzeR/input/ORF_gene_id_replaced_new.gtf.gz') #gtf with ORFanage predictions

# Step 1 of 2: importing GTF (this may take a while)...
# Step 2 of 2: Adding ORF...
# Added ORF info (incl info about isoforms annotated as not having an ORF) to 158844 isoforms.
# This correspond to 100% of isoforms in the switchAnalyzeRlist.
# Which includes NaN% of isoforms from annotated genes (novel isoforms not counted) in the switchAnalyzeRlist.
# Done.

#Pre-filtering isoforms
aSwitchList <- preFilter(
  switchAnalyzeRlist         = aSwitchList,
  geneExpressionCutoff       = 1,     # default (in FPKM/TPM/RPKM)
  isoformExpressionCutoff    = 0,     # default (removes completely unused isoforms)
  IFcutoff                   = 0.01,  # default (IF = isoform fraction - isoform must be used in at least 1 condition of a comparison)
  removeSingleIsoformGenes   = TRUE,  # default
  reduceToSwitchingGenes     = FALSE, # default 
  keepIsoformInAllConditions = T) # keep an isoform in all conditions, even if it only passes filters in some conditions

# The filtering removed 90043 ( 56.69% of ) transcripts. There is now 68801 isoforms left

saveRDS(aSwitchList, file = "./code/IsoformSwitchAnalyzeR/output/isoformswitch_v4.rds")



aSwitchList <- readRDS("./code/IsoformSwitchAnalyzeR/output/isoformswitch_v4.rds")



# DTE/DGE/DTU -------------------------------------------------------------
#DTE - already done in DESeq2 section
res_t04_vs_t00 <- read.csv("./code/expression/output/DESeq2_tr_t04_vs_t00.csv")

res_t30_vs_t00 <- read.csv("./code/expression/output/DESeq2_tr_t30_vs_t00.csv")

res_t30_vs_t04 <- read.csv("./code/expression/output/DESeq2_tr_t30_vs_t04.csv")

#Add gene_names
res_t04_vs_t00 <- merge(res_t04_vs_t00, classification[, c("isoform", "associated_gene", "structural_category")],
                        by.x = "pb_id", by.y = "isoform", all.x = T)
res_t04_vs_t00 <- res_t04_vs_t00 %>%
  mutate(padj = replace_na(padj, 1))
res_t30_vs_t00 <- merge(res_t30_vs_t00, classification[, c("isoform", "associated_gene", "structural_category")],
                        by.x = "pb_id", by.y = "isoform", all.x = T)
res_t30_vs_t00 <- res_t30_vs_t00 %>%
  mutate(padj = replace_na(padj, 1))
res_t30_vs_t04 <- merge(res_t30_vs_t04, classification[, c("isoform", "associated_gene", "structural_category")],
                        by.x = "pb_id", by.y = "isoform", all.x = T)
res_t30_vs_t04 <- res_t30_vs_t04 %>%
  mutate(padj = replace_na(padj, 1))

DTE <- res_t04_vs_t00[,c(1,3,6:9)]
colnames(DTE)[c(2:4)] <- c("log2FC_t04_v_t00", "pvalue_t04_v_t00", "padj_t04_v_t00")

colnames(res_t30_vs_t00)[c(3,6,7)] <- c("log2FC_t30_v_t00", "pvalue_t30_v_t00", "padj_t30_v_t00")
colnames(res_t30_vs_t04)[c(3,6,7)] <- c("log2FC_t30_v_t04", "pvalue_t30_v_t04", "padj_t30_v_t04")


DTE <- merge(DTE, res_t30_vs_t00[, c("pb_id",
                                               "log2FC_t30_v_t00",
                                               "pvalue_t30_v_t00",
                                               "padj_t30_v_t00")],
                  by.x = "pb_id", by.y = "pb_id", all.x = T)

DTE <- merge(DTE, res_t30_vs_t04[, c("pb_id",
                                               "log2FC_t30_v_t04",
                                               "pvalue_t30_v_t04",
                                               "padj_t30_v_t04")],
                  by.x = "pb_id", by.y = "pb_id", all.x = T)
DTE_results <- DTE
colnames(DTE_results)[1] <- "isoform_id"

#DGE
gene_sum_t04_vs_t00_df <- read.csv("./code/expression/output/DESeq2_gene_sum_t04_vs_t00.csv")
gene_sum_t04_vs_t00_df <- gene_sum_t04_vs_t00_df %>%
  mutate(padj = replace_na(padj, 1))

gene_sum_t30_vs_t00_df <- read.csv("./code/expression/output/DESeq2_gene_sum_t30_vs_t00.csv")
gene_sum_t30_vs_t00_df <- gene_sum_t30_vs_t00_df %>%
  mutate(padj = replace_na(padj, 1))

gene_sum_t30_vs_t04_df <- read.csv("./code/expression/output/DESeq2_gene_sum_t30_vs_t04.csv")
gene_sum_t30_vs_t04_df <- gene_sum_t30_vs_t04_df %>%
  mutate(padj = replace_na(padj, 1))


DGE <- gene_sum_t04_vs_t00_df[,c(2, 5:7)]
colnames(DGE)[c(1:3)] <- c("log2FC_t04_v_t00", "pvalue_t04_v_t00", "padj_t04_v_t00")

colnames(gene_sum_t30_vs_t00_df)[c(2,5,6)] <- c("log2FC_t30_v_t00", "pvalue_t30_v_t00", "padj_t30_v_t00")
colnames(gene_sum_t30_vs_t04_df)[c(2,5,6)] <- c("log2FC_t30_v_t04", "pvalue_t30_v_t04", "padj_t30_v_t04")


DGE <- merge(DGE, gene_sum_t30_vs_t00_df[, c("gene_name",
                                             "log2FC_t30_v_t00",
                                             "pvalue_t30_v_t00",
                                             "padj_t30_v_t00")],
                           by.x = "gene_name", by.y = "gene_name", all.x = T)

DGE <- merge(DGE, gene_sum_t30_vs_t04_df[, c("gene_name",
                                             "log2FC_t30_v_t04",
                                             "pvalue_t30_v_t04",
                                             "padj_t30_v_t04")],
                           by.x = "gene_name", by.y = "gene_name", all.x = T)
DGE_results <- DGE
colnames(DGE_results)[1] <- "gene_id"

#DTU
aSwitchList_part1 <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist         = aSwitchList,
    alpha = 0.05, #default
    dIFcutoff = 0.05, #CHANGED
    reduceToSwitchingGenes     = FALSE,
    reduceFurtherToGenesWithConsequencePotential = FALSE, 
    showProgress = T)
  
# Step 1 of 2: Testing each pairwise comparisons with DEXSeq (this might be a bit slow)...
# Estimated run time is: 115.7 min
# ======================================================================| 100%
# Step 2 of 2: Integrating result into switchAnalyzeRlist...
# Isoform switch analysis was performed for 31587 gene comparisons (100%).
# Total runtime: 10.23 min
# Done

                         
### Add DTE/DGE to switchList
#For DTE
DTE_results <- DTE_results %>%
  pivot_longer(
  cols = starts_with("padj"),           # Select columns that start with "padj"
  names_to = c("condition_1", "condition_2"), # Split the column names into two parts: condition_1 and condition_2
  names_sep = "_v_",                     # Split the column name by "_v_" to get the pair of conditions
  values_to = "padj_value"              # Put the padj values into a column named "padj_value"
  )

DTE_results$condition_1 <- gsub("padj_", "", DTE_results$condition_1)

colnames(DTE_results)[10:11] <- c("condition_2", "condition_1")

aSwitchList_part1$isoformFeatures <- merge(aSwitchList_part1$isoformFeatures, DTE_results[,c('isoform_id', 'condition_1', 'condition_2', 'padj_value')],
                   by.x = c("isoform_id", "condition_1", "condition_2"),
                   by.y = c("isoform_id", "condition_1", "condition_2"))
aSwitchList_part1$isoformFeatures$iso_q_value <- NULL
colnames(aSwitchList_part1$isoformFeatures)[colnames(aSwitchList_part1$isoformFeatures) == "padj_value"] <- "iso_q_value"

#For DGE
DGE_results <- DGE_results %>%
  pivot_longer(
    cols = starts_with("padj"),           # Select columns that start with "padj"
    names_to = c("condition_1", "condition_2"), # Split the column names into two parts: condition_1 and condition_2
    names_sep = "_v_",                     # Split the column name by "_v_" to get the pair of conditions
    values_to = "padj_value"              # Put the padj values into a column named "padj_value"
  )

DGE_results$condition_1 <- gsub("padj_", "", DGE_results$condition_1)

colnames(DGE_results)[8:9] <- c("condition_2", "condition_1")

aSwitchList_part1$isoformFeatures <- merge(aSwitchList_part1$isoformFeatures, DGE_results[,c('gene_id', 'condition_1', 'condition_2', 'padj_value')],
                                           by.x = c("gene_id", "condition_1", "condition_2"),
                                           by.y = c("gene_id", "condition_1", "condition_2"))
aSwitchList_part1$isoformFeatures$gene_q_value <- NULL
colnames(aSwitchList_part1$isoformFeatures)[colnames(aSwitchList_part1$isoformFeatures) == "padj_value"] <- "gene_q_value"




###Save
saveRDS(aSwitchList_part1, file = "./code/IsoformSwitchAnalyzeR/output/isoformswitch_part1_v4.rds")

aSwitchList_part1 <- readRDS("./code/IsoformSwitchAnalyzeR/output/isoformswitch_part1_v4.rds")


summary(aSwitchList_part1)

### Export DTU results
supp_table <- aSwitchList_part1$isoformSwitchAnalysis %>%
  dplyr::rename(
    DTU_dIF    = "dIF",
    DTU_pval   = "pvalue",
    DTU_qval   = "padj")

write.csv(supp_table, "./code/IsoformSwitchAnalyzeR/output/DTU_table.csv")



# ALTERNATIVE SPLICING & FUNC CONSEQUENCES-------------------------------------
aSwitchList_part1 <- readRDS("./code/IsoformSwitchAnalyzeR/output/isoformswitch_part1_v4.rds")

###Setting all FDR < 0.05, and dIF = 0.05

#Analyze alternative splicing
aSwitchList_part1 <- analyzeAlternativeSplicing(aSwitchList_part1, 
                                                onlySwitchingGenes=F, #get splicing for all isoforms
                                                alpha=0.05, #FDR corrected pval
                                                dIFcutoff = 0.05, #min change in abs. isoform usage
                                                showProgress=TRUE)

#Add PFAM, IUPRED, SignalP data, DeepLoc2
#Add PFAM data
aSwitchList_part2 <- analyzePFAM(
  switchAnalyzeRlist   = aSwitchList_part1,
  pathToPFAMresultFile = "./code/IsoformSwitchAnalyzeR/input/full_pfam_scan_results_isoform_edited.txt",
  showProgress=T)

#Add IUPRED data
aSwitchList_part2 <- analyzeIUPred2A(
  switchAnalyzeRlist        = aSwitchList_part2,
  pathToIUPred2AresultFile  = "./code/IsoformSwitchAnalyzeR/input/iupred2a_processed_result.txt", 
  showProgress = T)

#Add SignalP data
aSwitchList_part2 <- analyzeSignalP(
  switchAnalyzeRlist = aSwitchList_part2,
  pathToSignalPresultFile = "./code/IsoformSwitchAnalyzeR/input/processed.signalp5",
  minSignalPeptideProbability = 0.5,
  quiet=FALSE
) 

#Add DeepLoc data
aSwitchList_part2 <- analyzeDeepLoc2(
  switchAnalyzeRlist = aSwitchList_part2,
  pathToDeepLoc2resultFile = "./v4/deeploc2.csv",
  enforceProbabilityCutoff = TRUE,
  probabilityCutoff = NULL,
  quiet = FALSE
)


#Extract sequences
aSwitchList_part2 <- extractSequence(
  aSwitchList_part2,
  genomeObject = NULL,
  onlySwitchingGenes = F,
  alpha = 0.0,
  dIFcutoff = 0.0,
  extractNTseq = TRUE,
  extractAAseq = TRUE,
  removeShortAAseq = F,
  removeLongAAseq = FALSE,
  alsoSplitFastaFile = FALSE,
  removeORFwithStop=F,
  addToSwitchAnalyzeRlist = TRUE,
  writeToFile = F,
  quiet=FALSE
)

#Save temp file:
#saveRDS(aSwitchList_part2, file = "./isoformswitch_part2_temp.rds")

#Analyze more SwitchConsequences:
source('./code/IsoformSwitchAnalyzeR/new_functions.R') #fixed TTS issues

aSwitchList_part2 <- analyzeSwitchConsequences_new(aSwitchList_part2, 
                                                   onlySigIsoforms = T, 
                                                   alpha = 0.05,
                                                   dIFcutoff = 0.05,
                                                   consequencesToAnalyze = c(
                                                     # Transcript
                                                     'tss',
                                                     'tts',
                                                     'last_exon',
                                                     'isoform_length',
                                                     'exon_number',
                                                     # 'intron_structure',
                                                     # 'intron_retention',
                                                     # 'isoform_class_code',
                                                     # cpat
                                                     #'coding_potential',
                                                     # ORF
                                                     # 'ORF_genomic',
                                                     'ORF_length',
                                                     '5_utr_length',
                                                     '3_utr_length',
                                                     # seq similarity
                                                     # 'isoform_seq_similarity',
                                                     # 'ORF_seq_similarity',
                                                     # '5_utr_seq_similarity',
                                                     # '3_utr_seq_similarity',
                                                     # ORF
                                                     'NMD_status',
                                                     # pfam
                                                     'domains_identified',
                                                     # 'genomic_domain_position',
                                                     'domain_length',
                                                     'domain_isotype',
                                                     
                                                     # SignalIP
                                                     'signal_peptide_identified',
                                                     
                                                     # IDR
                                                     'IDR_identified',
                                                     'IDR_length',
                                                     'IDR_type',
                                                     
                                                     # sub cell
                                                     'sub_cell_location',
                                                     'sub_cell_shift_to_cell_membrane',
                                                     'sub_cell_shift_to_cytoplasm',
                                                     'sub_cell_shift_to_nucleus',
                                                     'sub_cell_shift_to_Extracellular'
                                                     
                                                     # topology
                                                     #  'isoform_topology',
                                                     #  'extracellular_region_count',
                                                     #  'intracellular_region_count',
                                                     #  'extracellular_region_length',
                                                     #  'intracellular_region_length'
                                                   ))

#Export and save:
saveRDS(aSwitchList_part2, file = "./code/IsoformSwitchAnalyzeR/output/isoformswitch_part2_v4.rds")





