library(IsoformSwitchAnalyzeR)
library("BSgenome.Hsapiens.UCSC.hg38")
bsg <- BSgenome.Hsapiens.UCSC.hg38

IsoseqsSwitchList <- readRDS("results/long_read/IsoseqsSwitchList.rds")

IsoseqsSwitchList <- extractSequence(
    IsoseqsSwitchList,
    bsg,
    onlySwitchingGenes = FALSE,
    removeShortAAseq = FALSE,
    removeORFwithStop = FALSE
)

IsoseqsSwitchList <- analyzeSwitchConsequences(
    IsoseqsSwitchList,
    consequencesToAnalyze = c('coding_potential','NMD_status','domains_identified','ORF_seq_similarity'),
    dIFcutoff = 0.1, # very high cutoff for fast runtimes - you should use the default (0.1)
    showProgress=FALSE
)

saveRDS(IsoseqsSwitchList, "results/long_read/IsoseqsSwitchList_final.rds")