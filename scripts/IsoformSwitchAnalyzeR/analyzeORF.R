library("BSgenome.Hsapiens.UCSC.hg38")
library(IsoformSwitchAnalyzeR)
bsg <- BSgenome.Hsapiens.UCSC.hg38

IsoseqsSwitchList_part1 <- readRDS("proc/IsoseqsSwitchList_part1.rds")
IsoseqsSwitchList_part1_Analyzed <- analyzeORF(IsoseqsSwitchList_part1, bsg)
saveRDS(IsoseqsSwitchList_part1_Analyzed, "IsoseqsSwitchList_part1_Analyzed.rds")