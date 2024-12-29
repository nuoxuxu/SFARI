library(IsoformSwitchAnalyzeR)

IsoseqsSwitchList <- readRDS("results/long_read/IsoseqsSwitchList.rds")

IsoseqsSwitchList <- analyzeAlternativeSplicing(IsoseqsSwitchList, onlySwitchingGenes=FALSE)

saveRDS(IsoseqsSwitchList, "results/long_read/full_analyzeAlternativeSplicing.rds")