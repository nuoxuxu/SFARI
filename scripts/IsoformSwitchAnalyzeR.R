library(IsoformSwitchAnalyzeR)

aSwitchList <- readRDS("data/katherine/isoformswitch_inProg.rds")

aSwitchList_part2 <- analyzeIUPred2A(
  switchAnalyzeRlist       = aSwitchList,
  pathToIUPred2AresultFile = "export/iupred2a_processed_result.txt",
  showProgress = T)