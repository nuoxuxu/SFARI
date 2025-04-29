library(IsoformSwitchAnalyzeR)

aSwitchList <- readRDS("data/katherine/isoformswitch_inProg.rds")

aSwitchList_part2 <- analyzeIUPred2A(
  switchAnalyzeRlist       = aSwitchList,
  pathToIUPred2AresultFile = "export/iupred2a_processed_result.txt",
  showProgress = T)

aSwitchList <- analyzeSignalP(aSwitchList, pathToSignalPresultFile = "data/katherine/orfanage_peptide_summary.signalp5")

analyzeDeepLoc2(aSwitchList, pathToDeepLoc2resultFile = "export/results_20250325-202643.csv")

analyzeDeepTMHMM(
  switchAnalyzeRlist   = aSwitchList,
  pathToDeepTMHMMresultFile = "export/DeepTMHMM_output.txt",
  showProgress=FALSE
)
