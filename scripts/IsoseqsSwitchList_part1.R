library(IsoformSwitchAnalyzeR)

IsoseqsSwitchList <- readRDS("IsoseqsSwitchList.rds")
IsoseqsSwitchList <- preFilter(
switchAnalyzeRlist         = IsoseqsSwitchList,
geneExpressionCutoff       = 1, # default
isoformExpressionCutoff    = 0, # default
IFcutoff                   = 0.01, # default
removeSingleIsoformGenes   = TRUE, # default
reduceToSwitchingGenes     = FALSE, # default (we didn"t run DEXSeq yet)
keepIsoformInAllConditions = TRUE # we only have 2 conditions so doesn"t matter
)

IsoseqsSwitchList_part1 <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = IsoseqsSwitchList,
    reduceToSwitchingGenes = FALSE
)
saveRDS(IsoseqsSwitchList_part1, "proc/IsoseqsSwitchList_part1.rds")