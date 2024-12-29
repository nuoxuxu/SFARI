library(tidyverse)
library(IsoformSwitchAnalyzeR)

Isoseq_Expression <- read.csv("results/long_read/Isoseq_Expression.csv")

sampleID <- colnames(Isoseq_Expression)[c(-1)]
time_point <- str_split(sampleID, "_", 2) %>% map_chr(~ .x[1])
time_point <- factor(time_point, levels = c("iPSC", "NPC", "CN"))
myDesign <- data.frame(
    sampleID = sampleID,
    condition = time_point
)

IsoseqsSwitchList <- importRdata(
    isoformCountMatrix = Isoseq_Expression,
    designMatrix = myDesign,
    isoformExonAnnoation = "proc/isoformExonAnnoation.gtf",
    isoformNtFasta = "SQANTI3_qc_corrected.fasta",
    addAnnotatedORFs = FALSE,
    fixStringTieAnnotationProblem = FALSE
)

IsoseqsSwitchList <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = IsoseqsSwitchList,
    reduceToSwitchingGenes = FALSE,
    showProgress = TRUE
)

saveRDS(IsoseqsSwitchList, "results/long_read/IsoseqsSwitchList_DEXSeq.rds")