library(IsoformSwitchAnalyzeR)
library("BSgenome.Hsapiens.UCSC.hg38")

IsoseqsSwitchList <- readRDS("results/long_read/full_analyzeAlternativeSplicing.rds")

bsg <- BSgenome.Hsapiens.UCSC.hg38

IsoseqsSwitchList <- addORFfromGTF(
    IsoseqsSwitchList,
    overwriteExistingORF=TRUE,
    "addORFfromGTF.gtf")

IsoseqsSwitchList <- analyzeNovelIsoformORF(
    IsoseqsSwitchList,
    TRUE,
    bsg,
    minORFlength = 0
)

message("analyzeNovelIsoformORF done")

extractSequence(IsoseqsSwitchList, bsg, outputPrefix = "full", onlySwitchingGenes = FALSE)

message("extractSequence done")

saveRDS(IsoseqsSwitchList, "results/long_read/full_analyzeAlternativeSplicing.rds")

# unique(IsoseqsSwitchList$isoformFeatures$isoform_id) %>%
#     writeLines("full_pbid_list.txt")

# system("seqkit grep -f full_pbid_list.txt SQANTI3_qc_corrected.fasta > full_cpc2.fasta")

# system('sbatch -o ${SCRATCH}/SFARI/slurm_logs/full_cpc2.out -J full_cpc2 -N 1 -n 1 -t 0-12:0 --wrap="python ./bin/CPC2.py -i /scratch/s/shreejoy/nxu/SFARI/full_cpc2.fasta -o /scratch/s/shreejoy/nxu/SFARI/full_CPC2result"')

# system('sbatch scripts/pfam_scan.sh')