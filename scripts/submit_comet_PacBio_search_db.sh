#!/bin/bash
#SBATCH --job-name=submit_comet_PacBio_search_db
#SBATCH --output=slurm_logs/submit_comet_PacBio_search_db.out
#SBATCH --time=0-1:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

~/tools/comet.linux.exe -Pproc/comet/comet.params.high-low data/tc-1154/*.mzXML