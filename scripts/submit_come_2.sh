#!/bin/bash
#SBATCH --job-name=submit_comet
#SBATCH --output=slurm_logs/submit_comet_2.out
#SBATCH --time=0-1:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

~/tools/comet.linux.exe -Pproc/comet/comet.params.transdecoder data/tc-1154/*.mzXML