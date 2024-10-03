#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --output=slurm_logs/salmon_index.out
#SBATCH --time=0-1:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

mamba activate patch_seq_spl
salmon index -p 30 -t proc/merged_collapsed.fasta -i proc/salmon_index_full