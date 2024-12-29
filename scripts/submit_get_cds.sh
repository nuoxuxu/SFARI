#!/usr/bin/env bash
#SBATCH --job-name=get_cds
#SBATCH --output=slurm_logs/get_cds.out
#SBATCH --time=0-12:0
#SBATCH -n 1
#SBATCH -N 1

mamba activate patch_seq_spl
python scripts/get_cds.py