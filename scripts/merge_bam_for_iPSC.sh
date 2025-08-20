#!/bin/bash
#SBATCH --job-name=merge_bam_iPSC
#SBATCH --output=slurm_logs/merge_bam_iPSC.out
#SBATCH --error=slurm_logs/merge_bam_iPSC.err
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 1
module load samtools/1.13
samtools merge data/long_read/LUO26876.20240514/iPSC*150PM_CELL1/outputs/mapped.bam -o data/long_read/iPSC_merged.bam -@ 40 -f