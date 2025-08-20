#!/bin/bash
#SBATCH --job-name=merge_bam_NPC
#SBATCH --output=slurm_logs/merge_bam_NPC.out
#SBATCH --error=slurm_logs/merge_bam_NPC.err
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 1
mamba activate /gpfs/fs0/scratch/s/shreejoy/nxu/SFARI/work/conda/env-0930b165b530f85f-a5e01a107d0bbdbd93d5689448b78b9f
samtools merge data/long_read/LUO26876.20240514/NPC*150PM_CELL1/outputs/mapped.bam -o data/long_read/NPC_merged.bam -@ 40 -f