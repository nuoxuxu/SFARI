#!/bin/bash
#SBATCH --job-name=STAR_CN_2_2_S9
#SBATCH --output=slurm_logs/STAR_logs/CN_2_2_S9_STAR.out
#SBATCH --error=slurm_logs/STAR_logs/CN_2_2_S9_STAR.err
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 1

mamba activate SQANTI3.env

STAR --runThreadN 40 \
    --genomeDir /project/s/shreejoy/Genomic_references//GENCODE/STAR_index_v47 \
    --readFilesIn data/SFARI_data/CN_2_2_S9_L001_R1_001.fastq.gz data/SFARI_data/CN_2_2_S9_L001_R2_001.fastq.gz \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix STAR_results/CN_2_2_S9_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts

