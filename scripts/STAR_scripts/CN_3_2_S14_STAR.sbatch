#!/bin/bash
#SBATCH --job-name=STAR_CN_3_2_S14
#SBATCH --output=slurm_logs/STAR_logs/CN_3_2_S14_STAR.out
#SBATCH --error=slurm_logs/STAR_logs/CN_3_2_S14_STAR.err
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 1

mamba activate SQANTI3.env

STAR --runThreadN 40 \
    --genomeDir /project/s/shreejoy/Genomic_references//GENCODE/STAR_index_v47 \
    --readFilesIn data/SFARI_data/CN_3_2_S14_L001_R1_001.fastq.gz data/SFARI_data/CN_3_2_S14_L001_R2_001.fastq.gz \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix STAR_results/CN_3_2_S14_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts

