#!/bin/bash
#SBATCH --job-name=STAR_NPC_2_2_S7
#SBATCH --output=slurm_logs/STAR_logs/NPC_2_2_S7_STAR.out
#SBATCH --error=slurm_logs/STAR_logs/NPC_2_2_S7_STAR.err
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 1

mamba activate SQANTI3.env

STAR --runThreadN 40 \
    --genomeDir /project/s/shreejoy/Genomic_references//GENCODE/STAR_index_v47 \
    --readFilesIn data/SFARI_data/NPC_2_2_S7_L001_R1_001.fastq.gz data/SFARI_data/NPC_2_2_S7_L001_R2_001.fastq.gz \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix STAR_results/NPC_2_2_S7_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts

