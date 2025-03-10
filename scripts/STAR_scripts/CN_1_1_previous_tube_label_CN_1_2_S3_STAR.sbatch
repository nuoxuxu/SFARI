#!/bin/bash
#SBATCH --job-name=STAR_CN_1_1_previous_tube_label_CN_1_2_S3
#SBATCH --output=slurm_logs/STAR_logs/CN_1_1_previous_tube_label_CN_1_2_S3_STAR.out
#SBATCH --error=slurm_logs/STAR_logs/CN_1_1_previous_tube_label_CN_1_2_S3_STAR.err
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 1

mamba activate SQANTI3.env

STAR --runThreadN 40 \
    --genomeDir /project/s/shreejoy/Genomic_references//GENCODE/STAR_index_v47 \
    --readFilesIn data/SFARI_data/CN_1_1_previous_tube_label_CN_1_2_S3_L001_R1_001.fastq.gz data/SFARI_data/CN_1_1_previous_tube_label_CN_1_2_S3_L001_R2_001.fastq.gz \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix STAR_results/CN_1_1_previous_tube_label_CN_1_2_S3_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts

