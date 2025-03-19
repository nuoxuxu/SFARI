#!/bin/bash
#SBATCH --job-name=STAR_NPC_1_2_previous_tube_label_NPC_1_3_S2
#SBATCH --output=slurm_logs/STAR_logs/NPC_1_2_previous_tube_label_NPC_1_3_S2_STAR.out
#SBATCH --error=slurm_logs/STAR_logs/NPC_1_2_previous_tube_label_NPC_1_3_S2_STAR.err
#SBATCH --time=02:00:00
#SBATCH -N 1
#SBATCH -n 1

mamba activate SQANTI3.env

STAR --runThreadN 40 \
    --genomeDir /project/s/shreejoy/Genomic_references//GENCODE/STAR_index_v47 \
    --readFilesIn data/SFARI_data/NPC_1_2_previous_tube_label_NPC_1_3_S2_L001_R1_001.fastq.gz data/SFARI_data/NPC_1_2_previous_tube_label_NPC_1_3_S2_L001_R2_001.fastq.gz \
    --readFilesCommand gunzip -c \
    --outFileNamePrefix STAR_results/NPC_1_2_previous_tube_label_NPC_1_3_S2_ \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts

