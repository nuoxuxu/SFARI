#!/bin/bash
#SBATCH --job-name=CN_1_1_previous_tube_label_CN_1_2_S3_L001
#SBATCH --output=slurm_logs/CN_1_1_previous_tube_label_CN_1_2_S3_L001.out
#SBATCH --time=0-5:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
mamba activate patch_seq_spl
salmon quant -i proc/salmon_index -l A -1 data/illumina/SFARI_data/CN_1_1_previous_tube_label_CN_1_2_S3_L001_R1_001.fastq.gz -2 data/illumina/SFARI_data/CN_1_1_previous_tube_label_CN_1_2_S3_L001_R2_001.fastq.gz -p 30 --validateMappings -o proc/salmon_quants/CN_1_1_previous_tube_label_CN_1_2_S3_L001_quant
