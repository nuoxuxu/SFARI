#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --output=slurm_logs/salmon_index.out
#SBATCH --time=0-5:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

sbatch -J "slurm_logs/${1}.out" -o slurm_logs/"${1}".out -t 0-1:0 -N 1 -n 1 \
mamba activate patch_seq_spl; salmon quant -i proc/decoy_transcriptome -l A -1 data/illumina/SFARI_data/${1}_R1_001.fastq.gz -2 data/illumina/SFARI_data/${1}_R2_001.fastq.gz -p 30 --validateMappings -o proc/salmon_quants/${1}_quant
