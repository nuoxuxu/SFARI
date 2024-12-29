#!/usr/bin/env bash
#SBATCH --job-name=TransDecoder_LongOrfs
#SBATCH --output=slurm_logs/TransDecoder_LongOrfs.out
#SBATCH --time=0-1:0
#SBATCH -n 1
#SBATCH -N 1

module apptainer
apptainer exec -e -B="/scratch/s/shreejoy/nxu/SFARI:/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg TransDecoder.LongOrfs -S -t proc/merge_collapsed.filtered.fasta