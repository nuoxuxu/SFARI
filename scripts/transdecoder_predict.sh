#!/usr/bin/env bash
#SBATCH --job-name=TransDecoder_predict
#SBATCH --output=slurm_logs/TransDecoder_predict.out
#SBATCH --time=0-1:0
#SBATCH -n 1
#SBATCH -N 1

module apptainer
apptainer exec -e -B="/scratch/s/shreejoy/nxu/SFARI:/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg TransDecoder.Predict --single_best_only -t proc/merge_collapsed.filtered.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6