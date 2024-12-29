#!/usr/bin/env bash
#SBATCH --job-name=hmmsearch_predict
#SBATCH --output=slurm_logs/hmmsearch_predict.out
#SBATCH --time=0-12:0
#SBATCH -n 1
#SBATCH -N 1

module apptainer
apptainer exec -e -B="/scratch/s/shreejoy/nxu/SFARI:/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg hmmsearch --cpu 40 -E 1e-10 --domtblout pfam.domtblout Pfam-A.hmm full_nt.fasta.transdecoder_dir/longest_orfs.pep
apptainer exec -e -B="/scratch/s/shreejoy/nxu/SFARI:/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg TransDecoder.Predict --single_best_only -t full_nt.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6