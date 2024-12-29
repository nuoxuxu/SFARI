#!/usr/bin/env bash
#SBATCH --job-name=hmmsearch
#SBATCH --output=slurm_logs/hmmsearch.out
#SBATCH --time=0-12:0
#SBATCH -n 1
#SBATCH -N 1

module load apptainer
apptainer exec -e -B="/scratch/s/shreejoy/nxu/SFARI:/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg hmmsearch --cpu 40 -E 1e-10 --domtblout pfam.domtblout Pfam-A.hmm full_nt.fasta.transdecoder_dir/longest_orfs.pep