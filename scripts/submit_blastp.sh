#!/usr/bin/env bash
#SBATCH --job-name=blastp
#SBATCH --output=slurm_logs/blastp.out
#SBATCH --time=0-4:0
#SBATCH -n 1
#SBATCH -N 1

module apptainer
apptainer exec -e -B="/scratch/s/shreejoy/nxu/SFARI:/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg blastp -query full_nt.fasta.transdecoder_dir/longest_orfs.pep -db uniprotkb_proteome_UP000005640_AND_revi_2024_10_07.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 40 > blastp.outfmt6