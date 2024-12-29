#!/usr/bin/env bash
#SBATCH --job-name=gtf_to_alignment_gff3
#SBATCH --output=slurm_logs/gtf_to_alignment_gff3.out
#SBATCH --time=0-2:0
#SBATCH -n 1
#SBATCH -N 1

module apptainer
apptainer exec -e -B="/scratch/s/shreejoy/nxu/SFARI:/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg /usr/local/bin/util/gtf_to_alignment_gff3.pl proc/merged_collapsed.sorted.filtered_lite.gff > proc/merged_collapsed.sorted.filtered_lite.gff3

