#!/bin/bash
#SBATCH --job-name=collapse_bam
#SBATCH --output=slurm_logs/collapse_bam.out
#SBATCH --time=10:00:00
#SBATCH -N 1
#SBATCH -n 1

bam_dir="/scratch/s/shreejoy/ehogan/sfari/riboseq/intermediate/aligned_BAM"
mamba activate patch_seq_spl
# samtools merge ${bam_dir}/merged.bam ${bam_dir}/*Aligned.sortedByCoord.out.bam

bedtools bamtobed -i ${bam_dir}/merged.bam \
  | sort -k1,1 -k2,2n -k3,3n -k6,6 \
  | bedtools groupby -i - -g 1,2,3,6 -c 1 -o count \
  > ${bam_dir}/slurm_collapsed.bed