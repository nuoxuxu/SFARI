#!/bin/bash
#SBATCH --job-name=salmon_generateDecoy
#SBATCH --output=slurm_logs/salmon_generateDecoy.out
#SBATCH --time=0-2:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

mamba activate SQANTI3.env
./scripts/generateDecoyTranscriptome.sh \
-j 40 \
-g "${GENOMIC_DATA_DIR}"/GENCODE/GRCh38.primary_assembly.genome.fa \
-t proc/merged_collapsed_filtered.fasta \
-a proc/merged_collapsed.filtered.gff \
-o proc/decoy_transcriptome