#!/usr/bin/env bash
#SBATCH --job-name=cdna_alignment_orf_to_genome_orf
#SBATCH --output=cdna_alignment_orf_to_genome_orf.out
#SBATCH --time=0-2:0
#SBATCH -n 1
#SBATCH -N 1

# module apptainer
# apptainer exec -e -B="/scratch/s/shreejoy/nxu/SFARI:/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg /usr/local/bin/util/cdna_alignment_orf_to_genome_orf.pl \
#     merge_collapsed.filtered.fasta.transdecoder.gff3 \
#     proc/merged_collapsed.sorted.filtered_lite.gff3 \
#     proc/merge_collapsed.filtered.fasta > proc/merge_collapsed.filtered.fasta.transdecoder.genome.gff3

apptainer exec -e -B="${PWD}:${PWD},/gpfs/fs0/scratch/s/shreejoy/nxu/SFARI:/gpfs/fs0/scratch/s/shreejoy/nxu/SFARI" ~/tools/transdecoder.v5.7.1.simg /usr/local/bin/util/cdna_alignment_orf_to_genome_orf.pl \
    merge_collapsed.filtered.fasta.transdecoder.gff3 \
    transcripts.gff3 \
    transcripts.fasta > merge_collapsed.filtered.fasta.transdecoder.genome.gff3