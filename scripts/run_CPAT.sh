#!/bin/sh

mamba activate SQANTI3.env
cpat \
  -x "data/Human_Hexamer.tsv" \
  -d "data/Human_logitModel.RData" \
  -g "full_nt.fasta" \
  --min-orf=50 \
  --top-orf=50 \
  -o SFARI \
  1> SFARI_cpat.output \
  2> SFARI_cpat.error

scripts/orf_calling.py \
  --orf_coord "SFARI.ORF_prob.tsv" \
  --orf_fasta "SFARI.ORF_seqs.fa" \
  --gencode "/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.annotation.gtf" \
  --sample_gtf "proc/isoformExonAnnoation.gtf" \
  --pb_gene "pb_gene.txt" \
  --classification "merged_collapsed_classification.filtered_lite_classification.txt" \
  --sample_fasta "full_nt.fasta" \
  --num_cores 1 \
  --output SFARI_best_orf.tsv