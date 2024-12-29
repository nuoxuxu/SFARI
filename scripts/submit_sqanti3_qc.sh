#!/bin/bash
#SBATCH --job-name=sqanti3_qc
#SBATCH --output=sqanti3_qc.out
#SBATCH --time=0-12:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

conda activate SQANTI3.env
python /home/s/shreejoy/nxu/tools/SQANTI3/sqanti3_qc.py \
-t 40 \
--skipORF \
--output SQANTI3_qc \
--CAGE_peak "${GENOMIC_DATA_DIR}"/Human_hg38_Gencode_v39/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
--polyA_motif_list "${GENOMIC_DATA_DIR}"/Human_hg38_Gencode_v39/polyA.list.txt \
proc/merged_collapsed.filtered.gff "${GENOMIC_DATA_DIR}"/Human_hg38_Gencode_v39/gencode.v39.annotation.sorted.gtf "${GENOMIC_DATA_DIR}"/GENCODE/GRCh38.primary_assembly.genome.fa