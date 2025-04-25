#!/bin/bash
#SBATCH --job-name=sqanti3_qc
#SBATCH --output=slurm_logs/sqanti3_qc.out
#SBATCH --time=0-6:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

mamba activate SQANTI3.env
python /home/s/shreejoy/nxu/tools/SQANTI3/sqanti3_qc.py \
    -t 40 \
    --skipORF \
    --output SQANTI3_qc \
    --CAGE_peak "${GENOMIC_DATA_DIR}"/Human_hg38_Gencode_v39/refTSS_v3.3_human_coordinate.hg38.sorted.bed \
    --polyA_motif_list "${GENOMIC_DATA_DIR}"/Human_hg38_Gencode_v39/polyA.list.txt \
    "nextflow_results/V47/merged_collapsed.sorted.filtered_lite.gff" "${GENOMIC_DATA_DIR}/GENCODE/gencode.v47.annotation.sorted.gtf" "${GENOMIC_DATA_DIR}/GENCODE/GRCh38.primary_assembly.genome.fa"