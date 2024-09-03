#!/bin/bash
#SBATCH --job-name=sqanti3_qc
#SBATCH --output=sqanti3_qc.out
#SBATCH -t 1-0:0
#SBATCH --nodes=1
conda activate SQANTI3.env
python /home/s/shreejoy/nxu/tools/SQANTI3-5.2.2/sqanti3_qc.py --skipORF proc/merged_collapsed.filtered.gff "${GENOMIC_DATA_DIR}"/GENCODE/gencode.v33.annotation.sorted.gtf "${GENOMIC_DATA_DIR}"/GENCODE/GRCh38.primary_assembly.genome.fa