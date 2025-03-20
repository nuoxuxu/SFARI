#!/bin/bash
#SBATCH --job-name=pfam_scan
#SBATCH --output=pfam_scan.out
#SBATCH --time=2-0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=249G
#SBATCH --mail-user=nuoxu.xu@mail.utoronto.ca
#SBATCH --mail-type=FAIL,END

export PERL5LIB=/scratch/nxu/SFARI/PfamScan:$PERL5LIB
module load hmmer/3.4
pfam_scan/pfam_scan.py \
    -out full_pfam_scan_results.txt \
    -cpu 64 \
    orfanage_peptide.fasta \
    Pfam_flat_files