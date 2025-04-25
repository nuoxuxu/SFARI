#!/bin/bash
#SBATCH --job-name=iupred
#SBATCH --output=slurm_logs/iupred.out
#SBATCH --time=0-2:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

mamba activate patch_seq_spl
python /home/s/shreejoy/nxu/tools/iupred2a/iupred2a.py -a nextflow_results/V47/orfanage/orfanage_peptide.fasta long > export/iupred2a_result.txt