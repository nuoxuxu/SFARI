#!/bin/bash
#SBATCH --job-name=deeploc
#SBATCH --output=slurm_logs/deeploc.out
#SBATCH --time=1-0:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

mamba deactivate
module load NiaEnv/2019b python/3.11.5
source .virtualenvs/deeploc/bin/activate

.virtualenvs/deeploc/bin/deeploc2 \
    -f /scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/orfanage/orfanage_peptide.fasta \
    -o export/