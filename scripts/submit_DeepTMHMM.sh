#!/bin/bash
#SBATCH --job-name=DeepTMHMM
#SBATCH --output=slurm_logs/DeepTMHMM.out
#SBATCH --time=0-5:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

mamba deactivate
module load NiaEnv/2019b python/3.11.5
source .virtualenvs/DeepTMHMM/bin/activate

# This does not work in compute nodes because of the lack of internet access
# biolib run DTU/DeepTMHMM --fasta /scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/orfanage/orfanage_peptide.fasta

biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta /scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/orfanage/orfanage_peptide.fasta