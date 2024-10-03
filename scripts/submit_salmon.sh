#!/bin/bash
#SBATCH --job-name=salmon_illumina
#SBATCH --output=salmon_illumina.out
#SBATCH --time=0-5:0
#SBATCH --nodes=1
#SBATCH --ntasks=1

conda activate SQANTI3.env
salmon 
