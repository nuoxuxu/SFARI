#!/usr/bin/env bash
#SBATCH -J getSingleCellObject
#SBATCH -o slurm_logs/getSingleCellObject.out
#SBATCH -t 0-1:0
#SBATCH -n 1
#SBATCH -N 1

bin/get_SingleCell_object.py \
    --id_to_sample proc/id_to_sample.txt \
    --classification proc/merged_collapsed_classification.filtered_lite_classification.txt \
    --read_stat proc/merged_collapsed.read_stat.txt \
    --output proc/pbid.h5ad