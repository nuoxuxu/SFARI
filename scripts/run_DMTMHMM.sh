#!/bin/sh
mamba deactivate
module load NiaEnv/2019b python/3.11.5
source .virtualenvs/DeepTMHMM/bin/activate

DIRECTORY="DeepTMHMM_input"
i=0

for file in "$DIRECTORY"/*
do
  i=$((i+1))  # Increment the counter
  echo "File #$i: $file"
  if [ "$i" -gt 6 ]; then
    if [ "$(hostname)" == "xunuos-MacBook-Air.local" ]; then
      biolib run --local 'DTU/DeepTMHMM:1.0.24' --fasta "${file}"
      mv "biolib_results" "export/DeepTMHMM_outputs/${i}_output"
    else
      biolib run DTU/DeepTMHMM --fasta "${file}"
      mv "biolib_results" "export/DeepTMHMM_outputs/${i}_output"
    fi
  else
    echo "Skipping $file (i <= 6)"
  fi    
done