#!/bin/sh
if [ "$(hostname)" == "xunuos-MacBook-Air.local" ]; then
  mamba activate deeptmhmm
else
  mamba deactivate
  module load NiaEnv/2019b python/3.11.5
  source .virtualenvs/DeepTMHMM/bin/activate
fi

DIRECTORY="assets/DeepTMHMM_input"

for file in "$DIRECTORY"/*
do
  i="${file##*_}"
  echo "File #$i: $file"
  if [ -d "export/DeepTMHMM_outputs/${i}_output" ]; then
    echo "Skipping $file (alread processed)"
  else
    biolib run DTU/DeepTMHMM --fasta "${file}"
    mv "biolib_results" "export/DeepTMHMM_outputs/${i}_output"
  fi    
done