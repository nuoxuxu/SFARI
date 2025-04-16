#!/usr/bin/env bash

# Usage check
if [ $# -ne 2 ]; then
  echo "Usage: $0 <fasta_file> <entries_per_file>"
  exit 1
fi

# Input arguments
input_file="$1"
x="$2"

# Remove trailing .fasta if it exists to get a base name
# (If your file doesn’t end with .fasta, basename without the second argument 
# will just remove the path, which is still fine.)
base=$(basename "$input_file" .fasta)

# Make a directory named after the base
mkdir -p "$base"

# Use awk to do the splitting
awk -v x="$x" -v base="$base" '
BEGIN {
  filecount = 1
  entrycount = 0
  # Construct the first output file name
  outfilename = base "/" base "_" filecount ".fasta"
}

/^>/ {
  # We encountered a new FASTA entry
  if (entrycount == x) {
    # If weve already written x entries, move to a new file
    close(outfilename)
    filecount++
    entrycount = 0
    outfilename = base "/" base "_" filecount ".fasta"
  }
  entrycount++
}

# Append the current line to the current output file
{ 
  print $0 >> outfilename
}
' "$input_file"

echo "Splitting complete. Check the '$base' directory."
