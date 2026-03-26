#!/bin/bash

# 1. Detect total available CPUs (using nproc)
TOTAL_CPUS=$(nproc)

# 2. Get the list of unique primer IDs
# We capture this into a variable to count lines and pass to parallel
PRIMER_IDS=$(ls *.bam | cut -d '.' -f 2 | sort | uniq)

# 3. Count how many unique primer groups there are
NUM_GROUPS=$(echo "$PRIMER_IDS" | wc -l)

# 4. Calculate threads per job
# Math: Total CPUs / Number of Groups
if [ "$NUM_GROUPS" -gt 0 ]; then
    THREADS_PER_JOB=$(( TOTAL_CPUS / NUM_GROUPS ))
else
    echo "No BAM files found."
    exit 1
fi

# 5. Logic adjustment for Samtools flags
# Samtools -@ flag specifies *additional* threads. 
# If we calculated 1 thread total, we want -@ 0. If 4 total, we want -@ 3.
if [ "$THREADS_PER_JOB" -gt 1 ]; then
    SAMTOOLS_THREADS=$(( THREADS_PER_JOB - 1 ))
else
    # If the division results in 0 (more groups than CPUs), use minimal resources
    SAMTOOLS_THREADS=0
fi

echo "------------------------------------------------"
echo "Detected $TOTAL_CPUS CPUs and $NUM_GROUPS primer groups."
echo "Running $NUM_GROUPS jobs in parallel."
echo "Allocating $SAMTOOLS_THREADS *additional* threads per samtools job."
echo "------------------------------------------------"

# 6. Run with GNU Parallel
# -j $NUM_GROUPS : Forces parallel to run ALL groups at once
# {}             : The placeholder for the primer ID being processed
echo "$PRIMER_IDS" | parallel -j "$NUM_GROUPS" \
    "samtools merge -@ $SAMTOOLS_THREADS -o {}.flnc.bam *.{}.flnc.bam"