#!/bin/bash

# List of read names
read_names=("NPC_2_2_transcript/4807745" "CN_2_1_transcript/5187203" "NPC_3_3_transcript/4128037" "CN_1_3_transcript/5264342" "iPSC_3_transcript/2707976" "NPC_3_1_transcript/3850976" "CN_2_2_transcript/4371872" "iPSC_1_transcript/2708362" "iPSC_2_transcript/3284227" "CN_3_2_transcript/5516681" "CN_1_2_transcript/5129410" "NPC_1_1_transcript/3824895" "NPC_1_3_transcript/4245305" "iPSC_3_transcript/2778395" "iPSC_1_transcript/2847360")

# Input BAM file
bam_file="proc/merged.bam"

# Loop through the read names
for read_name in "${read_names[@]}"; do
    echo "Searching for read: $read_name"
    # Search for the read in the BAM file using samtools view and the indexed BAM
    samtools view "$bam_file" | grep "$read_name" >> test.bam
done
