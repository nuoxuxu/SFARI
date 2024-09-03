#!/bin/bash
for file in data/long_read/LUO26876.20240514/*/outputs/flnc.bam; do
    if [[ -f "$file" ]]; then
        first_field=$(samtools view "$file" | head -1 | cut -f1)
        if [[ -n "$first_field" ]]; then
            echo -e "$first_field\t$file" >> proc/id_to_sample.txt
        else
            echo -e "No data\t$file" >> proc/id_to_sample.txt
        fi
    else
        echo -e "File not found\t$file" >> proc/id_to_sample.txt
    fi
done