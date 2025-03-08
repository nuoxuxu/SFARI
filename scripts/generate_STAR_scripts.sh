#!/bin/bash

# Define directories
input_dir="data/SFARI_data"  # Directory containing FASTQ files
genome_dir="${GENOMIC_DATA_DIR}/GENCODE/STAR_index_v47"  # STAR genome index directory
output_dir="STAR_results"  # Output directory for STAR results
sbatch_dir="scripts/STAR_scripts"  # Directory to store generated sbatch scripts
log_dir="slurm_logs/STAR_logs"  # Log directory

# Ensure output and sbatch directories exist
mkdir -p "$output_dir" "$sbatch_dir" "$log_dir"

# SLURM settings
time="02:00:00"

# Loop over R1 FASTQ files and generate sbatch scripts
for r1 in "$input_dir"/*R1_001.fastq.gz; do
    # Extract sample name
    sample=$(basename "$r1" | sed 's/_L001_R1_001.fastq.gz//')
    r2="${input_dir}/${sample}_L001_R2_001.fastq.gz"

    # Check if R2 file exists
    if [[ ! -f "$r2" ]]; then
        echo "Paired file missing for sample: $sample. Skipping..."
        continue
    fi

    # Create sbatch script
    sbatch_file="${sbatch_dir}/${sample}_STAR.sbatch"
    cat <<EOF > "$sbatch_file"
#!/bin/bash
#SBATCH --job-name=STAR_${sample}
#SBATCH --output=${log_dir}/${sample}_STAR.out
#SBATCH --error=${log_dir}/${sample}_STAR.err
#SBATCH --time=$time
#SBATCH -N 1
#SBATCH -n 1

mamba activate SQANTI3.env

STAR --runThreadN 40 \\
    --genomeDir $genome_dir \\
    --readFilesIn $r1 $r2 \\
    --readFilesCommand gunzip -c \\
    --outFileNamePrefix ${output_dir}/${sample}_ \\
    --outSAMtype BAM SortedByCoordinate \\
    --quantMode GeneCounts

EOF

    echo "Generated: $sbatch_file"
done

echo "All sbatch scripts generated in $sbatch_dir."
