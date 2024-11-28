#!/usr/bin/env nextflow

// listOfFiles = file('data/long_read/LUO26876.20240514/*/outputs/flnc.bam')

// println listOfFiles

params.outdir = "/scratch/s/shreejoy/nxu/SFARI/data/long_read/LUO26876.20240514"

process getIDToSample {
    publishDir "proc", mode: 'copy'
    
    input:
    val x

    output:
    path "id_to_sample.txt"

    script:
    """
    for file in $x; do
        if [[ -f "\$file" ]]; then
            first_field=\$(~/miniforge3/envs/SQANTI3.env/bin/samtools view "\$file" | head -1 | cut -f1)
            if [[ -n "\$first_field" ]]; then
                echo -e "\$first_field\t\$file" >> id_to_sample.txt
            else
                echo -e "No data\t\$file" >> id_to_sample.txt
            fi
        else
            echo -e "File not found\t\$file" >> id_to_sample.txt
        fi
    done        
    """
}

process mergeBamFiles {
    publishDir "proc", mode: 'copy'
    
    input:
    path "bam"

    output:
    path "bam_files.txt"

    script:
    """
    ls bam* > bam_files.txt
    """
}

workflow {
    getIDToSample(params.outdir + "/*/outputs/flnc.bam")
    Channel.fromPath(params.outdir + "/*/outputs/mapped.bam").collect().set { bamFiles }
    mergeBamFiles(bamFiles)
}