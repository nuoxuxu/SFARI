#!/usr/bin/env nextflow

params.datadir = "/scratch/s/shreejoy/nxu/SFARI/data/"
params.output_dir = "nextflow_results/"

process getIDToSample {
    publishDir "${params.output_dir}", mode: 'copy'
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
                echo -e "\$first_field\\t\$file" >> id_to_sample.txt
            else
                echo -e "No data\\t\$file" >> id_to_sample.txt
            fi
        else
            echo -e "File not found\\t\$file" >> id_to_sample.txt
        fi
    done
    """
}

process mergeBamFiles {

    label "short_slurm_job"
    
    input:
    path "bam"

    output:
    path "merged.bam"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/samtools merge bam* -o merged.bam -@ $task.cpus
    """
}

process isoseqCollapse {
    publishDir "${params.output_dir}", mode: 'copy'
    label "short_slurm_job"
    
    input:
    path merged_bam

    output:
    path "merged_collapsed.abundance.txt"
    path "merged_collapsed.fasta", emit: sample_fasta
    path "merged_collapsed.flnc_count.txt"
    path "merged_collapsed.gff", emit: isoform_gff
    path "merged_collapsed.group.txt"
    path "merged_collapsed.read_stat.txt", emit: read_stat
    path "merged_collapsed.report.json"

    script:
    """
    ~/miniforge3/envs/isoseq/bin/isoseq collapse -j $task.cpus $merged_bam merged_collapsed.gff
    """
}

workflow merge_and_collapse {
    take:
    flnc_bam
    mapped_bam
    
    main:
    getIDToSample(params.datadir + flnc_bam)
    Channel.fromPath(params.datadir + mapped_bam).collect().set { bamFiles }
    mergeBamFiles(bamFiles)
    isoseqCollapse(mergeBamFiles.out)    
}

workflow {
    flnc_bam = "long_read/LUO26876.20240514/*/outputs/flnc.bam"
    mapped_bam = "long_read/LUO26876.20240514/*/outputs/mapped.bam"
    merge_and_collapse(flnc_bam, mapped_bam)
}