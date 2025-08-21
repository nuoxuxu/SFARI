#!/usr/bin/env nextflow

process getIDToSample {
    publishDir "${params.output_dir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    
    input:
    val x

    output:
    path "id_to_sample.txt"

    script:
    """
    for file in $x; do
        if [[ -f "\$file" ]]; then
            first_field=\$(samtools view "\$file" | head -1 | cut -f1)
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

    conda "${moduleDir}/environment.yml"
    
    input:
    path "bam"

    output:
    path "merged.bam"

    script:
    """ 
    samtools merge bam* -o merged.bam -@ $task.cpus
    """
}

process isoseqCollapse {
    publishDir "${params.output_dir}", mode: 'copy'
    label "short_slurm_job"

    conda "${moduleDir}/environment.yml"
    
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
    isoseq collapse -j $task.cpus $merged_bam merged_collapsed.gff
    """
}

workflow merge_and_collapse {
    take:
    flnc_bam
    mapped_bam

    main:
    getIDToSample(flnc_bam)
    Channel.fromPath(mapped_bam).collect().set { bamFiles }
    mergeBamFiles(bamFiles)
    isoseqCollapse(mergeBamFiles.out)

    emit:
    isoform_gff = isoseqCollapse.out.isoform_gff
    id_to_sample = getIDToSample.out
    read_stat = isoseqCollapse.out.read_stat
}

workflow {
    getIDToSample(Channel.fromPath(params.flnc_bam))
    Channel.fromPath(params.mapped_bam).collect().set { bamFiles }
    mergeBamFiles(bamFiles)
    isoseqCollapse(mergeBamFiles.out)
}