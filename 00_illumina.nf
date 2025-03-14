process generateSalmonIndex {

    label "short_slurm_job"

    input:
    path final_sample_fasta

    output:
    path "salmon_index"
    script:
    """
    salmon index -p $task.cpus -t $final_sample_fasta -i salmon_index
    """
}

process generateDecoyTranscriptome {
    
    label "short_slurm_job"

    input:
    path genome_fasta
    path final_sample_gtf
    path final_sample_fasta

    output:
    path "decoy_transcriptome/"
    
    script:
    """
    generateDecoyTranscriptome.sh \\
        -j $task.cpus \\
        -g $genome_fasta \\
        -t $final_sample_fasta \\
        -a $final_sample_gtf \\
        -o decoy_transcriptome
    """
}
process runSalmon {
    input:
    path salmon_index
    path fastq_gz
    
    output:
    path "salmon_quants/${fastq_gz.baseName}_quant"

    script:
    """
    salmon quant -i proc/decoy_transcriptome -l A -1 data/illumina/SFARI_data/${1}_R1_001.fastq.gz -2 data/illumina/SFARI_data/${1}_R2_001.fastq.gz -p 30 --validateMappings -o proc/salmon_quants/${1}_quant    
    """
}

workflow salmon {
    take:

    main:
    generateSalmonIndex(params.transcriptome_fasta)
}