process convert_flnc_bam_to_fastqz {
    module "StdEnv/2023:samtools/1.22.1"
    label "short_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools/bambu"
    tag "${sample_name}"

    input:
    tuple val(sample_name), path(flnc_bam)
    output:
    path("${sample_name}.fastq.gz")

    script:
    """
    samtools fastq -@ ${task.cpus} -0 ${sample_name}.fastq.gz \\
    -c 9 $flnc_bam
    """
}

process minimap2_genome {
    module "StdEnv/2023:minimap2/2.28"
    label "short_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools/bambu"    
    input:
    path fastqz
    path ref_genome_fasta
    output:
    path("${fastqz.simpleName}.aligned.bam")
    script:
    """
    minimap2 -ax splice:hq -t ${task.cpus} -uf ${params.ref_genome_fasta} $fastqz | samtools sort -@ ${task.cpus} \\
        -o "${fastqz.simpleName}.aligned.bam"
    """
}

process bambu {
    conda "StdEnv/2023:gcc/12.3:r/4.5.0:r-bundle-bioconductor/3.21"
    label "mid_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools/bambu"
    input:
    path minimap_bam
    path annotation_gtf
    path ref_genome_fasta
    output:
    path("supportedTranscriptModels.gtf"), emit: supported_tx_gtf
    path("novelTranscripts.gtf"), emit: novel_tx_gtf
    script:
    """
    export R_LIBS="\${SCRATCH}/R/\${EBVERSIONR}"
    run_bambu.R ${params.annotation_gtf} "*minimap.bam" ${params.ref_genome_fasta} ${task.cpus}
    """
}

workflow BAMBU {
    take:
    flnc_bam

    main:
    // Can't run this on Nibi, not enough space
    flnc_bam
        .map { path ->
            def sample = path.parent.parent.name  // gets "CN_1_2_150PM_CELL1"
                            .replaceAll('_150PM_CELL1', '')  // gives "CN_1_2"
            [sample, path]
        }
        | convert_flnc_bam_to_fastqz
    minimap2_genome(convert_flnc_bam_to_fastqz.out)
    bambu(minimap2_genome.out.collect())

    emit:
    supported_tx_gtf = bambu.out.supported_tx_gtf
    novel_tx_gtf = bambu.out.novel_tx_gtf
}