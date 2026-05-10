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
    module "StdEnv/2023:minimap2/2.28:samtools/1.22.1"
    label "short_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools/bambu"    
    input:
    path fastqz
    path genome_fasta
    output:
    path("${fastqz.simpleName}.aligned.bam")
    script:
    """
    minimap2 -ax splice:hq -t ${task.cpus} -uf ${genome_fasta} $fastqz | samtools sort -@ ${task.cpus} \\
        -o "${fastqz.simpleName}.aligned.bam"
    """
}

process index_bam {
    module "StdEnv/2023:samtools/1.22.1"
    label "short_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools/bambu"
    tag "${bam.simpleName}"

    input:
    path(bam)
    output:
    path("${bam}.bai")

    script:
    """
    samtools index -@ ${task.cpus} $bam
    """
}

process bambu {
    module "StdEnv/2023:gcc/12.3:r/4.5.0:r-bundle-bioconductor/3.21"
    label "long_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools/bambu"
    input:
    path minimap_bam
    path annotation_gtf
    path genome_fasta
    output:
    path("bambu_result.rds"), emit: bambu_rds

    script:
    """
    export R_LIBS="\${SCRATCH}/R/\${EBVERSIONR}"
    run_bambu.R ${annotation_gtf} "*.aligned.bam" ${genome_fasta} 96 bambu_result.rds
    """
}

process wrangle_bambu_result {
    conda "/scratch/nxu/astrocytes/env"
    label "short_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools/bambu"
    input:
    path(bambu_rds)
    output:
    path("supportedTranscriptModels.gtf"), emit: supported_tx_gtf
    path("supportedTxClassification.txt"), emit: supported_tx_classification
    path("bambu_expression.csv"),           emit: bambu_expression

    script:
    """
    wrangle_bambu_result.R ${bambu_rds}
    """
}

process plot_lrs_comparison {
    module "StdEnv/2023:gcc/12.3:arrow/19.0.1:rust/1.85.0:python/3.11.5"
    label "short_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools"

    input:
    path(bambu_gtf)
    path(bambu_classification)
    path(bambu_expression)
    path(isoseq_gtf)
    path(isoseq_classification)
    path(isoseq_expression)
    path(isoseq_filtered_lite_classification)
    path(star_sj_files)
    path(cage_bed)
    path(polya_site)

    output:
    path("lrs_comparison.pdf")
    path("lrs_comparison.png")

    script:
    """
    source ${projectDir}/.venv/bin/activate
    plot_lrs_comparison.py \\
        --bambu-gtf ${bambu_gtf} \\
        --bambu-classification ${bambu_classification} \\
        --bambu-expression ${bambu_expression} \\
        --isoseq-gtf ${isoseq_gtf} \\
        --isoseq-classification ${isoseq_classification} \\
        --isoseq-expression ${isoseq_expression} \\
        --isoseq-filtered-lite-classification ${isoseq_filtered_lite_classification} \\
        --sj ${star_sj_files.join(' ')} \\
        --cage ${cage_bed} \\
        --polya ${polya_site} \\
        --output lrs_comparison.pdf
    """
}

workflow COMPARE_OTHER_TRANSCRIPT_DISCOVERY {
    take:
    flnc_bam
    annotation_gtf
    genome_fasta
    isoseq_gtf
    isoseq_classification
    isoseq_expression
    isoseq_filtered_lite_classification
    star_sj_files
    cage_bed
    polya_site

    main:
    flnc_bam
        .map { path ->
            def sample = path.parent.parent.name
                            .replaceAll('_150PM_CELL1', '')
            [sample, path]
        }
        | convert_flnc_bam_to_fastqz
    minimap2_genome(convert_flnc_bam_to_fastqz.out, genome_fasta)
    index_bam(minimap2_genome.out)
    bambu(minimap2_genome.out.collect(), annotation_gtf, genome_fasta)
    wrangle_bambu_result(bambu.out.bambu_rds)
    plot_lrs_comparison(
        wrangle_bambu_result.out.supported_tx_gtf,
        wrangle_bambu_result.out.supported_tx_classification,
        wrangle_bambu_result.out.bambu_expression,
        isoseq_gtf,
        isoseq_classification,
        isoseq_expression,
        isoseq_filtered_lite_classification,
        star_sj_files,
        cage_bed,
        polya_site,
    )

    emit:
    genome_bam = minimap2_genome.out
}