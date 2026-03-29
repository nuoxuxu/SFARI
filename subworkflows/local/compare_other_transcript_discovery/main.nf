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

process isoquant {
    module "StdEnv/2023:python/3.11.5:gcc/12.3:arrow/19.0.1:rust/1.85.0:minimap2/2.28:samtools/1.22.1"
    label "long_slurm_job"
    storeDir "nextflow_results/comapre_other_LRS_tools/isoquant"

    input:
    path(minimap_bam)
    path(bam_indexes)
    path(annotation_gtf)
    path(genome_fasta)

    output:
    path("isoquant_out")

    script:
    """
    source /scratch/nxu/SFARI/.virtualenvs/isoquant/bin/activate
    HOME="\$SLURM_TMPDIR" isoquant \\
        --reference ${genome_fasta} \\
        --genedb ${annotation_gtf} \\
        --complete_genedb \\
        --bam ${minimap_bam.join(' ')} \\
        --data_type pacbio_ccs \\
        --fl_data \\
        --threads ${task.cpus} \\
        --output isoquant_out
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
    path("supportedTranscriptModels.gtf"), emit: supported_tx_gtf
    path("novelTranscripts.gtf"), emit: novel_tx_gtf
    path("bambu_result.rds"), emit: bambu_rds
    script:
    """
    export R_LIBS="\${SCRATCH}/R/\${EBVERSIONR}"
    run_bambu.R ${annotation_gtf} "*.aligned.bam" ${genome_fasta} ${task.cpus} bambu_result.rds
    """
}

workflow COMPARE_OTHER_TRANSCRIPT_DISCOVERY {
    take:
    flnc_bam
    annotation_gtf
    genome_fasta

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
    isoquant(minimap2_genome.out.collect(), index_bam.out.collect(), annotation_gtf, genome_fasta)

    emit:
    supported_tx_gtf = bambu.out.supported_tx_gtf
    novel_tx_gtf     = bambu.out.novel_tx_gtf
    isoquant_out     = isoquant.out
}