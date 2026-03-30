process filterBySampleCount {
    label "short_slurm_job"
    storeDir "${params.output_dir}/compare_other_LRS_datasets/filtered"
    module 'StdEnv/2023:python/3.11.5:gcc/12.3:arrow/19.0.1:rust/1.85.0'

    input:
    path(bed_files)
    val(min_samples)

    output:
    path("joglekar_filtered.bed")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    filter_joglekar_by_sample_count.py ${min_samples} joglekar_filtered.bed
    """
}

process bedToGtf {
    label "short_slurm_job"
    storeDir "${params.output_dir}/compare_other_LRS_datasets/gtf"
    module 'StdEnv/2023:kent_tools/486'

    input:
    path(bed_file)

    output:
    path("joglekar.gtf")

    script:
    """
    bedToGenePred ${bed_file} joglekar.genePred
    genePredToGtf file joglekar.genePred joglekar.gtf
    """
}

process liftPatowaryGtf {
    label "short_slurm_job"
    storeDir "${params.output_dir}/compare_other_LRS_datasets"
    module 'StdEnv/2023:python/3.11.5:gcc/12.3:arrow/19.0.1:rust/1.85.0'

    input:
    path(gtf)
    path(chain_file)

    output:
    path("patowary_hg38.gtf")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    CrossMap gff ${chain_file} ${gtf} patowary_hg38.gtf
    """
}

process sqanti3QC {
    label "mid_slurm_job"
    storeDir "${params.output_dir}/compare_other_LRS_datasets"
    container "sqanti3_latest.sif"

    input:
    path(gtf)
    path annotation_gtf
    path genome_fasta
    path refTSS
    path polyA_site

    output:
    path("sqanti3_results")

    script:
    """
    export PATH=/conda/miniconda3/envs/sqanti3/bin:\$PATH
    sqanti3_qc.py \\
        -t ${task.cpus} \\
        --skipORF \\
        --output joglekar \\
        --dir sqanti3_results \\
        --CAGE_peak ${refTSS} \\
        --polyA_peak ${polyA_site} \\
        --isoforms ${gtf} \\
        --refGTF ${annotation_gtf} \\
        --refFasta ${genome_fasta}
    """
}

process makeVennDiagrams {
    label "short_slurm_job"
    storeDir "${params.output_dir}/compare_other_LRS_datasets/venn"
    module 'StdEnv/2023:python/3.11.5:gcc/12.3:arrow/19.0.1:rust/1.85.0'

    input:
    path sfari_classification
    path sfari_gtf
    path patowary_classification
    path patowary_corrected_gtf
    path joglekar_sqanti_dir
    path encode4_gtf
    path encode4_abundance

    output:
    path("upset_*.pdf")
    path("venn_*.pdf")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    comapre_LRS_datasets_venn.py \\
        --sfari_classification ${sfari_classification} \\
        --sfari_gtf ${sfari_gtf} \\
        --patowary_classification ${patowary_classification} \\
        --patowary_gtf ${patowary_corrected_gtf} \\
        --joglekar_classification ${joglekar_sqanti_dir}/joglekar_classification.txt \\
        --joglekar_gtf ${joglekar_sqanti_dir}/joglekar_corrected.gtf \\
        --encode4_gtf ${encode4_gtf} \\
        --encode4_abundance ${encode4_abundance}
    """
}

workflow COMPARE_OTHER_LRS_DATASETS {
    take:
    bed_files
    sfari_classification
    sfari_gtf
    annotation_gtf
    genome_fasta
    refTSS
    polyA_site
    patowary_classification
    patowary_gtf
    chain_file
    encode4_gtf
    encode4_abundance

    main:
    filterBySampleCount(bed_files.collect(), params.joglekar_min_samples)
    bedToGtf(filterBySampleCount.out)
    sqanti3QC(
        bedToGtf.out,
        annotation_gtf,
        genome_fasta,
        refTSS,
        polyA_site
    )
    liftPatowaryGtf(patowary_gtf, chain_file)

    makeVennDiagrams(
        sfari_classification,
        sfari_gtf,
        patowary_classification,
        liftPatowaryGtf.out,
        sqanti3QC.out,
        encode4_gtf,
        encode4_abundance
    )
}
