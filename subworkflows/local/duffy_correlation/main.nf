process duffyT30Correlation {
    label "short_slurm_job"
    storeDir "${params.output_dir}/duffy_correlation"
    module 'StdEnv/2023:python/3.11.5:gcc/12.3:arrow/19.0.1:rust/1.85.0'
    beforeScript "source /scratch/nxu/SFARI/.venv/bin/activate"

    input:
    path(salmon_dir)
    path(star_dir)

    output:
    path("duffy_t30_correlation.pdf")

    script:
    """
    duffy_t30_correlation.py \\
        --salmon_dir ${salmon_dir} \\
        --star_dir ${star_dir} \\
        --output_dir .
    """
}

workflow DUFFY {
    take:
    salmon_dir  // path to salmon riboseq results directory (contains SRR15175557–59 subdirs)
    star_dir    // path to STAR results directory (contains CN_*ReadsPerGene.out.tab files)

    main:
    duffyT30Correlation(salmon_dir, star_dir)

    emit:
    plot = duffyT30Correlation.out
}
