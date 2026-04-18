#!/usr/bin/env nextflow
process BAMTOBED {
    label "short_slurm_job"
    module "StdEnv/2023:bedtools/2.31.0:samtools/1.22.1:python:gcc:arrow/19.0.1:rust"
    storeDir "nextflow_results/alternative_start_codon"

    input:
    path(genome_bam)
    path(final_classification)
    path(orfanage_gtf)
    path(read_stat)

    output:
    path("${genome_bam.baseName}.bed")

    script:
    """
    source /scratch/nxu/astrocytes/pytorch/bin/activate
    samtools view -F 2304 -b $genome_bam \\
    | bedtools bamtobed -i stdin -split \\
    | filter_reads.py --final_classification $final_classification --orfanage_gtf $orfanage_gtf --read_stat $read_stat --out "${genome_bam.baseName}.bed"
    """
}

process SPEARMAN_PER_SAMPLE {
    label "short_slurm_job"
    module "StdEnv/2023:bedtools/2.31.0:python:gcc:arrow/19.0.1:rust"
    storeDir "nextflow_results/alternative_start_codon"

    input:
    path(bed_file)
    path(cage_peaks)
    path(polya_sites)

    output:
    path("${bed_file.baseName}_spearman.tsv")

    script:
    """
    source /scratch/nxu/SFARI/.venv/bin/activate
    spearman_per_sample.py \\
        --bed $bed_file \\
        --cage_peaks $cage_peaks \\
        --polya_sites $polya_sites \\
        --out "${bed_file.baseName}_spearman.tsv"
    """
}

process PLOT_PITA_COUPLING {
    label "short_slurm_job"
    module "StdEnv/2023:python:gcc:arrow/19.0.1:rust"
    storeDir "figures/supplementary"

    input:
    path(spearman_tsvs)

    output:
    path("pita_coupling.pdf")

    script:
    """
    source /scratch/nxu/SFARI/.venv/bin/activate
    plot_pita_coupling.py \\
        --input $spearman_tsvs \\
        --outdir .
    """
}

process SUBSET_PITA_BAM {
    label "short_slurm_job"
    module "StdEnv/2023:samtools/1.22.1:python:gcc:arrow/19.0.1:rust"
    storeDir "nextflow_results/alternative_start_codon"

    input:
    tuple path(bam), path(bed_file), path(spearman_tsv)

    output:
    tuple path("${bam.baseName}.pita.bam"), path("${bam.baseName}.pita.bam.bai")

    script:
    """
    source /scratch/nxu/SFARI/.venv/bin/activate
    get_pita_read_names.py \\
        --spearman $spearman_tsv \\
        --bed $bed_file \\
        --out pita_read_names.txt
    samtools view -N pita_read_names.txt -b $bam > "${bam.baseName}.pita.bam"
    samtools index "${bam.baseName}.pita.bam"
    """
}

process GVIZ_PLOT {
    label "short_slurm_job"
    module "StdEnv/2023:gcc/12.3:r/4.5.0:r-bundle-bioconductor/3.21"
    storeDir "figures/PITA"

    input:
    tuple path(pita_bam), path(pita_bai)
    path(orfanage_gtf)
    val(gene_name)
    path(gencode_gtf)

    output:
    path("${gene_name}.pdf")

    script:
    """
    export R_LIBS=\$SCRATCH/R/\$EBVERSIONR/
    PITA_Gviz.R \\
        $pita_bam \\
        $orfanage_gtf \\
        $gene_name \\
        $gencode_gtf \\
        "${gene_name}.pdf"
    """
}

workflow ALT_START_CODON {
    take:
    genome_bam
    final_classification
    orfanage_gtf
    read_stat
    cage_peaks
    polya_sites
    gencode_gtf
    gene_of_interest

    main:
    BAMTOBED(genome_bam, final_classification, orfanage_gtf, read_stat)
    SPEARMAN_PER_SAMPLE(BAMTOBED.out, cage_peaks, polya_sites)
    PLOT_PITA_COUPLING(SPEARMAN_PER_SAMPLE.out.collect())

    bam_bed_spearman = genome_bam
        .map { bam -> [bam.baseName, bam] }
        .join(BAMTOBED.out.map { bed -> [bed.baseName, bed] })
        .join(SPEARMAN_PER_SAMPLE.out.map { tsv -> [tsv.baseName.replace("_spearman", ""), tsv] })
        .map { key, bam, bed, tsv -> tuple(bam, bed, tsv) }
    bam_bed_spearman
        .filter { bam, bed, tsv -> bam.simpleName == "CN_1_3" }
        | SUBSET_PITA_BAM

    GVIZ_PLOT(SUBSET_PITA_BAM.out, orfanage_gtf, gene_of_interest, gencode_gtf)
}
