#!/usr/bin/env nextflow

process runORFanage {
    label "short_slurm_job"

    conda "$moduleDir/environment.yml"

    input:
    path genome_fasta
    path annotation_gtf
    path final_sample_gtf

    output:
    path "orfanage_with_gene_id.gtf", emit: orfanage_gtf
    path "orfanage.stats"

    script:
    """
    gffread \\
        -g $genome_fasta \\
        --adj-stop \\
        -T -F -J \\
        -o corrected.gtf \\
        $final_sample_gtf

    orfanage \\
        --reference $genome_fasta \\
        --query corrected.gtf \\
        --output orfanage_without_gene_id.gtf \\
        --threads $task.cpus \\
        --minlen 50 \\
        --stats orfanage.stats \\
        $annotation_gtf

    gffread \\
        -g $genome_fasta \\
        --adj-stop \\
        -T -F -J -C \\
        -o orfanage_with_gene_id.gtf \\
        orfanage_without_gene_id.gtf
    """
}

process fixORFanageFormat {
    publishDir "${params.output_dir}/orfanage", mode: 'copy'
    label "short_slurm_job"

    container "quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0"

    input:
    path orfanage_gtf
    path genome_fasta

    output:
    path "orfanage.gtf"
    
    script:
    """
    agat_sp_add_start_and_stop.pl --gff $orfanage_gtf --fasta $genome_fasta --out "added_codons_orfanage_with_gene_id.gff3"

    agat_convert_sp_gff2gtf.pl --gff "added_codons_orfanage_with_gene_id.gff3" -o "orfanage.gtf" --gtf_version 3
    """
}

process extractORFanageCdsFasta {
    label "short_slurm_job"

    container "quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0"

    input:
    path genome_fasta
    path orfanage_gtf

    output:
    path "orfanage_cds.fasta", emit: orfanage_cds

    script:
    """
    agat_sp_extract_sequences.pl -g $orfanage_gtf -f $genome_fasta -t cds -o orfanage_cds.fasta
    """
}

process extractORFanageTranslationFasta {
    publishDir "${params.output_dir}/orfanage", mode: 'copy'
    label "short_slurm_job"

    container "quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0"

    input:
    path genome_fasta
    path orfanage_gtf

    output:
    path "orfanage_peptide.fasta", emit: orfanage_peptide

    script:
    """
    agat_sp_extract_sequences.pl -g $orfanage_gtf -f $genome_fasta -t cds -p -o orfanage_peptide.fasta
    """
}

process getBestOrfCsv {
    publishDir "${params.output_dir}/orfanage", mode: 'copy'

    conda "$moduleDir/environment.yml"

    input:
    path final_sample_classification
    path final_sample_fasta
    path orfanage_cds

    output:
    path "best_orf.tsv"

    script:
    """
    get_best_orf.py \\
        --final_sample_classification $final_sample_classification \\
        --final_sample_fasta $final_sample_fasta \\
        --cds_fasta $orfanage_cds \\
        --output "best_orf.tsv"
    """
}

workflow ORFanage {
    take:
    final_sample_gtf
    final_sample_classification
    final_sample_fasta
    main:
    runORFanage(params.genome_fasta, params.annotation_gtf, final_sample_gtf)
    fixORFanageFormat(runORFanage.out.orfanage_gtf, params.genome_fasta)
    extractORFanageCdsFasta(params.genome_fasta, fixORFanageFormat.out)
    extractORFanageTranslationFasta(params.genome_fasta, fixORFanageFormat.out)
    getBestOrfCsv(final_sample_classification, final_sample_fasta, extractORFanageCdsFasta.out)
    emit:
    predicted_cds_gtf = fixORFanageFormat.out
    peptide_fasta = extractORFanageTranslationFasta.out
    best_orf = getBestOrfCsv.out
}

workflow {

    final_sample_gtf = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/final_transcripts.gtf" 
    final_sample_classification = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/final_classification.parquet"
    final_sample_fasta = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/final_transcripts.fasta"

    runORFanage(params.genome_fasta, params.annotation_gtf, final_sample_gtf)
    fixORFanageFormat(runORFanage.out.orfanage_gtf, params.genome_fasta)
    extractORFanageCdsFasta(params.genome_fasta, fixORFanageFormat.out)
    extractORFanageTranslationFasta(params.genome_fasta, fixORFanageFormat.out)
    getBestOrfCsv(final_sample_classification, final_sample_fasta, extractORFanageCdsFasta.out)
}