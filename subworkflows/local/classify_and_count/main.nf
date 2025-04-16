#!/usr/bin/env nextflow

process pigeonPrepare {
    conda "${moduleDir}/environment.yml"

    input:
    path isoform_gff
    path annotation_gtf
    path genome_fasta

    output:
    path "merged_collapsed.sorted.gff", emit: sorted_isoform_gff
    path "merged_collapsed.sorted.gff.pgi", emit: sorted_isoform_gff_pgi
    path "*.annotation.sorted.gtf", emit: sorted_annotation
    path "*.annotation.sorted.gtf.pgi", emit: sorted_annotation_gtf_pgi
    path "GRCh38.primary_assembly.genome.fa.fai", emit: reference_fasta_pgi

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon prepare $isoform_gff $annotation_gtf $genome_fasta
    """
}

process pigeonClassify {
    label "short_slurm_job"

    conda "${moduleDir}/environment.yml"
    
    input:
    path sorted_isoform_gff
    path sorted_isoform_gff_pgi
    path sorted_annotation
    path sorted_annotation_gtf_pgi
    path genome_fasta
    path reference_fasta_pgi

    output:
    path "merged_collapsed_classification.txt", emit: classification
    path "merged_collapsed_junctions.txt", emit: junctions
    path "merged_collapsed.report.json"
    path "merged_collapsed.summary.txt"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon classify $sorted_isoform_gff $sorted_annotation $genome_fasta
    """
}

process pigeonFilter {
    publishDir "${params.output_dir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"
    
    input:
    path classification
    path junctions
    path isoform_gff

    output:
    
    path "merged_collapsed_classification.filtered_lite_classification.txt", emit: filtered_classification
    path "merged_collapsed_classification.filtered_lite_junctions.txt"
    path "merged_collapsed.sorted.filtered_lite.gff", emit: filtered_gff
    path "merged_collapsed_classification.filtered_lite_reasons.txt"
    path "merged_collapsed_classification.filtered.summary.txt"
    path "merged_collapsed_classification.filtered.report.json"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon filter $classification --isoforms $isoform_gff
    """
}

process getFullExpression {
    label "short_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"

    input:
    path read_stat
    path id_to_sample

    output:
    path "full_expression.parquet", emit: full_expression

    script:
    """
    get_full_expression.py \\
        --read_stat $read_stat
        --id_to_sample $id_to_sample
    """
}

process filterByExpressionExternalSupport {
    label "short_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'

    conda "${moduleDir}/environment.yml"

    input:
    path classification
    path filtered_gff
    path full_expression
    val min_reads
    val min_n_sample
    path polyA_site
    path refTSS
    
    output:
    path "final_classification.parquet", emit: final_sample_classification
    path "final_transcripts.gtf", emit: final_sample_gtf
    path "final_expression.parquet", emit: final_sample_expression
    
    script:
    """
    filter_by_exp_ext.py \\
        --classification $classification \\
        --filtered_gff $filtered_gff \\
        --full_expression $full_expression \\        
        --min_reads $min_reads \\
        --min_n_sample $min_n_sample \\
        --polyA_site $polyA_site \\
        --refTSS $refTSS \\
    """
}

process extractTranscriptFasta {
    conda "${moduleDir}/environment.yml"

    input:
    path final_sample_gtf
    path genome_fasta

    output:
    path "final_transcripts.fasta", emit: final_sample_fasta

    script:
    """
    gffread -w final_transcripts.fasta -g $genome_fasta $final_sample_gtf
    """
}

workflow classify_and_count {
    take:
    isoform_gff
    id_to_sample
    read_stat

    main:
    pigeonPrepare(isoform_gff, params.annotation_gtf, params.genome_fasta)
    pigeonClassify(pigeonPrepare.out.sorted_isoform_gff, pigeonPrepare.out.sorted_isoform_gff_pgi, pigeonPrepare.out.sorted_annotation, pigeonPrepare.out.sorted_annotation_gtf_pgi, params.genome_fasta, pigeonPrepare.out.reference_fasta_pgi)
    pigeonFilter(pigeonClassify.out.classification, pigeonClassify.out.junctions, pigeonPrepare.out.sorted_isoform_gff)
    getFullExpression(read_stat, id_to_sample)
    filterByExpressionExternalSupport(pigeonFilter.out.filtered_classification, pigeonFilter.out.filtered_gff, getFullExpression.out, params.min_reads, params.min_n_sample, params.polyA_site, params.refTSS)
    extractTranscriptFasta(filterByExpressionExternalSupport.out.final_sample_gtf, params.genome_fasta)

    emit:
    final_sample_gtf = filterByExpressionExternalSupport.out.final_sample_gtf
    final_sample_classification = filterByExpressionExternalSupport.out.final_sample_classification
    final_sample_fasta = extractTranscriptFasta.out.final_sample_fasta
}

workflow {
    pigeonPrepare(params.isoform_gff, params.annotation_gtf, params.genome_fasta)
    pigeonClassify(pigeonPrepare.out.sorted_isoform_gff, pigeonPrepare.out.sorted_isoform_gff_pgi, pigeonPrepare.out.sorted_annotation, pigeonPrepare.out.sorted_annotation_gtf_pgi, params.genome_fasta, pigeonPrepare.out.reference_fasta_pgi)
    pigeonFilter(pigeonClassify.out.classification, pigeonClassify.out.junctions, pigeonPrepare.out.sorted_isoform_gff)
    filterByExpressionExternalSupport(pigeonFilter.out.filtered_classification, pigeonFilter.out.filtered_gff, getFullExpression.out, params.min_reads, params.min_n_sample, params.polyA_site, params.refTSS)
    extractTranscriptFasta(filterByExpressionExternalSupport.out.final_sample_gtf, params.genome_fasta)
}