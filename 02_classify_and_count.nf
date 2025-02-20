#!/usr/bin/env nextflow

process pigeonPrepare {
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

process filterByExpression {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"
    label "short_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path annotation_gtf
    path genome_fasta
    path id_to_sample
    path classification
    path read_stat
    path filtered_gff
    val min_reads
    val min_n_sample
    
    output:
    path "final_classification.parquet", emit: final_sample_classification
    path "final_transcripts.gtf", emit: final_sample_gtf
    path "final_expression.parquet", emit: final_sample_expression
    path "fulll_expression.parquet"
    path "final_transcripts.fasta", emit: final_sample_fasta
    
    script:
    """
    filter_by_expression.py \\
        --annotation_gtf $annotation_gtf \\
        --id_to_sample $id_to_sample \\
        --classification $classification \\
        --read_stat $read_stat \\
        --filtered_gff $filtered_gff \\
        --min_reads $min_reads \\
        --min_n_sample $min_n_sample \\
        --classification_output "final_classification.parquet" \\
        --gtf_output "final_transcripts.gtf" \\
        --expression_output "final_expression.parquet"

    ~/miniforge3/envs/SQANTI3.env/bin/gffread -w final_transcripts.fasta -g $genome_fasta "final_transcripts.gtf"
    """
}

workflow {
    isoform_gff = Channel.fromPath("proc/merged_collapsed.gff")
    getIDToSample = channel.fromPath("proc/id_to_sample.txt")
    read_stat = channel.fromPath("proc/merged_collapsed.read_stat.txt")
    pigeonPrepare(isoform_gff, params.annotation_gtf, params.genome_fasta)
    pigeonClassify(pigeonPrepare.out.sorted_isoform_gff, pigeonPrepare.out.sorted_isoform_gff_pgi, pigeonPrepare.out.sorted_annotation, pigeonPrepare.out.sorted_annotation_gtf_pgi, params.genome_fasta, pigeonPrepare.out.reference_fasta_pgi)
    pigeonFilter(pigeonClassify.out.classification, pigeonClassify.out.junctions, pigeonPrepare.out.sorted_isoform_gff)
    filterByExpression(params.annotation_gtf, params.genome_fasta, getIDToSample, pigeonFilter.out.filtered_classification, read_stat, pigeonFilter.out.filtered_gff, params.min_reads, params.min_n_sample)    
}