#!/usr/bin/env nextflow

process renameCdsToExon {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    label "short_slurm_job"

    input:
    path predicted_cds_gtf
    path annotation_gtf

    output:
    path "SFARI.cds_renamed_exon.gtf", emit: sample_cds_renamed
    path "SFARI.transcript_exons_only.gtf", emit: sample_transcript_exon_only
    path "gencode.cds_renamed_exon.gtf", emit: ref_cds_renamed
    path "gencode.transcript_exons_only.gtf", emit: ref_transcript_exon_only

    script:
    """
    rename_cds_to_exon.py \\
        --sample_gtf $predicted_cds_gtf \\
        --sample_name "SFARI" \\
        --reference_gtf $annotation_gtf \\
        --reference_name "gencode" \\
        --num_cores $task.cpus
    """
}

process sqantiProtein {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    label "short_slurm_job"

    input:
    path sample_cds_renamed
    path sample_transcript_exon_only
    path ref_cds_renamed
    path ref_transcript_exon_only
    path best_orf
    path predicted_cds_gtf
    
    output:
    path "${predicted_cds_gtf}.sqanti_protein_classification.tsv"

    script:

    """
    sqanti3_protein.py \\
        "SFARI.transcript_exons_only.gtf" \\
        "SFARI.cds_renamed_exon.gtf" \\
        $best_orf \\
        "gencode.transcript_exons_only.gtf" \\
        "gencode.cds_renamed_exon.gtf" \\
        -d ./ \\
        -p $predicted_cds_gtf
    """

}

process fivePrimeUtr {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    label "short_slurm_job"

    input:
    path annotation_gtf
    path predicted_cds_gtf
    path protein_classification

    output:
    path "${predicted_cds_gtf}.sqanti_protein_classification._w_5utr_info.tsv"

    script:
    """
    1_get_gc_exon_and_5utr_info.py \\
        --gencode_gtf $annotation_gtf \\
        --odir ./

    2_classify_5utr_status.py \\
        --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \\
        --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \\
        --sample_cds_gtf $predicted_cds_gtf \\
        --odir ./
    
    3_merge_5utr_info_to_pclass_table.py \\
        --utr_info pb_5utr_categories.tsv \\
        --sqanti_protein_classification $protein_classification \\
        --odir ./  
    """
}

process proteinClassification {

    publishDir "${params.output_dir}/${params.orf_prediction}", mode: 'copy'
    label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    path protein_classification_w_5utr

    output:
    path "SFARI.protein_classification.tsv"

    script:
    """
    protein_classification.py \\
        --sqanti_protein $protein_classification_w_5utr \\
        --name "SFARI"
    """
}

process makeProteinSearchDatabase {
    publishDir "${params.output_dir}/${params.orf_prediction}", mode: 'copy'

    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    val mode
    path protein_classification
    path filtered_predicted_cds_gtf
    path peptide_fasta
    path translation_fasta

    output:
    path "*.fasta"

    script:
    """
    make_protein_search_database.py \\
        --mode $mode \\
        --protein_classification $protein_classification \\
        --filtered_predicted_cds_gtf $filtered_predicted_cds_gtf \\
        --peptide_fasta $peptide_fasta \\
        --translation_fasta $translation_fasta \\
        --output "${mode}.fasta"
    """
}

workflow proteoform_classification {
    take:
    predicted_cds_gtf
    peptide_fasta
    best_orf
    main:
    renameCdsToExon(predicted_cds_gtf, params.annotation_gtf)
    sqantiProtein(renameCdsToExon.out.sample_cds_renamed, renameCdsToExon.out.sample_transcript_exon_only, renameCdsToExon.out.ref_cds_renamed, renameCdsToExon.out.ref_transcript_exon_only, best_orf, predicted_cds_gtf)
    fivePrimeUtr(params.annotation_gtf, predicted_cds_gtf, sqantiProtein.out)
    proteinClassification(fivePrimeUtr.out)
    makeProteinSearchDatabase(params.searchDB, proteinClassification.out, predicted_cds_gtf, peptide_fasta, params.translation_fasta)
    emit:
    protein_database = makeProteinSearchDatabase.out
}