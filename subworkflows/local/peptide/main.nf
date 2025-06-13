process peptideTrackUCSC {
    
    conda "$moduleDir/python.yml"
    
    input:
    val mode
    path annotation_gtf
    path final_sample_classification
    path predicted_cds_gtf
    path protein_search_database
    path peptides

    output:
    path "peptides_${mode}.gtf"

    script:
    """
    make_peptide_gtf_file.py \\
        --annotation_gtf $annotation_gtf \\
        --final_sample_classification $final_sample_classification \\
        --predicted_cds_gtf $predicted_cds_gtf \\
        --protein_search_database $protein_search_database \\
        --peptides $peptides \\
        --output "peptides_${mode}.gtf"
    """
}

process addPeptideAnnotation {
    publishDir "${params.output_dir}/${params.orf_prediction}/UCSC_tracks", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/envs/r_env"

    input:
    val mode
    path annotation_gtf
    path predicted_cds_gtf
    path peptides_gtf
    
    output:
    path "annot_peptides_${mode}.gtf"
    
    script:
    """
    get_novel_peptides.R \\
        $annotation_gtf \\
        $predicted_cds_gtf \\
        $peptides_gtf \\
        "annot_peptides_${mode}.gtf"
    """
}

process filterPeptidesUCSC {
    publishDir "${params.output_dir}/${params.orf_prediction}/UCSC_tracks", mode: 'copy'
    
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"

    input:
    path annot_peptides_gtf
    
    output:
    path "detected_peptides.gtf"
    
    script:
    """
    filter_peptides_UCSC.py \\
        --annot_peptides_gtf $annot_peptides_gtf \\
        --output "detected_peptides.gtf" \\
    """
}
process peptideMapping {
    publishDir "${params.output_dir}/${params.orf_prediction}", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/envs/r_env"

    input:
    path annotation_gtf
    path predicted_cds_gtf
    path peptides_gtf
    
    output:
    path "peptide_mapping.parquet"
    
    script:
    """
    peptide_mapping.R \\
        $annotation_gtf \\
        $predicted_cds_gtf \\
        $peptides_gtf \\
        "peptide_mapping.parquet"
    """
}

workflow peptide {
    take:
    final_sample_classification
    predicted_cds_gtf
    protein_search_database
    peptides
    main:
    peptideTrackUCSC(params.searchDB, params.annotation_gtf, final_sample_classification, predicted_cds_gtf, protein_search_database, peptides)
    addPeptideAnnotation(params.searchDB, params.annotation_gtf, predicted_cds_gtf, peptideTrackUCSC.out)
    filterPeptidesUCSC(addPeptideAnnotation.out)
    peptideMapping(params.annotation_gtf, predicted_cds_gtf, peptideTrackUCSC.out)
}

workflow {
    final_sample_classification = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/final_classification.parquet"
    predicted_cds_gtf = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/orfanage/orfanage.gtf"
    protein_search_database = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/orfanage/hybrid.fasta"
    peptides = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/orfanage/hybrid_percolator.tsv"

    peptideTrackUCSC(params.searchDB, params.annotation_gtf, final_sample_classification, predicted_cds_gtf, protein_search_database, peptides)
    addPeptideAnnotation(params.searchDB, params.annotation_gtf, predicted_cds_gtf, peptideTrackUCSC.out)
    filterPeptidesUCSC(addPeptideAnnotation.out)
    peptideMapping(params.annotation_gtf, predicted_cds_gtf, peptideTrackUCSC.out)
}