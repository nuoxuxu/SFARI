process peptideTrackUCSC {
    storeDir "${params.output_dir}/peptide"
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

process peptideMapping {
    storeDir "${params.output_dir}/peptide"
    module 'StdEnv/2023:gcc/12.3:r/4.5.0:r-bundle-bioconductor/3.21'

    input:
    path annotation_gtf
    path predicted_cds_gtf
    path peptides_gtf

    output:
    path "peptide_mapping.parquet"

    script:
    """
    export R_LIBS=\$SCRATCH/R/\$EBVERSIONR
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
    peptideMapping(params.annotation_gtf, predicted_cds_gtf, peptideTrackUCSC.out)
}