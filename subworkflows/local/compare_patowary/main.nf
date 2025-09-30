process liftOverTalonGtf {
    container "crossmap.sif"
    publishDir "${params.output_dir}/compare", mode: 'copy'
    
    input:
    path talon_gtf
    path chain_file

    output:
    path "cp_vz_0.75_min_7_recovery_talon_hg38.gtf", emit: talon_gtf_lifted

    script:
    """
    CrossMap gff $chain_file $talon_gtf cp_vz_0.75_min_7_recovery_talon_hg38.gtf
    """
}

process getTalonGencodeV47Refmap {
    publishDir "${params.output_dir}/compare", mode: 'copy'

    input:
    path annotation_gtf
    path talon_gtf_lifted

    output:
    path 'TALON_GENCODE_V47.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap'

    script:
    """
    module load StdEnv/2023
    module load gffcompare/0.12.6
    gffcompare -r $annotation_gtf $talon_gtf_lifted -o TALON_GENCODE_V47
    """
}

process getTalonSfariRefmap {
    publishDir "${params.output_dir}/compare", mode: 'copy'

    input:
    path final_transcripts_gtf
    path talon_gtf_lifted

    output:
    path 'TALON_SFARI.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap'

    script:
    """
    module load StdEnv/2023
    module load gffcompare/0.12.6    
    gffcompare -r $final_transcripts_gtf $talon_gtf_lifted -o TALON_SFARI
    """
}

process addPAtowaryColumnToClassification {
    publishDir "${params.output_dir}/compare", mode: 'copy'

    input:
    path classification
    path TALON_SFARI_refmap

    output:
    path 'classification_with_patowary.txt'

    script:
    """
    compare_to_patowary.py --classification $classification --TALON_SFARI_refmap $TALON_SFARI_refmap --output classification_with_patowary.txt
    """
}

workflow compare_patowary {
    take:
    talon_gtf
    chain_file
    annotation_gtf
    final_transcripts_gtf

    main:
    liftOverTalonGtf(talon_gtf, chain_file)
    getTalonGencodeV47Refmap(annotation_gtf, liftOverTalonGtf.out)
    getTalonSfariRefmap(final_transcripts_gtf, liftOverTalonGtf.out)
    addPAtowaryColumnToClassification(params.classification, getTalonSfariRefmap.out)
}

workflow {
    liftOverTalonGtf(params.patowary_talon_gtf, params.chain_file)
    getTalonGencodeV47Refmap(params.annotation_gtf, liftOverTalonGtf.out)
    getTalonSfariRefmap(params.final_transcripts_gtf, liftOverTalonGtf.out)
    addPAtowaryColumnToClassification(params.classification, getTalonSfariRefmap.out)
}