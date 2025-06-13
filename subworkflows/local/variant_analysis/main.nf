process novelExonicRegions {
    publishDir "${params.output_dir}/${params.orf_prediction}/UCSC_tracks", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/envs/r_env"

    input:
    path annotation_gtf
    path predicted_cds_gtf

    output:
    path "novel_exonic_regions.gtf"

    script:
    """
    novel_exonic_regions.R \\
        $annotation_gtf \\
        $predicted_cds_gtf
    """
}

process novelSpliceSites {
    publishDir "${params.output_dir}/${params.orf_prediction}/splice_sites", mode: 'copy'
    
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"

    input:
    path annotation_gtf
    path predicted_cds_gtf

    output:
    path "novel_splice_sites_all.csv"
    path "novel_splice_sites_cds.csv"
    path "known_splice_sites_all.csv"
    path "known_splice_sites_cds.csv"

    script:
    """
    novel_splice_sites.py \\
        --annotation_gtf $annotation_gtf \\
        --predicted_cds_gtf $predicted_cds_gtf
    """
}

process novelCDS {
    publishDir "${params.output_dir}/${params.orf_prediction}/UCSC_tracks", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/envs/r_env"

    input:
    path annotation_gtf
    path predicted_cds_gtf

    output:
    path "novel_CDS.gtf"

    script:
    """
    novel_CDS.R \\
        $annotation_gtf \\
        $predicted_cds_gtf
    """
}

workflow UCSCTracks {
    take:
    predicted_cds_gtf
    annotation_gtf

    main:
    novelExonicRegions(predicted_cds_gtf, annotation_gtf)
    novelSpliceSites(predicted_cds_gtf, annotation_gtf)
    novelCDS(predicted_cds_gtf, annotation_gtf)
}

workflow {
    predicted_cds_gtf = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/orfanage/orfanage.gtf"
    annotation_gtf = params.annotation_gtf

    novelExonicRegions(annotation_gtf, predicted_cds_gtf)
    novelSpliceSites(annotation_gtf, predicted_cds_gtf)
    novelCDS(annotation_gtf, predicted_cds_gtf)
}