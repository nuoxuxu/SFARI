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
        $predicted_cds_gtf \\
        "novel_exonic_regions.gtf"
    """
}

process novelSpliceSites {
    publishDir "${params.output_dir}/${params.orf_prediction}/splice_sites", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/env"

    input:
    path annotation_gtf
    path predicted_cds_gtf

    output:
    path "known_splice_sites_cds.csv", emit: known_splice_sites_cds_csv
    path "novel_splice_sites_cds.csv", emit: novel_splice_sites_cds_csv
    path "known_splice_sites_all.csv", emit: known_splice_sites_all_csv
    path "novel_splice_sites_all.csv", emit: novel_splice_sites_all_csv

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
        $predicted_cds_gtf \\
        "novel_CDS.gtf"
    """
}

process riboseqTrack {
    publishDir "${params.output_dir}/${params.orf_prediction}/UCSC_tracks", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/envs/r_env"

    input:
    path riboseq_file

    output:
    path "riboseq.gtf"

    script:
    """
    riboseq_track.R \\
        $riboseq_file
    """
}

process formatExonicRegions {
    publishDir "${params.output_dir}/${params.orf_prediction}/exonic_regions", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/envs/r_env"

    input:
    path annotation_gtf
    path predicted_cds_gtf
    path novel_CDS_gtf

    output:
    path "CDS_regions.csv", emit: CDS_regions_csv
    path "UTR_regions.csv", emit: UTR_regions_csv

    script:
    """
    format_exonic_region_files.R \\
        $annotation_gtf \\
        $predicted_cds_gtf \\
        $novel_CDS_gtf
    """
}


process addPhyloPToExonicRegions {
    publishDir "${params.output_dir}/${params.orf_prediction}/exonic_regions", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/env"

    input:
    path utr_regions_csv
    path cds_regions_csv
    path phyloP_bigwig

    output:
    path "UTR_regions_with_phyloP.csv"
    path "CDS_regions_with_phyloP.csv"

    script:
    """
    add_phyloP_to_exonic_regions.py \\
        --cds_regions_csv $cds_regions_csv \\
        --utr_regions_csv $utr_regions_csv \\
        --phyloP_bigwig $phyloP_bigwig
    """
}

process addPhyloPToSpliceSite {
    publishDir "${params.output_dir}/${params.orf_prediction}/splice_sites", mode: 'copy'
    
    conda "/scratch/s/shreejoy/nxu/SFARI/env"

    input:
    path known_splice_sites_cds_csv
    path novel_splice_sites_cds_csv
    path known_splice_sites_all_csv
    path novel_splice_sites_all_csv
    path phyloP_bigwig

    output:
    path "known_splice_sites_cds_phyloP.csv"
    path "novel_splice_sites_cds_phyloP.csv"
    path "known_splice_sites_all_phyloP.csv"
    path "novel_splice_sites_all_phyloP.csv"

    script:
    """
    add_phyloP_to_splice_sites.py \\
        --known_splice_sites_cds $known_splice_sites_cds_csv \\
        --novel_splice_sites_cds $novel_splice_sites_cds_csv \\
        --known_splice_sites_all $known_splice_sites_all_csv \\
        --novel_splice_sites_all $novel_splice_sites_all_csv \\
        --phyloP_bigwig $phyloP_bigwig
    """
}

workflow UCSCTracks {
    take:
    predicted_cds_gtf
    annotation_gtf

    main:
    novelExonicRegions(predicted_cds_gtf, annotation_gtf)
    novelSpliceSites(annotation_gtf, predicted_cds_gtf)
    novelCDS(predicted_cds_gtf, annotation_gtf)
    riboseqTrack(params.riboseq_file)
}

workflow PreprocessFigure6Files {
    take:
    predicted_cds_gtf
    annotation_gtf    

    main:
    novelExonicRegions(predicted_cds_gtf, annotation_gtf)
    novelSpliceSites(annotation_gtf, predicted_cds_gtf)
    novelCDS(predicted_cds_gtf, annotation_gtf)    
    formatExonicRegions(annotation_gtf, predicted_cds_gtf, novelCDS.out)
    addPhyloPToExonicRegions(formatExonicRegions.out.UTR_regions_csv, formatExonicRegions.out.CDS_regions_csv, params.phyloP_bigwig)
    addPhyloPToSpliceSite(
        novelSpliceSites.out.known_splice_sites_cds_csv, 
        novelSpliceSites.out.novel_splice_sites_cds_csv,
        novelSpliceSites.out.known_splice_sites_all_csv,
        novelSpliceSites.out.novel_splice_sites_all_csv,
        params.phyloP_bigwig)
}

workflow {
    predicted_cds_gtf = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/orfanage/orfanage.gtf"
    annotation_gtf = params.annotation_gtf

    novelExonicRegions(annotation_gtf, predicted_cds_gtf)
    novelSpliceSites(annotation_gtf, predicted_cds_gtf)
    novelCDS(annotation_gtf, predicted_cds_gtf)
    riboseqTrack(params.riboseq_file)
    formatExonicRegions(annotation_gtf, predicted_cds_gtf, novelCDS.out)
    addPhyloPToExonicRegions(formatExonicRegions.out.UTR_regions_csv, formatExonicRegions.out.CDS_regions_csv, params.phyloP_bigwig)
    addPhyloPToSpliceSite(
        novelSpliceSites.out.known_splice_sites_cds_csv, 
        novelSpliceSites.out.novel_splice_sites_cds_csv,
        novelSpliceSites.out.known_splice_sites_all_csv,
        novelSpliceSites.out.novel_splice_sites_all_csv,
        params.phyloP_bigwig)
}