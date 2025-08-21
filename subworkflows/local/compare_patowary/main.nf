process liftOverTalonGtf {
    conda '/scratch/s/shreejoy/nxu/SFARI/env'
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
    conda '/scratch/s/shreejoy/nxu/SFARI/env'
    publishDir "${params.output_dir}/compare", mode: 'copy'

    input:
    path annotation_gtf
    path talon_gtf_lifted

    output:
    path 'TALON_GENCODE_V47.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap'

    script:
    """
    gffcompare -r $annotation_gtf $talon_gtf_lifted -o TALON_GENCODE_V47
    """
}

process getTalonSfariRefmap {
    conda '/scratch/s/shreejoy/nxu/SFARI/env'
    publishDir "${params.output_dir}/compare", mode: 'copy'

    input:
    path final_transcripts_gtf
    path talon_gtf_lifted

    output:
    path 'TALON_SFARI.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap'

    script:
    """
    gffcompare -r $final_transcripts_gtf $talon_gtf_lifted -o TALON_SFARI
    """
}
// process readRefmap {
//     conda '/scratch/s/shreejoy/nxu/SFARI/env'

//     input:
//     path annotation_gtf
//     path refmap
//     path talon_gtf_lifted
//     path final_transcripts_gtf

//     output:
//     path "TALON_GENCODE_V39.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap.reads.csv", emit: reads_csv

// }

workflow compare_patowary {
    take:
    talon_gtf
    chain_file
    annotation_gtf

    main:
    liftOverTalonGtf(talon_gtf, chain_file)
    getTalonGencodeV47Refmap(annotation_gtf, liftOverTalonGtf.out)
}

workflow {
    talon_gtf = channel.fromPath("/scratch/s/shreejoy/nxu/SFARI/data/cp_vz_0.75_min_7_recovery_talon.gtf")
    chain_file = channel.fromPath("/scratch/s/shreejoy/nxu/SFARI/data/hg19ToHg38.over.chain.gz")
    final_transcripts_gtf = channel.fromPath("/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/final_transcripts.gtf")
    liftOverTalonGtf(talon_gtf, chain_file)
    getTalonGencodeV47Refmap(params.annotation_gtf, liftOverTalonGtf.out)
    getTalonSfariRefmap(final_transcripts_gtf, liftOverTalonGtf.out)
}