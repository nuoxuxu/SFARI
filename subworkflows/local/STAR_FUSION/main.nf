process run_star {
    module "StdEnv/2023:star/2.7.11b:samtools/1.22.1"
    label "short_slurm_job"
    storeDir "nextflow_results/STAR"
    tag "${prefix}"
    
    input:
    tuple val(prefix), path(fastq_gz), path(star_genomeDir)
    
    output:
    tuple path("${prefix}.Aligned.sortedByCoord.out.bam"), path("${prefix}.Aligned.sortedByCoord.out.bam.bai"), emit: genome_bam
    path("${prefix}.SJ.out.tab"), emit: sj_tab
    path("${prefix}.Log.final.out"), emit: log_final_out
    path("${prefix}.Log.out"), emit: log_out
    path("${prefix}.Chimeric.out.junction"), emit: chimeric_junction
    
    script:
    """
    STAR --genomeDir $star_genomeDir \\
        --readFilesIn ${fastq_gz} \\
        --outReadsUnmapped None \\
        --twopassMode Basic \\
        --readFilesCommand "gunzip -c" \\
        --outSAMstrandField intronMotif \\
        --outSAMunmapped Within \\
        --chimSegmentMin 12 \\
        --chimJunctionOverhangMin 8 \\
        --chimOutJunctionFormat 1 \\
        --alignSJDBoverhangMin 10 \\
        --alignMatesGapMax 100000 \\
        --alignIntronMax 100000 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --outSAMattrRGline ID:GRPundef \\
        --chimMultimapScoreRange 3 \\
        --chimScoreJunctionNonGTAG -4 \\
        --chimMultimapNmax 20 \\
        --chimNonchimScoreDropMin 10 \\
        --peOverlapNbasesMin 12 \\
        --peOverlapMMp 0.1 \\
        --alignInsertionFlush Right \\
        --alignSplicedMateMapLminOverLmate 0 \\
        --alignSplicedMateMapLmin 30 \\
        --outSAMtype BAM SortedByCoordinate \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix "${prefix}."

    samtools index -@ ${task.cpus} "${prefix}.Aligned.sortedByCoord.out.bam"
    """
}

// process star_fusion {
//     container "star-fusion.v1.15.1.simg"
//     label "short_slurm_job"
//     script:
//     """
//     STAR-Fusion --genome_lib_dir /path/to/your/CTAT_resource_lib \\
//         -J Chimeric.out.junction \\
//         --output_dir star_fusion_outdir

//     """
// }

workflow STAR_FUSION {
    channel.value(file("/project/rrg-shreejoy/Genomic_references/GENCODE/STAR_index_v47")).set { star_genomeDir }

    channel
        .fromFilePairs('data/short_read/*_R{1,2}_001.fastq.gz')
        .combine(star_genomeDir)
        | run_star
}