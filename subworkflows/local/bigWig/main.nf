process samtools_merge {
    module "StdEnv/2023:samtools/1.22.1"
    label "short_slurm_job"
    storeDir "nextflow_results/bigWig"

    input:
    tuple val(time_point), path("mapped_bam")

    output:
    tuple val(time_point), path("${time_point}_merged.bam")
    
    script:
    """
    samtools merge mapped_bam* -o "${time_point}_merged.bam" -@ 40 -f
    """
}

process tobigWig {
    module "StdEnv/2023:bedtools/2.31.0:kent_tools/486"
    label "short_slurm_job"
    storeDir "nextflow_results/bigWig"

    input:
    tuple val(time_point), path(merged_bam), path(genome_sizes)

    output:
    path("${time_point}_merged.bw")
    
    script:
    """
    genomeCoverageBed -ibam $merged_bam -bg -split -g genome.sizes > output.bedgraph

    bedGraphToBigWig output.bedgraph $genome_sizes "${time_point}_merged.bw"
    """    
}

workflow bigWig {

    channel.value(file("/scratch/nxu/SFARI/data/hg38.p13.chrom.sizes")).set { genome_sizes }

    t00_bam_ch = channel.fromPath("data/long_read/7T3PQ0T/LUO26876.20240514/iPSC_*/outputs/mapped.bam").collect()
        .collect()
        .map { bams -> ["t00", bams] }
    t04_bam_ch = channel.fromPath("data/long_read/7T3PQ0T/LUO26876.20240514/NPC_*/outputs/mapped.bam").collect()
        .collect()
        .map { bams -> ["t04", bams] }
    t30_bam_ch = channel.fromPath("data/long_read/7T3PQ0T/LUO26876.20240514/CN_*/outputs/mapped.bam")
        .collect()
        .map { bams -> ["t30", bams] }

    t00_bam_ch.mix(t04_bam_ch).mix(t30_bam_ch)
        | samtools_merge
    samtools_merge.out.combine(genome_sizes)
        | tobigWig
}