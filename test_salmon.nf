include {SALMON_INDEX} from './modules/nf-core/salmon/index'

params.genome_fasta = "/project/s/shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa"
params.final_transcript_gtf = "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/V47/final_transcripts.gtf"

process extractFinalTranscriptsFasta {
    label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path genome_fasta
    path final_transcript_gtf
    output:
    path "final_transcripts.fasta"
    script:
    """
    gffread $final_transcript_gtf -g $genome_fasta -w final_transcripts.fasta
    """

}
workflow {
    extractFinalTranscriptsFasta(params.genome_fasta, params.final_transcript_gtf)
    SALMON_INDEX(params.genome_fasta, extractFinalTranscriptsFasta.out)
}