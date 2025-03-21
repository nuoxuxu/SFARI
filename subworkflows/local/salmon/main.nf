include {SALMON_INDEX} from '../../../modules/nf-core/salmon/index'
include {SALMON_QUANT} from '../../../modules/nf-core/salmon/quant'
include {samplesheetToList} from 'plugin/nf-schema'

params.final_sample_gtf = "nextflow_results/V47/final_transcripts.gtf"

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
    extractFinalTranscriptsFasta(params.genome_fasta, params.final_sample_gtf)
    SALMON_INDEX(params.genome_fasta, extractFinalTranscriptsFasta.out)
    input_reads = Channel
        .fromList(samplesheetToList("assets/samplesheet.csv", "${projectDir}/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
    SALMON_QUANT(input_reads, SALMON_INDEX.out.index, params.final_sample_gtf, extractFinalTranscriptsFasta.out, false, "")
}