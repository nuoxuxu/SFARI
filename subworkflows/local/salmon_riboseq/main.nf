include {SALMON_INDEX} from '../../../modules/nf-core/salmon/index'
include {SALMON_QUANT} from '../../../modules/nf-core/salmon/quant'
include {samplesheetToList} from 'plugin/nf-schema' 

workflow {
    SALMON_INDEX(params.genome_fasta, params.transcript_fasta)
    input_reads = Channel
        .fromList(samplesheetToList("/scratch/s/shreejoy/nxu/SFARI/assets/riboseq.csv", "/scratch/s/shreejoy/nxu/SFARI/assets/schema_input.json"))
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
    SALMON_QUANT(input_reads, SALMON_INDEX.out.index, params.annotation_gtf, params.transcript_fasta, false, "")
}