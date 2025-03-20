include {merge_and_collapse} from 'subworkflows/local/merge_and_collapse'
include {classify_and_count} from 'subworkflows/local/classify_and_count'
include {ORFanage} from 'subworkflows/local/ORFanage'
include {proteoform_classification} from 'subworkflows/local/proteoform_classification'
include {proteomic} from 'subworkflows/local/proteomic'
include {peptide} from 'subworkflows/local/peptide'
include {SALMON_INDEX} from './modules/nf-core/salmon/index'
include {SALMON_QUANT} from './modules/nf-core/salmon/quant'
include {samplesheetToList} from 'plugin/nf-schema'

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

    merge_and_collapse(params.flnc_bam, params.mapped_bam)
    classify_and_count(merge_and_collapse.out.isoform_gff, merge_and_collapse.out.id_to_sample, merge_and_collapse.out.read_stat)
    ORFanage(classify_and_count.out.final_sample_gtf, classify_and_count.out.final_sample_classification, classify_and_count.out.final_sample_fasta)
    proteoform_classification(ORFanage.out.predicted_cds_gtf, ORFanage.out.peptide_fasta, ORFanage.out.best_orf)
    proteomic(proteoform_classification.out.protein_database, params.mzXMLfiles)
    peptide(classify_and_count.out.final_sample_classification, ORFanage.out.predicted_cds_gtf, proteoform_classification.out.protein_database, proteomic.out.peptides)

    extractFinalTranscriptsFasta(params.genome_fasta, classify_and_count.out.final_sample_gtf)
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
    SALMON_QUANT(input_reads, SALMON_INDEX.out.index, classify_and_count.out.final_sample_gtf, extractFinalTranscriptsFasta.out, false, "")    
}    