include {merge_and_collapse} from './subworkflows/local/merge_and_collapse'
include {classify_and_count} from './subworkflows/local/classify_and_count'
include {ORFanage} from './subworkflows/local/ORFanage'
include {proteoform_classification} from './subworkflows/local/proteoform_classification'
include {proteomic} from './subworkflows/local/proteomic'
include {peptide} from './subworkflows/local/peptide'
include {UCSCTracks} from './subworkflows/local/variant_analysis'
include {PreprocessFigure6Files} from './subworkflows/local/variant_analysis'
include { COMPARE_OTHER_ORF_CALLING } from './subworkflows/local/compare_other_ORF_calling'
include { COMPARE_OTHER_TRANSCRIPT_DISCOVERY } from './subworkflows/local/compare_other_transcript_discovery'
include { bigWig } from './subworkflows/local/bigWig'
include { COMPARE_OTHER_LRS_DATASETS } from './subworkflows/local/compare_other_LRS_datasets'
include { ALT_START_CODON } from './subworkflows/local/alternative_start_codon'
include { DUFFY } from './subworkflows/local/duffy_correlation'

workflow {
    channel.fromPath("data/long_read/7T3PQ0T/LUO26876.20240514/*/outputs/flnc.bam").set { flnc_bams }
    merge_and_collapse(params.flnc_bam, params.mapped_bam)
    classify_and_count(merge_and_collapse.out.isoform_gff, merge_and_collapse.out.id_to_sample, merge_and_collapse.out.read_stat)
    ORFanage(classify_and_count.out.final_sample_gtf, classify_and_count.out.final_sample_classification, classify_and_count.out.final_sample_fasta)
    // proteoform_classification(ORFanage.out.predicted_cds_gtf, ORFanage.out.peptide_fasta, ORFanage.out.best_orf)
    // proteomic(classify_and_count.out.final_sample_classification, ORFanage.out.predicted_cds_gtf, proteoform_classification.out.protein_database)
    // peptide(classify_and_count.out.final_sample_classification, ORFanage.out.predicted_cds_gtf, proteoform_classification.out.protein_database, proteomic.out.peptides)
    // UCSCTracks(ORFanage.out.predicted_cds_gtf, params.annotation_gtf)
    // PreprocessFigure6Files(channel.fromPath("nextflow_results/V47/orfanage/orfanage.gtf"), params.annotation_gtf)

    // channel.value(file("data/Human_Hexamer.tsv")).set{ Human_Hexamer }
    // channel.value(file("data/Human_logitModel.RData")).set{ Human_logitModel }
    COMPARE_OTHER_ORF_CALLING(classify_and_count.out.final_sample_fasta, ORFanage.out.orfanage_cds, Human_Hexamer, Human_logitModel)
    channel.fromPath("nextflow_results/STAR/*.SJ.out.tab").collect().set { star_sj_files }
    COMPARE_OTHER_TRANSCRIPT_DISCOVERY(flnc_bams, params.annotation_gtf, params.genome_fasta, classify_and_count.out.filtered_lite_gff, classify_and_count.out.final_sample_classification, classify_and_count.out.full_expression, classify_and_count.out.filtered_lite_classification, star_sj_files, params.refTSS, params.polyA_site)
    COMPARE_OTHER_LRS_DATASETS(
        channel.fromPath("data/joglekar_2024/other/*_corrected_reads.bed"),
        classify_and_count.out.final_sample_classification,
        classify_and_count.out.final_sample_gtf,
        params.annotation_gtf,
        params.genome_fasta,
        params.refTSS,
        params.polyA_site,
        params.patowary_sqanti_classification,
        params.patowary_sqanti_gtf,
        params.chain_file,
        params.encode4_gtf,
        params.encode4_abundance,
        channel.value(file("${params.output_dir}/classify_and_count/final_expression_snappy.parquet")),
        channel.value(file("${params.output_dir}/LR_patowary_snappy.parquet"))
    )
    DUFFY(
        channel.value(file("${params.output_dir}/salmon_riboseq")),
        channel.value(file("export/STAR_results"))
    )
    ALT_START_CODON(
        COMPARE_OTHER_TRANSCRIPT_DISCOVERY.out.genome_bam,
        classify_and_count.out.final_sample_classification,
        ORFanage.out.predicted_cds_gtf,
        merge_and_collapse.out.read_stat,
        channel.value(file(params.refTSS)),
        channel.value(file(params.polyA_site)),
        channel.value(file(params.annotation_gtf)),
        channel.value(params.gene_of_interest)
    )
    
}
