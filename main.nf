include {merge_and_collapse} from './subworkflows/local/merge_and_collapse'
include {classify_and_count} from './subworkflows/local/classify_and_count'
include {ORFanage} from './subworkflows/local/ORFanage'
include {proteoform_classification} from './subworkflows/local/proteoform_classification'
include {proteomic} from './subworkflows/local/proteomic'
include {peptide} from './subworkflows/local/peptide'
include {UCSCTracks} from './subworkflows/local/variant_analysis'
include {PreprocessFigure6Files} from './subworkflows/local/variant_analysis'

workflow {
    isoform_gff = channel.fromPath("nextflow_results/merged_collapsed.gff")
    id_to_sample = channel.fromPath("nextflow_results/id_to_sample.txt")
    read_stat = channel.fromPath("nextflow_results/merged_collapsed.read_stat.txt")
    merge_and_collapse(params.flnc_bam, params.mapped_bam)
    classify_and_count(merge_and_collapse.out.isoform_gff, merge_and_collapse.out.id_to_sample, merge_and_collapse.out.read_stat)
    ORFanage(classify_and_count.out.final_sample_gtf, classify_and_count.out.final_sample_classification, classify_and_count.out.final_sample_fasta)
    proteoform_classification(ORFanage.out.predicted_cds_gtf, ORFanage.out.peptide_fasta, ORFanage.out.best_orf)
    final_sample_classification = channel.fromPath("nextflow_results/V47/final_classification.parquet")
    predicted_cds_gtf = channel.fromPath("nextflow_results/V47/orfanage/orfanage.gtf")
    protein_database = channel.fromPath("nextflow_results/V47/orfanage/hybrid.fasta")    
    proteomic(final_sample_classification, predicted_cds_gtf, protein_database)
    peptide(classify_and_count.out.final_sample_classification, ORFanage.out.predicted_cds_gtf, proteoform_classification.out.protein_database, proteomic.out.peptides) 
    UCSCTracks(ORFanage.out.predicted_cds_gtf, params.annotation_gtf)
    PreprocessFigure6Files(channel.fromPath("nextflow_results/V47/orfanage/orfanage.gtf"), params.annotation_gtf)
}