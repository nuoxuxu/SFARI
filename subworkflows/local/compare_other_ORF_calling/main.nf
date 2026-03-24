process GeneMarkST {
    storeDir "nextflow_results/compare_other_ORF_calling"
    label "short_slurm_job"

    input:
    path(final_transcripts_fasta)

    output:
    path("${final_transcripts_fasta}.lst")

    script:
    """
    /home/nxu/tools/gmst/gmst.pl $final_transcripts_fasta
    """
}

process CPAT {
    module "StdEnv/2023:gcc/12.3:python/3.12.4:scipy-stack/2026a:r/4.5.0"
    label "short_slurm_job"
    storeDir "subworkflows/local/compare_other_ORF_calling"

    input:
    
    path(Human_Hexamer)
    path(Human_logitModel)
    path(final_transcripts_fasta)

    output:
    path("SFARI.ORF_seqs.fa"), emit: ORF_seqs
    path("SFARI.ORF_prob.tsv"), emit: ORF_prob
    path("SFARI.ORF_prob.best.tsv"), emit: ORF_prob_best
    path("SFARI.no_ORF.txt"), emit: no_ORF
    path("CPAT_run_info.log"), emit: log_file

    script:
    """
    source /scratch/nxu/SFARI/.virtualenvs/cpat/bin/activate
    cpat \
        -x $Human_Hexamer \
        -d $Human_logitModel \
        -g $final_transcripts_fasta \
        --min-orf=50 \
        --top-orf=50 \
        -o SFARI \
        1> SFARI_cpat.output \
        2> SFARI_cpat.error
    """
}

workflow compare_other_ORF_calling {
    channel.value(file("nextflow_results/classify_and_count/final_transcripts.fasta")).set{ final_transcripts_fasta }
    channel.value(file("data/Human_Hexamer.tsv")).set{ Human_Hexamer }
    channel.value(file("data/Human_logitModel.RData")).set{ Human_logitModel }
    GeneMarkST(final_transcripts_fasta)
    CPAT(Human_Hexamer, Human_logitModel, final_transcripts_fasta)
}
