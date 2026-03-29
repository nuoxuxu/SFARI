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
    storeDir "nextflow_results/compare_other_ORF_calling"

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

process extractGeneMarkCDS {
    label "short_slurm_job"
    storeDir "${params.output_dir}/compare_other_ORF_calling"

    input:
    path genemark_lst
    path transcripts_fasta

    output:
    path "genemark_coding_sequences.fasta", emit: genemark_cds

    script:
    """
    extract_cds.py \\
        --lst $genemark_lst \\
        --fasta $transcripts_fasta \\
        --output genemark_coding_sequences.fasta
    """
}

process subsetCpatOrfs {
    label "short_slurm_job"
    storeDir "${params.output_dir}/compare_other_ORF_calling"

    input:
    path orfs_fasta
    path best_tsv

    output:
    path "SFARI.ORF_seqs.best.fa", emit: cpat_best_orfs

    script:
    """
    subset_cpat_orfs.py \\
        --orfs_fasta $orfs_fasta \\
        --best_tsv $best_tsv \\
        --output SFARI.ORF_seqs.best.fa
    """
}

process vennOrfOverlap {
    label "short_slurm_job"
    storeDir "${params.output_dir}/compare_other_ORF_calling"

    input:
    path orfanage_cds
    path genemark_cds
    path cpat_best_orfs

    output:
    path "venn_orf_overlap.pdf"

    script:
    """
    venn_orf_overlap.py \\
        --orfanage $orfanage_cds \\
        --genemark $genemark_cds \\
        --cpat $cpat_best_orfs \\
        --output venn_orf_overlap.pdf
    """
}

workflow COMPARE_OTHER_ORF_CALLING {
    take:
    final_transcripts_fasta
    orfanage_cds
    Human_Hexamer
    Human_logitModel

    main:
    GeneMarkST(final_transcripts_fasta)
    CPAT(Human_Hexamer, Human_logitModel, final_transcripts_fasta)    
    extractGeneMarkCDS(GeneMarkST.out, final_transcripts_fasta)
    subsetCpatOrfs(CPAT.out.ORF_seqs, CPAT.out.ORF_prob_best)
    vennOrfOverlap(orfanage_cds, extractGeneMarkCDS.out.genemark_cds, subsetCpatOrfs.out.cpat_best_orfs)
}
