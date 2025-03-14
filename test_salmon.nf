include {SALMON_INDEX} from './modules/nf-core/salmon/index'

params.genome_fasta = "/Users/xunuo/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa"
params.transcript_fasta = "/Users/xunuo/Genomic_references/GENCODE/gencode.v47.annotation.gtf"
workflow {
    SALMON_INDEX(params.genome_fasta, params.transcript_fasta)
}