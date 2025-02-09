#!/usr/bin/env nextflow

// listOfFiles = file('data/long_read/LUO26876.20240514/*/outputs/flnc.bam')

// println listOfFiles

params.datadir = "/scratch/s/shreejoy/nxu/SFARI/data/"
params.comet_params = "/scratch/s/shreejoy/nxu/SFARI/data/comet.params"
params.hmmfile = "$projectDir/data/Pfam-A.hmm"
params.name = "SFARI"
params.min_junctions_after_stop_codon = 2

process pigeonPrepare {
    input:
    path isoform_gff
    path annotation_gtf
    path genome_fasta

    output:
    path "merged_collapsed.sorted.gff", emit: sorted_isoform_gff
    path "merged_collapsed.sorted.gff.pgi", emit: sorted_isoform_gff_pgi
    path "*.annotation.sorted.gtf", emit: sorted_annotation
    path "*.annotation.sorted.gtf.pgi", emit: sorted_annotation_gtf_pgi
    path "GRCh38.primary_assembly.genome.fa.fai", emit: reference_fasta_pgi

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon prepare $isoform_gff $annotation_gtf $genome_fasta
    """
}

process pigeonClassify {
    label "short_slurm_job"
    
    input:
    path sorted_isoform_gff
    path sorted_isoform_gff_pgi
    path sorted_annotation
    path sorted_annotation_gtf_pgi
    path genome_fasta
    path reference_fasta_pgi

    output:
    path "merged_collapsed_classification.txt", emit: classification
    path "merged_collapsed_junctions.txt", emit: junctions
    path "merged_collapsed.report.json"
    path "merged_collapsed.summary.txt"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon classify $sorted_isoform_gff $sorted_annotation $genome_fasta
    """
}

process pigeonFilter {
    publishDir "${params.output_dir}", mode: 'copy'
    
    input:
    path classification
    path junctions
    path isoform_gff

    output:
    
    path "merged_collapsed_classification.filtered_lite_classification.txt", emit: filtered_classification
    path "merged_collapsed_classification.filtered_lite_junctions.txt"
    path "merged_collapsed.sorted.filtered_lite.gff", emit: filtered_gff
    path "merged_collapsed_classification.filtered_lite_reasons.txt"
    path "merged_collapsed_classification.filtered.summary.txt"
    path "merged_collapsed_classification.filtered.report.json"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon filter $classification --isoforms $isoform_gff
    """
}

process filterByExpression {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"
    label "short_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path annotation_gtf
    path genome_fasta
    path id_to_sample
    path classification
    path read_stat
    path filtered_gff
    val min_reads
    val min_n_sample
    
    output:
    path "final_classification.parquet", emit: final_sample_classification
    path "final_transcripts.gtf", emit: final_sample_gtf
    path "final_expression.parquet", emit: final_sample_expression
    path "final_transcripts.fasta", emit: final_sample_fasta
    
    script:
    """
    filter_by_expression.py \\
        --annotation_gtf $annotation_gtf \\
        --id_to_sample $id_to_sample \\
        --classification $classification \\
        --read_stat $read_stat \\
        --filtered_gff $filtered_gff \\
        --min_reads $min_reads \\
        --min_n_sample $min_n_sample \\
        --classification_output "final_classification.parquet" \\
        --gtf_output "final_transcripts.gtf" \\
        --expression_output "final_expression.parquet"

    ~/miniforge3/envs/SQANTI3.env/bin/gffread -w final_transcripts.fasta -g $genome_fasta "final_transcripts.gtf"
    """
}

// process generatePlotsFigure1 {
//     publishDir "figures/figure_1", mode: 'copy'
//     input:

//     output:

//     script:
//     """
//     generate_figure_1.py \\


//     """
// }
process TransDecoderLongOrfs {
    label "short_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'

    input:

    path final_sample_fasta

    output:
    path "${final_sample_fasta}.transdecoder_dir/", emit: longest_orfs_dir
    path "${final_sample_fasta}.transdecoder_dir/longest_orfs.pep", emit: longest_orfs_pep
    path "${final_sample_fasta}.transdecoder_dir/longest_orfs.gff3", emit: longest_orfs_gff3

    script:
    """
    TransDecoder.LongOrfs -S -t $final_sample_fasta
    """
}

process blastpTransDecoder {
    label "short_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'

    input:

    path longest_orfs
    path protein_fasta

    output:
    path "blastp.outfmt6"

    script:
    """
    makeblastdb -dbtype prot -in $protein_fasta

    blastp -query "$longest_orfs" \\
        -db $protein_fasta -max_target_seqs 1 \\
        -outfmt 6 -evalue 1e-5 -num_threads $task.cpus > blastp.outfmt6
    """

    stub:
    """
    touch "blastp.outfmt6"
    """
}

process hmmSearch {
    label "mid_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path longest_orfs
    path hmmfile

    output:
    path "pfam.domtblout", emit : domtblout

    script:
    """
    hmmsearch --cpu $task.cpus -E 1e-10 --domtblout pfam.domtblout $hmmfile $longest_orfs
    """

    stub:
    """
    touch pfam.domtblout
    """
}

process transDecoderPredict {
    label "short_slurm_job"
    
    input:
    path longest_orfs_dir
    path final_sample_fasta
    path domtblout
    path blastpout
    
    output:
    path "${final_sample_fasta}.transdecoder.pep"
    path "${final_sample_fasta}.transdecoder.cds"
    path "${final_sample_fasta}.transdecoder.gff3", emit: transdecoder_gff3
    path "${final_sample_fasta}.transdecoder.bed"

    script:
    """
    rm -rf ${final_sample_fasta}.transdecoder_dir/__checkpoints_TDpredict
    TransDecoder.Predict --single_best_only -t ${final_sample_fasta} --retain_pfam_hits $domtblout --retain_blastp_hits $blastpout
    """

    stub:
    """
    touch "${final_sample_fasta}.transdecoder.pep"
    touch "${final_sample_fasta}.transdecoder.cds"
    touch "${final_sample_fasta}.transdecoder.gff3"
    touch "${final_sample_fasta}.transdecoder.bed"
    """
}

process cdnaAlignmentOrfToGenome {
    // label "short_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path transdecoder_gff3
    path final_sample_gtf
    path final_sample_fasta

    output:
    path "${final_sample_fasta}.transdecoder.genome.gff3", emit: genome_gff3

    script:
    """
    /usr/local/bin/util/gtf_to_alignment_gff3.pl $final_sample_gtf > transcripts.gff3

    /usr/local/bin/util/cdna_alignment_orf_to_genome_orf.pl \\
        $transdecoder_gff3 \\
        transcripts.gff3 \\
        $final_sample_fasta > "${final_sample_fasta}.transdecoder.genome.gff3"
    """

    stub:
    """
    cp -r /gpfs/fs0/scratch/s/shreejoy/nxu/SFARI/nextflow_results/transcripts_filtered.fasta.transdecoder_dir ./

    /usr/local/bin/util/gtf_to_alignment_gff3.pl $final_sample_gtf > transcripts.gff3

    /usr/local/bin/util/cdna_alignment_orf_to_genome_orf.pl \\
        /gpfs/fs0/scratch/s/shreejoy/nxu/SFARI/nextflow_results/transcripts_filtered.fasta.transdecoder.gff3 \\
        transcripts.gff3 \\
        $final_sample_fasta > "${final_sample_fasta}.transdecoder.genome.gff3"
    """
}

// GFF to GTF conversion: Stop codon not removed in general #314 
// https://github.com/NBISweden/AGAT/issues/314

process convertGenomeGff3toGtf {
    publishDir "${params.output_dir}", mode: 'copy'
    label "short_slurm_job"
    input:
    path genome_gff3

    output:
    path "transdecoder.gtf"

    script:
    """
    reformat_transdecoder_gff3.py \\
        --genome_gff3 $genome_gff3 \\
        --output "reformted_${genome_gff3}"

    apptainer run \\
        -B \${PWD},\${GENOMIC_DATA_DIR} \\
        \${NXF_SINGULARITY_CACHEDIR}/agat_1.0.0--pl5321hdfd78af_0.sif \\
        agat_sp_add_start_and_stop.pl --gff "reformted_${genome_gff3}" --fasta "\${GENOMIC_DATA_DIR}/GENCODE/GRCh38.primary_assembly.genome.fa" --out "added_codons_${genome_gff3}"

    apptainer run \\
        -B \${PWD},\${GENOMIC_DATA_DIR} \\
        \${NXF_SINGULARITY_CACHEDIR}/agat_1.0.0--pl5321hdfd78af_0.sif \\
        agat_convert_sp_gff2gtf.pl --gff "added_codons_${genome_gff3}" -o transdecoder.gtf
    """

}

// Remove the transcripts whose gene_type is not protein_coding

process runORFanage {
    publishDir "${params.output_dir}", mode: 'copy'
    label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"

    input:
    path genome_fasta
    path annotation_gtf
    path final_sample_gtf

    output:
    path "orfanage.gtf", emit: orfanage_gtf
    path "orfanage.stats"

    script:
    """
    /home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env/bin/gffread \\
        -g $genome_fasta \\
        --adj-stop \\
        -T -F -J \\
        -o corrected.gtf \\
        $final_sample_gtf

    /home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/orfanage \\
        --reference $genome_fasta \\
        --query corrected.gtf \\
        --output orfanage_without_gene_id.gtf \\
        --threads $task.cpus \\
        --minlen 70 \\
        --stats orfanage.stats \\
        $annotation_gtf

    add_geneID_to_orfanage.py \\
        --orfanage_gtf orfanage_without_gene_id.gtf \\
        --output "orfanage.gtf"
    """
}

process extractORFanageTranslationFasta {
    publishDir "${params.output_dir}", mode: 'copy'
    label "short_slurm_job"

    input:
    path genome_fasta
    path orfanage_gtf

    output:
    path "orfanage_peptide.fasta", emit: orfanage_peptide
    path "orfanage_cds.fasta", emit: orfanage_cds

    script:
    """
    agat_sp_extract_sequences.pl -g $orfanage_gtf -f $genome_fasta -t cds -p -o orfanage_peptide.fasta
    agat_sp_extract_sequences.pl -g $orfanage_gtf -f $genome_fasta -t cds -o orfanage_cds.fasta
    """
}

process getBestOrfCsv {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"

    input:
    path final_sample_classification
    path final_sample_fasta
    path orfanage_cds

    output:
    path "best_orf.tsv"

    script:
    """
    get_best_orf.py \\
        --final_sample_classification $final_sample_classification \\
        --final_sample_fasta $final_sample_fasta \\
        --cds_fasta $orfanage_cds \\
        --output "best_orf.tsv"
    """
}

process filterOrfanage {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"

    input:
    path orfanage_gtf

    output:
    path "filtered_orfanage.gtf"

    script:
    """
    filter_orfanage_res.py \\
        --orfanage_gtf $orfanage_gtf \\
        --output "filtered_orfanage.gtf"
    """
}
process renameCdsToExon {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    label "short_slurm_job"

    input:
    path predicted_cds_gtf
    path annotation_gtf

    output:
    path "SFARI.cds_renamed_exon.gtf", emit: sample_cds_renamed
    path "SFARI.transcript_exons_only.gtf", emit: sample_transcript_exon_only
    path "gencode.cds_renamed_exon.gtf", emit: ref_cds_renamed
    path "gencode.transcript_exons_only.gtf", emit: ref_transcript_exon_only

    script:
    """
    rename_cds_to_exon.py \\
        --sample_gtf $predicted_cds_gtf \\
        --sample_name "SFARI" \\
        --reference_gtf $annotation_gtf \\
        --reference_name "gencode" \\
        --num_cores $task.cpus
    """
}

process sqantiProtein {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    label "short_slurm_job"

    input:
    path sample_cds_renamed
    path sample_transcript_exon_only
    path ref_cds_renamed
    path ref_transcript_exon_only
    path best_orf
    path predicted_cds_gtf
    
    output:
    path "${predicted_cds_gtf}.sqanti_protein_classification.tsv"

    script:

    """
    apptainer run \\
        -B "/gpfs/fs0/scratch/s/shreejoy/nxu/SFARI" \\
        \${NXF_SINGULARITY_CACHEDIR}/sqanti_protein.sif \\
        sqanti3_protein.py \\
        "SFARI.transcript_exons_only.gtf" \\
        "SFARI.cds_renamed_exon.gtf" \\
        $best_orf \\
        "gencode.transcript_exons_only.gtf" \\
        "gencode.cds_renamed_exon.gtf" \\
        -d ./ \\
        -p $predicted_cds_gtf
    """

}

process fivePrimeUtr {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    label "short_slurm_job"

    input:
    path annotation_gtf
    path predicted_cds_gtf
    path protein_classification

    output:
    path "${predicted_cds_gtf}.sqanti_protein_classification_w_5utr_info.tsv"

    script:
    """
    1_get_gc_exon_and_5utr_info.py \\
        --gencode_gtf $annotation_gtf \\
        --odir ./

    2_classify_5utr_status.py \\
        --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \\
        --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \\
        --sample_cds_gtf $predicted_cds_gtf \\
        --odir ./
    
    3_merge_5utr_info_to_pclass_table.py \\
        --utr_info pb_5utr_categories.tsv \\
        --sqanti_protein_classification $protein_classification \\
        --odir ./  
    """
}

process proteinClassification {

    publishDir "${params.output_dir}", mode: 'copy'
    label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    path protein_classification_w_5utr

    output:
    path "SFARI_unfiltered.protein_classification.tsv"

    script:
    """
    apptainer run \\
        -B \${PWD},\${GENOMIC_DATA_DIR} \\
        \${NXF_SINGULARITY_CACHEDIR}/sqanti_protein.sif \\
        protein_classification.py \\
        --sqanti_protein  $protein_classification_w_5utr \\
        --name "SFARI" \\
        --dest_dir ./
    """
}

process makePacBioDatabase {
    publishDir "${params.output_dir}", mode: 'copy'

    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    val mode
    path h5ad_file_orf
    path genome_gff3
    path protein_classification_unfiltered
    path gencode_fasta

    output:
    path "*.fasta"

    script:
    if( mode == 'hybrid' )
        """
        filter_genome_gff3.py \\
            --mode $mode \\
            --h5ad_file_orf /scratch/s/shreejoy/nxu/SFARI/nextflow_results/pbid_orf.h5ad \\
            --genome_gff3 $genome_gff3 \\
            --protein_classification_unfiltered $protein_classification_unfiltered \\
            --output "${mode}_novel_transcripts.gff3"

        apptainer run \\
            -B \${PWD},\${GENOMIC_DATA_DIR} \\
            \${NXF_SINGULARITY_CACHEDIR}/agat_1.0.0--pl5321hdfd78af_0.sif \\
            agat_sp_extract_sequences.pl -g ${mode}_novel_transcripts.gff3 -f \${GENOMIC_DATA_DIR}GENCODE/GRCh38.primary_assembly.genome.fa -t cds -p -o ${mode}.fasta

        cat $gencode_fasta >> ${mode}.fasta
        """

    else if ( mode == 'pacbio' )
        """
        filter_genome_gff3.py \\
            --mode $mode \\
            --h5ad_file_orf /scratch/s/shreejoy/nxu/SFARI/nextflow_results/pbid_orf.h5ad \\
            --genome_gff3 $genome_gff3 \\
            --protein_classification_unfiltered $protein_classification_unfiltered \\
            --output "${mode}_novel_transcripts.gff3"

        apptainer run \\
            -B \${PWD},\${GENOMIC_DATA_DIR} \\
            \${NXF_SINGULARITY_CACHEDIR}/agat_1.0.0--pl5321hdfd78af_0.sif \\
            agat_sp_extract_sequences.pl -g ${mode}_novel_transcripts.gff3 -f \${GENOMIC_DATA_DIR}GENCODE/GRCh38.primary_assembly.genome.fa -t cds -p -o ${mode}.fasta
        """
}

process cometSearch {
    // label "short_slurm_job"
    input:
    path comet_params
    path search_database
    path "file*.mzXML"

    output:
    path "file*.pin"

    script:
    """
    ~/tools/comet.linux.exe -P$comet_params -D$search_database file*.mzXML
    """
}

process runPercolator {
    publishDir "${params.output_dir}", mode: 'copy'
    // label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"
    
    input:
    val mode
    path "file*.pin"

    output:
    path "${mode}_percolator.tsv"

    script:
    """
    ls file*.pin | xargs -I {} tail -n +2 {} > pooled.pin
    
    echo "\$(head -1 file1.pin)" | cat - pooled.pin > pooled2.pin

    rm pooled.pin
    
    percolator pooled2.pin > percolator.tsv

    rm pooled2.pin

    awk '{
        for (i = 1; i <= NF; i++) {
            if (i <= 5) {
                printf "%s\\t", \$i;
            } else {
                printf "%s%s", \$i, (i < NF ? "," : "");
            }
        }
        printf "\\n";
    }' OFS="\\t" percolator.tsv > percolator_res.tsv

    mv percolator_res.tsv ${mode}_percolator.tsv
    """
}

process addPeptideSupport{
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    val mode
    path h5ad_file_orf
    path percolator_tsv

    output:
    path "pbid_orf_${mode}.h5ad"

    script:
    """
    add_peptide_support.py \\
        --h5ad_file $h5ad_file_orf \\
        --peptide_file $percolator_tsv \\
        --output "pbid_orf_${mode}.h5ad"
    """

    stub:
    """
    add_peptide_support.py \\
        --h5ad_file /scratch/s/shreejoy/nxu/SFARI/nextflow_results/pbid_orf.h5ad \\
        --peptide_file $percolator_tsv \\
        --output "pbid_orf_${mode}.h5ad"
    """
}

process prepareIsoformSwitchAnalysis {
    label "short_slurm_job"
    publishDir "${params.output_dir}", mode: 'copy'
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"

    input:
    path h5ad_file
    path filtered_gff
    path final_sample_fasta

    output:
    path "IsoseqsSwitchList.rds"

    script:
    """
    IsoformSwitchAnalyzeR.py \\
        --h5ad_file $h5ad_file \\
        --filtered_gff $filtered_gff \\
        --fasta_file $final_sample_fasta \\
        --rds_file "IsoseqsSwitchList.rds"
    """
}

process peptideTrackUCSC {
    publishDir "${params.output_dir}", mode: 'copy'
    input:
    val mode
    path genome_gff3_gtf
    path reference_gtf
    path peptides
    path h5ad_file
    path sample_fasta

    output:
    path "SFARI_peptides_${mode}.gtf"

    script:
    """
    make_peptide_gtf_file.py \\
        --filter \\
        --sample_gtf $genome_gff3_gtf \\
        --reference_gtf $reference_gtf \\
        --peptides $peptides \\
        --h5ad_file $h5ad_file \\
        --sample_fasta $sample_fasta \\
        --output "SFARI_peptides_${mode}.gtf"
    """

    stub:
    """
    make_peptide_gtf_file.py \\
        --filter \\
        --sample_gtf $genome_gff3_gtf \\
        --reference_gtf $reference_gtf \\
        --peptides $peptides \\
        --h5ad_file "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/pbid_orf.h5ad" \\
        --sample_fasta $sample_fasta \\
        --output "SFARI_peptides_${mode}.gtf"
    """    
}

process peptideTrackShinyApp {
    publishDir "${params.output_dir}", mode: 'copy'
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"
    input:
    val mode
    path genome_gff3_gtf
    path unfiltered_protein_classification
    path reference_gtf
    path peptides_gtf

    output:
    path "SFARI_peptides_${mode}.csv"
    path "genome_gff3_gtf_${mode}.csv"

    script:
    """
    add-protein_classification_base-and-gene_name.py \\
        --genome_gff3_gtf $genome_gff3_gtf \\
        --unfiltered_protein_classification $unfiltered_protein_classification \\
        --reference_gtf $reference_gtf \\
        --output "genome_gff3_gtf_processed.csv"

    get_novel_peptides.R \\
        $peptides_gtf \\
        $reference_gtf \\
        $genome_gff3_gtf \\
        "genome_gff3_gtf_processed.csv" \\
        "SFARI_peptides_${mode}.csv" \\
        "genome_gff3_gtf_${mode}.csv"
    """
}

process peptideTrackAll {
    publishDir "${params.output_dir}", mode: 'copy'
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"
    input:
    val mode
    path genome_gff3_gtf
    path reference_gtf
    path peptides
    path h5ad_file
    path sample_fasta
    
    output:
    path "all_peptides_annotated_${mode}.csv"

    script:
    """
    make_peptide_gtf_file.py \\
        --sample_gtf $genome_gff3_gtf \\
        --reference_gtf $reference_gtf \\
        --peptides $peptides \\
        --h5ad_file "/scratch/s/shreejoy/nxu/SFARI/nextflow_results/pbid_orf.h5ad" \\
        --sample_fasta $sample_fasta \\
        --output "all_peptides.gtf"
    
    get_novel_peptides.R \\
        all_peptides.gtf \\
        $reference_gtf \\
        $genome_gff3_gtf \\
        "genome_gff3_gtf_processed.csv" \\
        "all_peptides_annotated_${mode}.csv"
    """
}

workflow {
    isoform_gff = Channel.fromPath("proc/merged_collapsed.gff")
    getIDToSample = channel.fromPath("proc/id_to_sample.txt")
    read_stat = channel.fromPath("proc/merged_collapsed.read_stat.txt")
    pigeonPrepare(isoform_gff, params.annotation_gtf, params.genome_fasta)
    pigeonClassify(pigeonPrepare.out.sorted_isoform_gff, pigeonPrepare.out.sorted_isoform_gff_pgi, pigeonPrepare.out.sorted_annotation, pigeonPrepare.out.sorted_annotation_gtf_pgi, params.genome_fasta, pigeonPrepare.out.reference_fasta_pgi)
    pigeonFilter(pigeonClassify.out.classification, pigeonClassify.out.junctions, pigeonPrepare.out.sorted_isoform_gff)
    filterByExpression(params.annotation_gtf, params.genome_fasta, getIDToSample, pigeonFilter.out.filtered_classification, read_stat, pigeonFilter.out.filtered_gff, params.min_reads, params.min_n_sample)
    runORFanage(params.genome_fasta, params.annotation_gtf, filterByExpression.out.final_sample_gtf)
    extractORFanageTranslationFasta(params.genome_fasta, runORFanage.out.orfanage_gtf)
    getBestOrfCsv(filterByExpression.out.final_sample_classification, filterByExpression.out.final_sample_fasta, extractORFanageTranslationFasta.out.orfanage_cds)
    filterOrfanage(runORFanage.out.orfanage_gtf)
    renameCdsToExon(filterOrfanage.out, params.annotation_gtf)
    sqantiProtein(renameCdsToExon.out.sample_cds_renamed, renameCdsToExon.out.sample_transcript_exon_only, renameCdsToExon.out.ref_cds_renamed, renameCdsToExon.out.ref_transcript_exon_only, getBestOrfCsv.out, filterOrfanage.out)
    fivePrimeUtr(params.annotation_gtf, filterOrfanage.out, sqantiProtein.out)
    proteinClassification(fivePrimeUtr.out)
    // makePacBioDatabase(params.searchDB, addORFPredictions.out.h5ad_file_orf, cdnaAlignmentOrfToGenome.out.genome_gff3, proteinClassification.out, params.translation_fasta)
    // Channel.fromPath(params.datadir + "tc-1154/*.mzXML").collect().set { mzXMLiles }
    // cometSearch(params.comet_params, makePacBioDatabase.out, mzXMLiles)
    // runPercolator(params.searchDB, cometSearch.out)
    // addPeptideSupport(params.searchDB, addORFPredictions.out, runPercolator.out)
    // peptideTrackUCSC(params.searchDB, convertGenomeGff3toGtf.out, params.annotation_gtf, runPercolator.out, addORFPredictions.out, makePacBioDatabase.out)
    // peptideTrackShinyApp(params.searchDB, convertGenomeGff3toGtf.out, proteinClassification.out, params.annotation_gtf, peptideTrackUCSC.out)
    // peptideTrackAll(params.searchDB, convertGenomeGff3toGtf.out, params.annotation_gtf, runPercolator.out, addORFPredictions.out, makePacBioDatabase.out)
    // prepareIsoformSwitchAnalysis(getSingleCellObject.out, pigeonFilter.out.filtered_gff, filterSampleFasta.out)
    // Channel.fromPath(params.datadir + "/*/outputs/flnc.report.csv").collect().set { flncReports }
    // collectPolyATailLength(flncReports, isoseqCollapse.out.read_stat)
}