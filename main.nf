#!/usr/bin/env nextflow

// listOfFiles = file('data/long_read/LUO26876.20240514/*/outputs/flnc.bam')

// println listOfFiles

params.datadir = "/scratch/s/shreejoy/nxu/SFARI/data/"
params.comet_params = "/scratch/s/shreejoy/nxu/SFARI/data/comet.params"

params.human_hexamer = "$projectDir/data/Human_Hexamer.tsv"
params.human_logit_model = "$projectDir/data/Human_logitModel.RData"
params.protein_fasta = "/project/s/shreejoy/Genomic_references/UniProt/UP000005640_9606.fasta"
params.hmmfile = "$projectDir/data/Pfam-A.hmm"
params.name = "SFARI"
params.min_junctions_after_stop_codon = 2

process getIDToSample {
    
    input:
    val x

    output:
    path "id_to_sample.txt"

    script:
    """
    for file in $x; do
        if [[ -f "\$file" ]]; then
            first_field=\$(~/miniforge3/envs/SQANTI3.env/bin/samtools view "\$file" | head -1 | cut -f1)
            if [[ -n "\$first_field" ]]; then
                echo -e "\$first_field\\t\$file" >> id_to_sample.txt
            else
                echo -e "No data\\t\$file" >> id_to_sample.txt
            fi
        else
            echo -e "File not found\\t\$file" >> id_to_sample.txt
        fi
    done
    """
}

process mergeBamFiles {
    label "short_slurm_job"
    
    input:
    path "bam"

    output:
    path "merged.bam"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/samtools merge bam* -o merged.bam -@ $task.cpus
    """
}

process isoseqCollapse {
    label "short_slurm_job"
    
    input:
    path merged_bam

    output:
    path "merged_collapsed.abundance.txt"
    path "merged_collapsed.fasta", emit: sample_fasta
    path "merged_collapsed.flnc_count.txt"
    path "merged_collapsed.gff", emit: isoform_gff
    path "merged_collapsed.group.txt"
    path "merged_collapsed.read_stat.txt", emit: read_stat
    path "merged_collapsed.report.json"

    script:
    """
    ~/miniforge3/envs/isoseq/bin/isoseq collapse -j $task.cpus $merged_bam merged_collapsed.gff
    """
}

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
    publishDir "nextflow_results", mode: 'copy'
    
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
    label "short_slurm_job"
    
    input:
    path id_to_sample
    path classification
    path read_stat
    path filtered_gff
    val min_reads
    val min_n_sample
    path genome_fasta

    output:
    path "final_classification.parquet", emit: final_sample_classification
    path "final_transcripts.gtf", emit: final_sample_gtf
    path "final_transcripts.fasta", emit: final_sample_fasta
    
    script:
    """
    filter_by_expression.py \\
        --id_to_sample $id_to_sample \\
        --classification $classification \\
        --read_stat $read_stat \\
        --filtered_gff $filtered_gff \\
        --min_reads $min_reads \\
        --min_n_sample $min_n_sample \\
        --classification_output "final_classification.parquet" \\
        --gtf_output  "final_transcripts.gtf"

    
    """
}

process TransDecoderLongOrfs {
    label "short_slurm_job"
    publishDir "nextflow_results", mode: 'copy'

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
    publishDir "nextflow_results", mode: 'copy'

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
    publishDir "nextflow_results", mode: 'copy'

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
    publishDir "nextflow_results", mode: 'copy'

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
    publishDir "nextflow_results", mode: 'copy'
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

process runORFanage {
    publishDir "nextflow_results", mode: 'copy'

    input:
    path genome_fasta
    path annotation_gtf
    path final_sample_gtf

    output:
    path "orfanage.gtf", emit: orfanage_gtf

    script:
    """
    /home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env/bin/gffread \
        -g $genome_fasta \
        --adj-stop \
        -T -F -J \
        -o corrected.gtf \
        $final_sample_gtf

    /home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/orfanage \
        --reference $genome_fasta \
        --query corrected.gtf \
        --output orfanage.gtf \
        $annotation_gtf
    """
}

process extractORFanageTranslationFasta {
    publishDir "nextflow_results", mode: 'copy'

    input:
    path genome_fasta
    path annotation_gtf
    path orfanage_gtf

    output:
    path "orphanage_peptide.fasta"

    script:    """
    apptainer run \\
        -B "\${PWD}","\${GENOMIC_DATA_DIR}" \\
        "\${NXF_SINGULARITY_CACHEDIR}/agat_1.0.0--pl5321hdfd78af_0.sif" \\
        agat_sp_extract_sequences.pl -g $orfanage_gtf -f $genome_fasta -t cds -p -o orfanage.fasta
    """
}

process collectPolyATailLength {
    label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path "flnc_report"
    path read_stat

    output:
    path "output.tsv"

    script:
    """
    head -n 1 flnc_report1 > all.flnc.report.csv
    tail -q -n +2 flnc_report* >> all.flnc.report.csv
    
    Collect_A-tail_lens_v0.3.py \\
        -r $read_stat \\
        -f all.flnc.report.csv -o output.tsv
    """
}

// Running Predictions.py from the command line gives out segementation fault error
// This process adds the following columns to the SingleCell object
// at_least_one_orf, predicted_orf, CDS_genomic_start, CDS_genomic_end, base_isoform

process addORFPredictions {

    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    path h5ad_object
    path longest_orfs_gff3
    path transdecoder_gff3

    output:
    path "pbid_orf.h5ad", emit: h5ad_file_orf

    script:
    """
    addORFPredictions.py \\
        --h5ad_file $h5ad_object \\
        --longest_orfs_gff3 $longest_orfs_gff3 \\
        --pred_orfs_gff3 $transdecoder_gff3 \\
        --output "pbid_orf.h5ad"
    """

    stub:
    """
    touch "pbid_orf.h5ad"
    """
}

process proteinClassification {
    publishDir "nextflow_results", mode: 'copy'
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    path final_sample_classification
    path predicted_cds_gtf
    path annotation_gtf

    output:
    path "${predicted_cds_gtf}.sqanti_protein_classification.tsv"

    script:
    """
    get_best_orf.py \\
        --final_sample_classification $final_sample_classification \\
        --pred_orfs_gff3 $predicted_cds_gtf \\
        --output "SFARI.best_orf.tsv"

    rename_cds_to_exon.py \\
        --sample_gtf $predicted_cds_gtf \\
        --sample_name "SFARI" \\
        --reference_gtf $annotation_gtf \\
        --reference_name "gencode"

    sqanti3_protein.py \\
        "SFARI.cds_renamed_exon.gtf" \\
        "SFAR.transcript_exons_only.gtf" \\
        "SFARI.best_orf.tsv" \\
        "gencode.cds_renamed_exon.gtf" \\
        "gencode.transcript_exons_only.gtf" \\
        -d ./ \\
        -p $predicted_cds_gtf

    1_get_gc_exon_and_5utr_info.py \\
        --gencode_gtf $annotation_gtf \\
        --odir ./

    2_classify_5utr_status.py \\
        --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \\
        --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \\
        --sample_cds_gtf "SFARI.cds_renamed_exon.gtf" \\
        --odir ./ 

    
    3_merge_5utr_info_to_pclass_table.py \\
        --utr_info pb_5utr_categories.tsv \\
        --sqanti_protein_classification "${predicted_cds_gtf}.sqanti_protein_classification.tsv" \\
        --odir ./  

    protein_classification.py \\
        --sqanti_protein "${predicted_cds_gtf}.sqanti_protein_classification._w_5utr_info.tsv" \\
        --name "${predicted_cds_gtf}" \\
        --dest_dir ./
    """
}

process makePacBioDatabase {
    publishDir "nextflow_results", mode: 'copy'

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
    publishDir "nextflow_results", mode: 'copy'
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
    publishDir "nextflow_results", mode: 'copy'

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
    publishDir "nextflow_results", mode: 'copy'
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
    publishDir "nextflow_results", mode: 'copy'
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
    publishDir "nextflow_results", mode: 'copy'
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
    publishDir "nextflow_results", mode: 'copy'
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
    getIDToSample(params.datadir + "long_read/LUO26876.20240514/*/outputs/flnc.bam")
    Channel.fromPath(params.datadir + "long_read/LUO26876.20240514/*/outputs/mapped.bam").collect().set { bamFiles }
    mergeBamFiles(bamFiles)
    isoseqCollapse(mergeBamFiles.out)
    pigeonPrepare(isoseqCollapse.out.isoform_gff, params.annotation_gtf, params.genome_fasta)
    pigeonClassify(pigeonPrepare.out.sorted_isoform_gff, pigeonPrepare.out.sorted_isoform_gff_pgi, pigeonPrepare.out.sorted_annotation, pigeonPrepare.out.sorted_annotation_gtf_pgi, params.genome_fasta, pigeonPrepare.out.reference_fasta_pgi)
    pigeonFilter(pigeonClassify.out.classification, pigeonClassify.out.junctions, pigeonPrepare.out.sorted_isoform_gff)
    filterByExpression(getIDToSample.out, pigeonFilter.out.filtered_classification, isoseqCollapse.out.read_stat, pigeonFilter.out.filtered_gff, params.min_reads, params.min_n_sample, params.genome_fasta)
    TransDecoderLongOrfs(filterByExpression.out.final_sample_fasta)
    blastpTransDecoder(TransDecoderLongOrfs.out.longest_orfs_pep, params.protein_fasta)
    hmmSearch(TransDecoderLongOrfs.out.longest_orfs_pep, params.hmmfile)
    transDecoderPredict(TransDecoderLongOrfs.out.longest_orfs_dir, filterByExpression.out.final_sample_fasta, hmmSearch.out, blastpTransDecoder.out)
    cdnaAlignmentOrfToGenome(transDecoderPredict.out.transdecoder_gff3, pigeonFilter.out.filtered_gff, filterByExpression.out.final_sample_fasta)
    // addORFPredictions(getSingleCellObject.out, TransDecoderLongOrfs.out.longest_orfs_gff3, transDecoderPredict.out.transdecoder_gff3)
    convertGenomeGff3toGtf(cdnaAlignmentOrfToGenome.out.genome_gff3)
    runORFanage(params.genome_fasta, params.annotation_gtf, filterByExpression.out.final_sample_classification)
    proteinClassification(filterByExpression.out.final_sample_classification, convertGenomeGff3toGtf.out.mix(runORFanage.out), params.annotation_gtf)
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