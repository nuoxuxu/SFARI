#!/usr/bin/env nextflow

// listOfFiles = file('data/long_read/LUO26876.20240514/*/outputs/flnc.bam')

// println listOfFiles

params.outdir = "/scratch/s/shreejoy/nxu/SFARI/data/long_read/LUO26876.20240514"
params.annotation_gtf = "/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.annotation.gtf"
params.transcripts_fasta = "/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.transcripts.fa"
params.reference_fasta = "/project/s/shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa"
params.gencode_translation_fasta = "/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.pc_translations.fa"

params.human_hexamer = "$projectDir/data/Human_Hexamer.tsv"
params.human_logit_model = "$projectDir/data/Human_logitModel.RData"
params.protein_fasta = "$projectDir/data/uniprotkb_proteome_UP000005640_AND_revi_2024_10_07.fasta"
params.sq_out = "/scratch/s/shreejoy/nxu/SFARI/proc/SQANTI3_qc_classification.txt"
params.hmmfile = "$projectDir/data/Pfam-A.hmm"
params.name = "SFARI"
params.min_junctions_after_stop_codon = 2

params.lower_kb = 1
params.upper_kb = 4
params.lower_cpm = 3

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

    stub: 
    """
    touch merged.bam
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

    stub:
    """
    touch merged_collapsed.abundance.txt
    touch merged_collapsed.fasta
    touch merged_collapsed.flnc_count.txt
    touch merged_collapsed.gff
    touch merged_collapsed.group.txt
    touch merged_collapsed.read_stat.txt
    touch merged_collapsed.report.json
    """
}

process pigeonPrepare {
    input:
    path isoform_gff
    path annotation_gtf
    path reference_fasta

    output:
    path "merged_collapsed.sorted.gff", emit: sorted_isoform_gff
    path "merged_collapsed.sorted.gff.pgi", emit: sorted_isoform_gff_pgi
    path "*.annotation.sorted.gtf", emit: sorted_annotation
    path "*.annotation.sorted.gtf.pgi", emit: sorted_annotation_gtf_pgi
    path "GRCh38.primary_assembly.genome.fa.fai", emit: reference_fasta_pgi

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon prepare $isoform_gff $annotation_gtf $reference_fasta
    """

    stub:
    """
    touch merged_collapsed.sorted.gff
    touch merged_collapsed.sorted.gff.pgi
    touch gencode.v39.annotation.sorted.gtf
    touch gencode.v39.annotation.sorted.gtf.pgi
    touch GRCh38.primary_assembly.genome.fa.fai
    """
}

process pigeonClassify {
    label "short_slurm_job"
    
    input:
    path sorted_isoform_gff
    path sorted_isoform_gff_pgi
    path sorted_annotation
    path sorted_annotation_gtf_pgi
    path reference_fasta
    path reference_fasta_pgi

    output:
    path "merged_collapsed_classification.txt", emit: classification
    path "merged_collapsed_junctions.txt", emit: junctions
    path "merged_collapsed.report.json"
    path "merged_collapsed.summary.txt"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon classify $sorted_isoform_gff $sorted_annotation $reference_fasta
    """

    stub:
    """
    touch merged_collapsed_classification.txt
    touch merged_collapsed_junctions.txt
    touch merged_collapsed.report.json
    touch merged_collapsed.summary.txt
    """
}

process pigeonFilter {
    
    input:
    path classification
    path isoform_gff

    output:
    
    path "*.filtered_lite_classification.txt", emit: filtered_classification
    path "*.filtered_lite_junctions.txt"
    path "*.filtered_lite_reasons.txt"
    path "*.filtered_lite.gff"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/pigeon filter $classification --isoforms $isoform_gff
    """

    stub:
    """
    touch merged_collapsed.filtered_lite_classification.txt
    touch merged_collapsed.filtered_lite_junctions.txt
    touch merged_collapsed.filtered_lite_reasons.txt
    touch merged_collapsed.sorted.filtered_lite.gff
    """
}

process getSingleCellObject {
    label "short_slurm_job"
    input:
    path id_to_sample
    path classification
    path read_stat

    output:
    path "pbid.h5ad"

    script:
    """
    get_SingleCell_object.py \\
        --id_to_sample $id_to_sample \\
        --classification $classification \\
        --read_stat $read_stat \\
        --output pbid.h5ad
    """

    stub:
    """
    touch pbid.h5ad
    """
}

process writePbidList {

    input:
    path filtered_h5ad_file

    output:
    path "pbid_list.txt"

    script:
    """
    #!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
    from src.single_cell import SingleCell

    lr_bulk = SingleCell("${filtered_h5ad_file}")
    lr_bulk.var["pbid"].to_frame().write_csv("pbid_list.txt", include_header=False)
    """

    stub:
    """
    touch pbid_list.txt
    """
}

process cleanSampleFasta {

    input:
    path sample_fasta

    output:
    path "cleaned.fasta"

    script:
    """
    while IFS= read -r line; do
        if [[ \$line == ">"* ]]; then
            # Extract the part before the first "|"
            echo "\${line%%|*}" >> "cleaned.fasta"
        else
            # Write the line unchanged
            echo "\$line" >> "cleaned.fasta"
        fi
    done < "$sample_fasta"
    """
    stub:
    """
    touch cleaned.fasta
    """
}

process filterSampleFasta {
    
    input:
    path pbid_list
    path cleaned_fasta

    output:
    path "cleaned_filtered.fasta"

    script:
    """
    ~/miniforge3/envs/patch_seq_spl/bin/seqkit grep -f $pbid_list $cleaned_fasta > "cleaned_filtered.fasta"
    """
    stub:
    """
    touch cleaned_filtered.fasta
    """
}

process filterSampleGTF {
    label "short_slurm_job"
    input:
    path isoform_gff
    path filtered_h5ad_file

    output:
    path "filtered_SFARI.gtf"

    script:
    """
    #!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
    import polars as pl
    from src.single_cell import SingleCell
    from src.utils import read_gtf

    lr_bulk = SingleCell('${filtered_h5ad_file}')
    gtf = read_gtf('${isoform_gff}')

    gtf\\
        .filter(pl.col("transcript_id").is_in(lr_bulk.var["pbid"]))\\
        .drop("transcript_id")\\
        .write_csv("filtered_SFARI.gtf", separator="\\t", quote_style="never", include_header=False)
    """
    stub:
    """
    touch filtered_SFARI.gtf
    """
}

process CPAT {
    label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path cleaned_filtered_fasta

    output:
    path "SFARI_cpat.output"
    path "SFARI_cpat.error"
    path "SFARI.ORF_seqs.fa", emit: orf_fasta //The top 50 ORF sequences (at least 50 nucleotides long) in FASTA format.
    path "SFARI.ORF_prob.tsv", emit: orf_coord //ORF information (strand, frame, start, end, size, Fickett TESTCODE score, Hexamer score) and coding probability)
    path "SFARI.ORF_prob.best.tsv" //The information of the best ORF. This file is a subset of "SFARI.ORF_prob.tsv"
    path "SFARI.no_ORF.txt" //Sequence IDs or BED entried with no ORF found.
    path "SFARI.r"  //Rscript file.

    script:
    """
    cpat \\
        -x $params.human_hexamer \\
        -d $params.human_logit_model \\
        -g $cleaned_filtered_fasta \\
        --min-orf=50 \\
        --top-orf=50 \\
        -o SFARI \\
        1> SFARI_cpat.output \\
        2> SFARI_cpat.error
    """
    stub:
    """
    touch "SFARI_cpat.output"
    touch "SFARI_cpat.error"
    touch SFARI.ORF_seqs.fa
    touch SFARI.ORF_prob.tsv
    touch SFARI.ORF_prob.best.tsv
    touch SFARI.no_ORF.txt
    touch SFARI.r
    """
}

process generateReferenceTables {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path annotation_gtf
    path transcripts_fasta

    output:
    path "ensg_gene.tsv", emit: ensg_gene
    path "enst_isoname.tsv", emit: enst_isoname
    path "gene_ensp.tsv"
    path "gene_isoname.tsv", emit: gene_isoname
    path "isoname_lens.tsv"
    path "gene_lens.tsv", emit: gene_lens
    path "protein_coding_genes.txt", emit: protein_coding_genes

    script:
    """
    prepare_reference_tables.py \
        --gtf $annotation_gtf \
        --fa $transcripts_fasta \
        --ensg_gene ensg_gene.tsv \
        --enst_isoname enst_isoname.tsv \
        --gene_ensp gene_ensp.tsv \
        --gene_isoname gene_isoname.tsv \
        --isoname_lens isoname_lens.tsv \
        --gene_lens gene_lens.tsv \
        --protein_coding_genes protein_coding_genes.txt
    """
}

process transcriptomeSummary {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    
    input:
    path sq_out
    path ensg_to_gene
    path enst_to_isoname
    
    output:
    path "pb_gene.tsv"

    script:
    """
    transcriptome_summary.py \\
        --sq_out $sq_out
        --ensg_to_gene $ensg_to_gene
        --enst_to_isoname $enst_to_isoname
    """
    stub:
    """
    touch pb_gene.tsv
    """
}
process getPBGene {

    input:
    path filtered_h5ad_file

    output:
    path "pb_gene.txt"

    script:
    """
    #!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python

    from src.single_cell import SingleCell

    lr_bulk = SingleCell("${filtered_h5ad_file}")

    lr_bulk.var[["pbid", "associated_gene"]]\\
        .rename({"pbid": "pb_acc", "associated_gene": "gene"})\\
        .write_csv("pb_gene.txt", separator = "\\t")
    """
    stub:
    """
    touch pb_gene.txt
    """
}

process orf_calling {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    label "short_slurm_job"

    input:
    path orf_coord
    path orf_fasta
    path annotation_gtf
    path filtered_sample_gtf
    path pb_gene
    path classification
    path sample_fasta

    output:
    path "SFARI_best_orf.tsv", emit: best_orf

    script:
    """
    orf_calling.py \\
        --orf_coord $orf_coord \\
        --orf_fasta $orf_fasta \\
        --gencode $annotation_gtf \\
        --sample_gtf $filtered_sample_gtf \\
        --pb_gene "pb_gene.txt" \\
        --classification $classification \\
        --sample_fasta $sample_fasta \\
        --num_cores $task.cpus \\
        --output SFARI_best_orf.tsv
    """
    stub:
    """
    touch SFARI_best_orf.tsv
    """
}

process refineOrfDatabase {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    
    input:
    path best_orf
    path filtered_sample_fasta

    output:
    path "SFARI_orf_refined.tsv", emit: refined_tsv
    path "SFARI_orf_refined.fasta", emit: refined_fasta
    
    script:
    """
    refine_orf_database.py \\
        --name "SFARI" \\
        --orfs $best_orf \\
        --pb_fasta $filtered_sample_fasta \\
        --coding_score_cutoff 0
    """

    stub:
    """
    touch SFARI_orf_refined.tsv
    touch SFARI_orf_refined.fasta
    """
}

process makePacbioCDSGTF {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path filtered_sample_gtf
    path agg_orfs
    path refined_orfs
    path pb_gene

    output:
    path "LRP_CDS.gtf", emit: pb_cds_gtf

    script:
    """
    make_pacbio_cds_gtf.py \\
        --sample_gtf $filtered_sample_gtf \\
        --agg_orfs /scratch/s/shreejoy/nxu/SFARI/results/long_read/SFARI_orf_refined.tsv \\
        --refined_orfs $refined_orfs \\
        --pb_gene $pb_gene \\
        --output_cds "LRP_CDS.gtf"
    """
    stub:
    """
    touch LRP_CDS.gtf
    """
}

process TransDecoderLongOrfs {
    label "short_slurm_job"
    input:

    path filtered_sample_fasta
    path protein_fasta

    output:
    path "${filtered_sample_fasta}.transdecoder_dir/", emit: longest_orfs_dir
    path "${filtered_sample_fasta}.transdecoder_dir/longest_orfs.pep", emit: longest_orfs_pep

    script:
    """
    TransDecoder.LongOrfs -S -t $filtered_sample_fasta
    """
    stub:
    """
    mkdir -p ${filtered_sample_fasta}.transdecoder_dir/
    touch ${filtered_sample_fasta}.transdecoder_dir/longest_orfs.pep
    """
}

process blastpTransDecoder {
    label "short_slurm_job"

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
        -outfmt 6 -evalue 1e-5 -num_threads 40 > blastp.outfmt6
    """
    stub:
    """
    touch blastp.outfmt6
    """
}

process hmmSearch {
    label "short_slurm_job"

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
    path filtered_sample_fasta
    path domtblout
    path blastpout
    
    output:
    path "${filtered_sample_fasta}.transdecoder.pep"
    path "${filtered_sample_fasta}.transdecoder.cds"
    path "${filtered_sample_fasta}.transdecoder.gff3", emit: transdecoder_gff3
    path "${filtered_sample_fasta}.transdecoder.bed"

    script:
    """
    TransDecoder.Predict --single_best_only -t $filtered_sample_fasta --retain_pfam_hits $domtblout --retain_blastp_hits $blastpout
    """

    stub:
    """
    touch ${filtered_sample_fasta}.transdecoder.pep
    touch ${filtered_sample_fasta}.transdecoder.cds
    touch ${filtered_sample_fasta}.transdecoder.gff3
    touch ${filtered_sample_fasta}.transdecoder.bed
    """
}

process cdnaAlignmentOrfToGenome {
    input:
    path transdecoder_gff3
    path filtered_sample_gtf
    path filtered_sample_fasta

    output:
    path "SFARI.transdecoder.genome.gff3"

    script:
    """
    /usr/local/bin/util/gtf_to_alignment_gff3.pl $filtered_sample_gtf > transcripts.gff3

    /usr/local/bin/util/cdna_alignment_orf_to_genome_orf.pl \\
        $transdecoder_gff3 \\
        transcripts.gff3 \\
        $filtered_sample_fasta > "SFARI.transdecoder.genome.gff3"
    """

    // stub:
    // """
    // /usr/local/bin/util/gtf_to_alignment_gff3.pl $filtered_sample_gtf > transcripts.gff3

    // /usr/local/bin/util/cdna_alignment_orf_to_genome_orf.pl \\
    //     "/gpfs/fs0/scratch/s/shreejoy/nxu/SFARI/proc/cleaned_filtered.fasta.transdecoder.gff3" \\
    //     transcripts.gff3 \\
    //     $filtered_sample_fasta > "SFARI.transdecoder.genome.gff3"
    // """
    stub:
    """
    touch SFARI.transdecoder.genome.gff3
    """
}

process reformatTransDecoderGFF3 {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"

    input:
    path transdecoder_genome_gff3

    output:
    path "formatted_${transdecoder_genome_gff3}"

    script:
    """
    reformat_transdecoder_gff3.py $transdecoder_genome_gff3 formatted_${transdecoder_genome_gff3}
    """
    stub:
    """
    touch formatted_${transdecoder_genome_gff3}
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

process renameCDSToExon {
    label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    path filtered_sample_gtf
    path reference_gtf

    output:
    path "SFARI.transcript_exons_only.gtf", emit: sample_transcript_exon_only
    path "SFARI.cds_renamed_exon.gtf", emit: sample_cds_renamed
    path "gencode.transcript_exons_only.gtf", emit: ref_transcript_exon_only
    path "gencode.cds_renamed_exon.gtf", emit: ref_cds_renamed
    
    script:
    """
    rename_cds_to_exon.py \\
        --sample_gtf $filtered_sample_gtf \\
        --sample_name SFARI \\
        --reference_gtf $reference_gtf \\
        --reference_name gencode \\
        --num_cores $task.cpus
    """
    stub:
    """
    touch SFARI.transcript_exons_only.gtf
    touch SFARI.cds_renamed_exon.gtf
    touch gencode.transcript_exons_only.gtf
    touch gencode.cds_renamed_exon.gtf
    """
}

process SQANTIProtein {

    input:
    path sample_exon
    path sample_cds
    path reference_exon
    path reference_cds
    path best_orf

    output:
    path "${params.name}.sqanti_protein_classification.tsv", emit: sqanti_protein_classification

    script:
    """
    sqanti3_protein.py \\
        $sample_exon \\
        $sample_cds \\
        $best_orf \\
        $reference_exon \\
        $reference_cds \\
        -d ./ \\
        -p ${params.name}
    """
    stub:
    """
    touch ${params.name}.sqanti_protein_classification.tsv
    """
}

process fivePrimeUtr {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    path reference_gtf
    path pb_cds_gtf
    path sqanti_protein_classification

    output:
    path "${params.name}.sqanti_protein_classification_w_5utr_info.tsv", emit: sqanti_protein_classification_w_5utr

    script:
    """
    1_get_gc_exon_and_5utr_info.py \\
        --gencode_gtf $reference_gtf \\
        --odir ./

    2_classify_5utr_status.py \\
        --gencode_exons_bed gencode_exons_for_cds_containing_ensts.bed \\
        --gencode_exons_chain gc_exon_chain_strings_for_cds_containing_transcripts.tsv \\
        --sample_cds_gtf $pb_cds_gtf \\
        --odir ./ 
    
    3_merge_5utr_info_to_pclass_table.py \\
        --name ${params.name} \\
        --utr_info pb_5utr_categories.tsv \\
        --sqanti_protein_classification $sqanti_protein_classification \\
        --odir ./
    """
    stub:
    """
    touch ${params.name}.sqanti_protein_classification_w_5utr_info.tsv
    """
}

process proteinClassification {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path protein_classification
    path best_orf
    path refined_info
    path ensg_gene

    output:
    path "${params.name}_unfiltered.protein_classification.tsv", emit: protein_classification_unfiltered
    path "${params.name}_genes.tsv", emit: pr_genes

    script:
    """
    protein_classification_add_meta.py \\
        --protein_classification $protein_classification \\
        --best_orf $best_orf \\
        --refined_meta $refined_info \\
        --ensg_gene $ensg_gene \\
        --name ${params.name} \\
        --dest_dir ./


    protein_classification.py \\
        --sqanti_protein ${params.name}.protein_classification_w_meta.tsv \\
        --name ${params.name}_unfiltered \\
        --dest_dir ./
    """
    stub:
    """
    protein_classification_add_meta.py \\
        --protein_classification /scratch/s/shreejoy/nxu/SFARI/nextflow_results/SFARI.sqanti_protein_classification_w_5utr_info.tsv \\
        --best_orf /scratch/s/shreejoy/nxu/SFARI/nextflow_results/SFARI_best_orf.tsv \\
        --refined_meta /scratch/s/shreejoy/nxu/SFARI/nextflow_results/SFARI_orf_refined.tsv \\
        --ensg_gene $ensg_gene \\
        --name ${params.name} \\
        --dest_dir ./


    protein_classification.py \\
        --sqanti_protein ${params.name}.protein_classification_w_meta.tsv \\
        --name ${params.name}_unfiltered \\
        --dest_dir ./    
    """
}

process proteinGeneRename {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path protein_genes
    path sample_cds
    path refined_fasta
    path refined_info

    output:
    path "${params.name}_with_cds_refined.gtf", emit: pr_renamed_refined_cds
    path "${params.name}.protein_refined.fasta", emit: pr_renamed_refined_fasta
    path "${params.name}_orf_refined_gene_update.tsv", emit: pr_renamed_refined_info

    script:
    """
    protein_gene_rename.py \\
        --sample_gtf $sample_cds \\
        --sample_protein_fasta $refined_fasta \\
        --sample_refined_info $refined_info \\
        --pb_protein_genes $protein_genes \\
        --name ${params.name}
    """
    stub:
    """
    protein_gene_rename.py \\
        --sample_gtf /scratch/s/shreejoy/nxu/SFARI/nextflow_results/LRP_CDS.gtf \\
        --sample_protein_fasta /scratch/s/shreejoy/nxu/SFARI/nextflow_results/SFARI_orf_refined.fasta \\
        --sample_refined_info /scratch/s/shreejoy/nxu/SFARI/nextflow_results/SFARI_orf_refined.tsv \\
        --pb_protein_genes $protein_genes \\
        --name ${params.name}    
    """
}

process filterProtein {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path reference_gtf
    path protein_classification
    path protein_fasta
    path sample_cds

    output:
    path "${params.name}.classification_filtered.tsv"
    path "${params.name}.filtered_protein.fasta", emit: filtered_protein_fasta
    path "${params.name}_with_cds_filtered.gtf", emit: filtered_cds

    script:
    """
    protein_filter.py \\
        --protein_classification $protein_classification \\
        --gencode_gtf $reference_gtf \\
        --protein_fasta $protein_fasta \\
        --sample_cds_gtf $sample_cds \\
        --min_junctions_after_stop_codon ${params.min_junctions_after_stop_codon} \\
        --name ${params.name}
    """
}
process  makeGencodeDatabase {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path gencode_translation_fasta
    path protein_coding_genes
    
    output:
    path "gencode_protein.fasta", emit: gc_fasta
    path "gencode_isoname_clusters.tsv"

    script:
    """
    make_gencode_database.py \\
        --gencode_fasta $gencode_translation_fasta \\
        --protein_coding_genes $protein_coding_genes \\
        --output_fasta gencode_protein.fasta \\
        --output_cluster gencode_isoname_clusters.tsv
    """
}

// Add CPM to classification_filtered.tsv

// from src.single_cell import SingleCell
// import polars as pl
// classification_filtered = pl.read_csv("nextflow_results/SFARI.classification_filtered.tsv", separator="\t")

// lr_bulk = SingleCell("results/long_read/pbid_filtered.h5ad")

// classification_filtered\
//     .rename({"pb":"pbid"})\
//     .join(lr_bulk.CPM().to_frame().with_columns(mean = pl.mean_horizontal("*"))[["pbid", "mean"]],on="pbid",how="left")\
//     .with_columns(
//         CPM = pl.when(pl.col("mean").is_null()).then(0).otherwise(pl.col("mean"))
//     )\
//     .drop("mean")\
//     .rename({"pbid": "pb"})\
//     .write_csv("nextflow_results/SFARI.classification_filtered2.tsv", separator="\t")

process makeHybridDatabase {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path protein_classification
    path gene_lens
    path pb_fasta
    path gc_fasta
    path refined_info
    path sample_cds
    
    output:
    path "*"
    path "${params.name}_cds_high_confidence.gtf", emit: high_confidence_cds
    path "${params.name}_hybrid.fasta", emit: hybrid_fasta
    path "${params.name}_refined_high_confidence.tsv"
    script:
    """
    make_hybrid_database.py \\
        --protein_classification $protein_classification \\
        --gene_lens $gene_lens \\
        --pb_fasta $pb_fasta \\
        --gc_fasta $gc_fasta \\
        --refined_info $refined_info \\
        --pb_cds_gtf $sample_cds \\
        --name ${params.name} \\
        --lower_kb ${params.lower_kb} \\
        --upper_kb ${params.upper_kb} \\
        --lower_cpm ${params.lower_cpm} \\
    """
    stub:
    """
    make_hybrid_database.py \\
        --protein_classification /scratch/s/shreejoy/nxu/SFARI/nextflow_results/SFARI.classification_filtered2.tsv \\
        --gene_lens $gene_lens \\
        --pb_fasta $pb_fasta \\
        --gc_fasta $gc_fasta \\
        --refined_info $refined_info \\
        --pb_cds_gtf $sample_cds \\
        --name ${params.name} \\
        --lower_kb ${params.lower_kb} \\
        --upper_kb ${params.upper_kb} \\
        --lower_cpm ${params.lower_cpm}
    """
}

process massSpecRawConvert {
    label "short_slurm_job"

    input:
    path raw

    output:
    path "*.mzML"

    script:
    """
    wine msconvert $raw --filter "peakPicking true 1-"
    """
}

process metamorpheusWithGencodeDatabase {
    label "short_slurm_job"
    label "metamorpheus"
    input:
    path gc_fasta
    path mass_spec
    path toml_file

    output:
    path "toml/*"
    path "search_results/*"
    path "search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv", emit: gencode_peptides
    path "search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv", emit: gencode_protein_groups

    script:
    """
    dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
    dotnet /metamorpheus/CMD.dll -d $gc_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t toml/SearchTask.toml -v normal --mmsettings settings -o ./search_results

    mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv
    mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.Gencode.tsv
    """
}

process metamorpheusWithUniprotDatabase {
    label "short_slurm_job"
    label "metamorpheus"
    
    input:
    path uniprot_fasta
    path mass_spec
    path toml_file

    output:
    path "toml/*"
    path "search_results/*"
    path "search_results/Task1SearchTask/AllPeptides.Gencode.psmtsv"
    path "search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv"
    
    script:
    """
    dotnet /metamorpheus/CMD.dll -g -o ./toml --mmsettings ./settings
    dotnet /metamorpheus/CMD.dll -d $uniprot_fasta settings/Contaminants/MetaMorpheusContaminants.xml -s $mass_spec -t $toml_file -v normal --mmsettings settings -o ./search_results

    mv search_results/Task1SearchTask/AllPeptides.psmtsv search_results/Task1SearchTask/AllPeptides.UniProt.psmtsv
    mv search_results/Task1SearchTask/AllQuantifiedProteinGroups.tsv search_results/Task1SearchTask/AllQuantifiedProteinGroups.UniProt.tsv    
    """
}

process peptideAnalysis {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    publishDir "nextflow_results", mode: 'copy'

    input:
    path gencode_peptides
    path gene_isoname
    path refined_fasta
    path filtered_fasta
    path hybrid_fasta
    path pb_gene

    output:
    path "*"

    script:
    """
    peptide_analysis.py \\
        -gmap $gene_isoname \\
        --gencode_peptides $gencode_peptides \\
        --pb_refined_fasta $refined_fasta \\
        --pb_filtered_fasta $filtered_fasta \\
        --pb_hybrid_fasta $hybrid_fasta \\
        --pb_gene $pb_gene \\
        -odir ./    
    """
    
    stub:
    """
    peptide_analysis.py \\
        -gmap $gene_isoname \\
        --gencode_peptides $gencode_peptides \\
        --pb_refined_fasta $refined_fasta \\
        --pb_filtered_fasta $filtered_fasta \\
        --pb_hybrid_fasta $hybrid_fasta \\
        --pb_gene /scratch/s/shreejoy/nxu/SFARI/nextflow_results/pb_gene.tsv \\
        -odir ./    
    """
}

process gencodeTrackVisualization {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    input:
    path annotation_gtf
    
    output:
    path "gencode_shaded.bed12"
    path "gencode.filtered.gtf", emit: gencode_filtered_gtf

    script:
    """
    gencode_filter_protein_coding.py \\
        --reference_gtf $annotation_gtf

    gtfToGenePred gencode.filtered.gtf gencode.filtered.genePred
    genePredToBed gencode.filtered.genePred gencode.filtered.bed12

    gencode_add_rgb_to_bed.py \\
        --gencode_bed gencode.filtered.bed12 \\
        --rgb 0,0,140 \\
        --version V35
    """
}

process proteinTrackVisualization{
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"
    publishDir "nextflow_results/${params.name}/track_visualization/refined/protein", pattern: "*_refined_*"
    publishDir "nextflow_results/${params.name}/track_visualization/filtered/protein", pattern: "*_filtered_*"
    publishDir "nextflow_results/${params.name}/track_visualization/hybrid/protein", pattern: "*_hybrid_*"
  
    input:
    path refined_cds
    path filtered_cds
    path hybrid_cds
    
    output:
    path "*_shaded_*"
    
    script:
    """
    #************************************
    # Convert GTF to Bed
    #************************************
    #--------------------------
    # Refined
    #--------------------------
    gtfToGenePred $refined_cds ${params.name}_refined_cds.genePred
    genePredToBed ${params.name}_refined_cds.genePred ${params.name}_refined_cds.bed12
    #--------------------------
    # Filtered
    #--------------------------
    gtfToGenePred $filtered_cds ${params.name}_filtered_cds.genePred
    genePredToBed ${params.name}_filtered_cds.genePred ${params.name}_filtered_cds.bed12

    #--------------------------
    # Hybrid
    #--------------------------
    gtfToGenePred $hybrid_cds ${params.name}_hybrid_cds.genePred
    genePredToBed ${params.name}_hybrid_cds.genePred ${params.name}_hybrid_cds.bed12


    #************************************
    # Add RGB colors
    #************************************
    #--------------------------
    # Refined
    #--------------------------
    track_add_rgb_colors_to_bed.py \\
        --name ${params.name}_refined \\
        --bed_file ${params.name}_refined_cds.bed12

    #--------------------------
    # Filtered
    #--------------------------
    track_add_rgb_colors_to_bed.py \\
        --name ${params.name}_filtered \\
        --bed_file ${params.name}_filtered_cds.bed12

    #--------------------------
    # Hybrid
    #--------------------------
    track_add_rgb_colors_to_bed.py \\
        --name ${params.name}_hybrid \\
        --bed_file ${params.name}_hybrid_cds.bed12
    """
}

process peptideTrackVisualization{
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    publishDir "nextflow_results/${params.name}/track_visualization/refined/peptide", pattern: "*_refined_*"
    publishDir "nextflow_results/${params.name}/track_visualization/filtered/peptide", pattern: "*_filtered_*"
    publishDir "nextflow_results/${params.name}/track_visualization/hybrid/peptide", pattern: "*_hybrid_*"
  
    input:
    path sample_gtf
    path reference_gtf
    path refined_peptides
    path filtered_peptides
    path hybrid_peptides
    path pb_gene
    path gene_isoname
    path hybrid_fasta
    path refined_fasta

    output:
    path "*.gtf"
    path "*.bed12"

    script:
    """
    #************************************
    # Make peptide gtf files
    #************************************
    #--------------------------
    # Refined
    #--------------------------
    make_peptide_gtf_file.py \\
        --name ${params.name}_refined \\
        --sample_gtf $sample_gtf \\
        --reference_gtf $reference_gtf \\
        --peptides $refined_peptides \\
        --pb_gene $pb_gene \\
        --gene_isoname $gene_isoname \\
        --refined_fasta $refined_fasta

    #--------------------------
    # Filtered
    #--------------------------
    make_peptide_gtf_file.py \\
        --name ${params.name}_filtered \\
        --sample_gtf $sample_gtf \\
        --reference_gtf $reference_gtf \\
        --peptides $filtered_peptides \\
        --pb_gene $pb_gene \\
        --gene_isoname $gene_isoname \\
        --refined_fasta $refined_fasta

    #--------------------------
    # Hybrid
    #--------------------------
    make_peptide_gtf_file.py \\
        --name ${params.name}_hybrid \\
        --sample_gtf $sample_gtf \\
        --reference_gtf $reference_gtf \\
        --peptides $hybrid_peptides \\
        --pb_gene $pb_gene \\
        --gene_isoname $gene_isoname \\
        --refined_fasta $hybrid_fasta

    #************************************
    # Convert GTF to bed and add RGB
    #************************************
    #--------------------------
    # Refined
    #--------------------------
    gtfToGenePred ${params.name}_refined_peptides.gtf ${params.name}_refined_peptides.genePred
    genePredToBed ${params.name}_refined_peptides.genePred ${params.name}_refined_peptides.bed12

    # add rgb to colorize specific peptides 
    finalize_peptide_bed.py \\
        --bed ${params.name}_refined_peptides.bed12 \\
        --name ${params.name}_refined

    #--------------------------
    # Filtered
    #--------------------------
    gtfToGenePred ${params.name}_filtered_peptides.gtf ${params.name}_filtered_peptides.genePred
    genePredToBed ${params.name}_filtered_peptides.genePred ${params.name}_filtered_peptides.bed12
    # add rgb to colorize specific peptides 
    finalize_peptide_bed.py \\
        --bed ${params.name}_filtered_peptides.bed12 \\
        --name ${params.name}_filtered
    #--------------------------
    # High Confidence
    #--------------------------
    gtfToGenePred ${params.name}_hybrid_peptides.gtf ${params.name}_hybrid_peptides.genePred
    genePredToBed ${params.name}_hybrid_peptides.genePred ${params.name}_hybrid_peptides.bed12
    # add rgb to colorize specific peptides 
    finalize_peptide_bed.py \\
        --bed ${params.name}_hybrid_peptides.bed12 \\
        --name ${params.name}_hybrid
    """
}

workflow {
    getIDToSample(params.outdir + "/*/outputs/flnc.bam")
    Channel.fromPath(params.outdir + "/*/outputs/mapped.bam").collect().set { bamFiles }
    mergeBamFiles(bamFiles)
    isoseqCollapse(mergeBamFiles.out)
    pigeonPrepare(isoseqCollapse.out.isoform_gff, params.annotation_gtf, params.reference_fasta)
    pigeonClassify(pigeonPrepare.out.sorted_isoform_gff, pigeonPrepare.out.sorted_isoform_gff_pgi, pigeonPrepare.out.sorted_annotation, pigeonPrepare.out.sorted_annotation_gtf_pgi, params.reference_fasta, pigeonPrepare.out.reference_fasta_pgi)
    pigeonFilter(pigeonClassify.out.classification, isoseqCollapse.out.isoform_gff)
    getSingleCellObject(getIDToSample.out, pigeonFilter.out.filtered_classification, isoseqCollapse.out.read_stat)
    writePbidList(getSingleCellObject.out)
    cleanSampleFasta(isoseqCollapse.out.sample_fasta)
    filterSampleFasta(writePbidList.out, cleanSampleFasta.out)
    filterSampleGTF(isoseqCollapse.out.isoform_gff, getSingleCellObject.out)
    CPAT(filterSampleFasta.out)
    getPBGene(getSingleCellObject.out)
    orf_calling(
        CPAT.out.orf_coord, CPAT.out.orf_fasta, params.annotation_gtf, 
        filterSampleGTF.out, getPBGene.out, pigeonFilter.out.filtered_classification, filterSampleFasta.out)
    TransDecoderLongOrfs(filterSampleFasta.out, params.protein_fasta)
    blastpTransDecoder(TransDecoderLongOrfs.out.longest_orfs_pep, params.protein_fasta)
    hmmSearch(TransDecoderLongOrfs.out.longest_orfs_pep, params.hmmfile)
    Channel.fromPath("/scratch/s/shreejoy/nxu/SFARI/results/long_read/pfam.domtblout").set { domtblout_stub }
    transDecoderPredict(TransDecoderLongOrfs.out.longest_orfs_dir, filterSampleFasta.out, hmmSearch.out, blastpTransDecoder.out)
    cdnaAlignmentOrfToGenome(transDecoderPredict.out.transdecoder_gff3, filterSampleGTF.out, filterSampleFasta.out)
    reformatTransDecoderGFF3(cdnaAlignmentOrfToGenome.out)
    refineOrfDatabase(orf_calling.out, filterSampleFasta.out)
    makePacbioCDSGTF(filterSampleGTF.out, refineOrfDatabase.out.refined_tsv, orf_calling.out, getPBGene.out)
    generateReferenceTables(params.annotation_gtf, params.transcripts_fasta)
    transcriptomeSummary(params.sq_out, generateReferenceTables.out.ensg_gene, generateReferenceTables.out.enst_isoname)
    renameCDSToExon(makePacbioCDSGTF.out, params.annotation_gtf)
    SQANTIProtein(renameCDSToExon.out.sample_transcript_exon_only, renameCDSToExon.out.sample_cds_renamed, renameCDSToExon.out.ref_transcript_exon_only, renameCDSToExon.out.ref_cds_renamed, orf_calling.out)
    fivePrimeUtr(params.annotation_gtf, makePacbioCDSGTF.out, SQANTIProtein.out)
    proteinClassification(fivePrimeUtr.out, orf_calling.out, refineOrfDatabase.out.refined_tsv, generateReferenceTables.out.ensg_gene)
    proteinGeneRename(proteinClassification.out.pr_genes, makePacbioCDSGTF.out, refineOrfDatabase.out.refined_fasta, refineOrfDatabase.out.refined_tsv)
    filterProtein(params.annotation_gtf, proteinClassification.out.protein_classification_unfiltered, proteinGeneRename.out.pr_renamed_refined_fasta, proteinGeneRename.out.pr_renamed_refined_cds)
    makeGencodeDatabase(params.gencode_translation_fasta, generateReferenceTables.out.protein_coding_genes)
    makeHybridDatabase(proteinClassification.out.protein_classification_unfiltered, generateReferenceTables.out.gene_lens, filterProtein.out.filtered_protein_fasta, makeGencodeDatabase.out.gc_fasta, proteinGeneRename.out.pr_renamed_refined_info, filterProtein.out.filtered_cds)
    massSpecRawConvert(Channel.fromPath("/scratch/s/shreejoy/nxu/SFARI/data/tc-1154/*.raw")).collect().set { massSpecFiles }
    Channel.value("${projectDir}/assets/NO_TOML_FILE").set { toml_file }
    metamorpheusWithGencodeDatabase(makeGencodeDatabase.out.gc_fasta, massSpecFiles, toml_file)
    // metamorpheusWithUniprotDatabase(params.protein_fasta, massSpecFiles, toml_file)
    peptideAnalysis(metamorpheusWithGencodeDatabase.out.gencode_peptides, generateReferenceTables.out.gene_isoname, proteinGeneRename.out.pr_renamed_refined_fasta, filterProtein.out.filtered_protein_fasta, makeHybridDatabase.out.hybrid_fasta, transcriptomeSummary.out)
    gencodeTrackVisualization(params.annotation_gtf)
    proteinTrackVisualization(proteinGeneRename.out.pr_renamed_refined_cds, filterProtein.out.filtered_cds, makeHybridDatabase.out.high_confidence_cds)
    // Channel.fromPath(params.outdir + "/*/outputs/flnc.report.csv").collect().set { flncReports }
    // collectPolyATailLength(flncReports, isoseqCollapse.out.read_stat)
}