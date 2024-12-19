#!/usr/bin/env nextflow

// listOfFiles = file('data/long_read/LUO26876.20240514/*/outputs/flnc.bam')

// println listOfFiles

params.datadir = "/scratch/s/shreejoy/nxu/SFARI/data/"
params.annotation_gtf = "/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.annotation.gtf"
params.transcripts_fasta = "/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.transcripts.fa"
params.reference_fasta = "/project/s/shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa"
params.gencode_translation_fasta = "/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.pc_translations.fa"
params.comet_params = "/scratch/s/shreejoy/nxu/SFARI/data/comet.params"

params.human_hexamer = "$projectDir/data/Human_Hexamer.tsv"
params.human_logit_model = "$projectDir/data/Human_logitModel.RData"
params.protein_fasta = "$projectDir/data/uniprotkb_proteome_UP000005640_AND_revi_2024_10_07.fasta"
params.sq_out = "/scratch/s/shreejoy/nxu/SFARI/proc/SQANTI3_qc_classification.txt"
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

    stub:
    """
    touch merged_collapsed.filtered_lite_classification.txt
    touch merged_collapsed.filtered_lite_junctions.txt
    touch merged_collapsed.filtered_lite_reasons.txt
    touch merged_collapsed.sorted.filtered_lite.gff
    """
}

process getSingleCellObject {
    publishDir "nextflow_results", mode: 'copy'
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

process filterSampleFasta {
    publishDir "nextflow_results", mode: 'copy'
    input:
    path h5ad_file
    path filtered_gff
    path reference_fasta

    output:
    path "transcripts_filtered.fasta"

    script:
    """
    ~/miniforge3/envs/SQANTI3.env/bin/gffread -w transcripts.fasta -g $reference_fasta $filtered_gff

    get_pbid_list.py \\
        --h5ad_file $h5ad_file \\
        --output pbid_list.txt

    ~/miniforge3/envs/patch_seq_spl/bin/seqkit grep -f pbid_list.txt transcripts.fasta > transcripts_filtered.fasta
    """
    
    stub:
    """
    touch "transcripts_filtered.fasta"

    """
}

process TransDecoderLongOrfs {
    label "short_slurm_job"
    publishDir "nextflow_results", mode: 'copy'

    input:

    path filtered_sample_fasta

    output:
    path "${filtered_sample_fasta}.transdecoder_dir/", emit: longest_orfs_dir
    path "${filtered_sample_fasta}.transdecoder_dir/longest_orfs.pep", emit: longest_orfs_pep
    path "${filtered_sample_fasta}.transdecoder_dir/longest_orfs.gff3", emit: longest_orfs_gff3

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
    touch blastp.outfmt6
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
    publishDir "nextflow_results", mode: 'copy'
    
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
    rm -rf "${filtered_sample_fasta}.transdecoder_dir/__checkpoints_TDpredict"
    TransDecoder.Predict --single_best_only -t $filtered_sample_fasta --retain_pfam_hits $$domtblout --retain_blastp_hits $blastpout
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
    label "short_slurm_job"
    publishDir "nextflow_results", mode: 'copy'

    input:
    path transdecoder_gff3
    path filtered_sample_gtf
    path filtered_sample_fasta

    output:
    path "${filtered_sample_fasta}.transdecoder.genome.gff3"

    script:
    """
    /usr/local/bin/util/gtf_to_alignment_gff3.pl $filtered_sample_gtf > transcripts.gff3

    /usr/local/bin/util/cdna_alignment_orf_to_genome_orf.pl \\
        $transdecoder_gff3 \\
        transcripts.gff3 \\
        $filtered_sample_fasta > "${filtered_sample_fasta}.transdecoder.genome.gff3"
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

process addORFPredictions {

    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    path h5ad_object
    path longest_orfs_gff3
    path transdecoder_gff3

    output:
    path "pbid_orf.h5ad"

    script:
    """
    addORFPredictions.py \\
        --h5ad_object $h5ad_object \\
        --longest_orfs_gff3 $longest_orfs_gff3 \\
        --pred_orfs_gff3 $transdecoder_gff3 \\
        --output "pbid_orf.h5ad"
    """
}

process makePacBioDatabase {
    conda "/home/s/shreejoy/nxu/miniforge3/envs/SQANTI3.env"

    input:
    path h5ad_file_orf
    path gencode_fasta
    path filtered_sample_fasta

    output:
    path "PacBio_db.fasta"

    script:
    """
    get_novel_transcripts_txt.py \\
        --h5ad_file $h5ad_file_orf \\
        --output novel_transcripts.txt

    !~/miniforge3/envs/patch_seq_spl/bin/seqkit grep -f novel_transcripts.txt $filtered_sample_fasta > PacBio_db.fasta

    cat $gencode_fasta >> PacBio_db.fasta
    """
}

process cometSearch {
    label "short_slurm_job"
    input:
    path comet_params
    path search_database
    path "file*.mzXML"

    output:
    path "file*.pin"

    script:
    """
    ~/tools/comet.linux.exe -P $comet_params -D $search_database "file*.mzXML"
    """
}

process runPercolator {
    publishDir "nextflow_results", mode: 'copy'
    label "short_slurm_job"
    conda "/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl"
    
    input:
    path "file*.pin"

    output:
    path "percolator.tsv"

    script:
    """
    ls "file*.pin" | xargs -I {} tail -n +2 {} > pooled.pin
    
    echo "\$(head -1 file1.pin)" | cat - pooled.pin > temp && mv -f temp pooled.pin
    
    percolator pooled.pin > percolator.tsv
    """
}

workflow {
    getIDToSample(params.datadir + "long_read/LUO26876.20240514/*/outputs/flnc.bam")
    Channel.fromPath(params.datadir + "long_read/LUO26876.20240514/*/outputs/mapped.bam").collect().set { bamFiles }
    mergeBamFiles(bamFiles)
    isoseqCollapse(mergeBamFiles.out)
    pigeonPrepare(isoseqCollapse.out.isoform_gff, params.annotation_gtf, params.reference_fasta)
    pigeonClassify(pigeonPrepare.out.sorted_isoform_gff, pigeonPrepare.out.sorted_isoform_gff_pgi, pigeonPrepare.out.sorted_annotation, pigeonPrepare.out.sorted_annotation_gtf_pgi, params.reference_fasta, pigeonPrepare.out.reference_fasta_pgi)
    pigeonFilter(pigeonClassify.out.classification, pigeonClassify.out.junctions, pigeonPrepare.out.sorted_isoform_gff)
    getSingleCellObject(getIDToSample.out, pigeonFilter.out.filtered_classification, isoseqCollapse.out.read_stat)
    filterSampleFasta(getSingleCellObject.out, pigeonFilter.out.filtered_gff, params.reference_fasta)
    TransDecoderLongOrfs(filterSampleFasta.out)
    blastpTransDecoder(TransDecoderLongOrfs.out.longest_orfs_pep, params.protein_fasta)
    hmmSearch(TransDecoderLongOrfs.out.longest_orfs_pep, params.hmmfile)
    transDecoderPredict(TransDecoderLongOrfs.out.longest_orfs_dir, filterSampleFasta.out, hmmSearch.out, blastpTransDecoder.out)
    cdnaAlignmentOrfToGenome(transDecoderPredict.out.transdecoder_gff3, pigeonFilter.out.filtered_gff, filterSampleFasta.out)
    Channel.fromPath(params.datadir + "data/tc-1154/*.mzXML").collect().set { mzXMLiles }
    cometSearch(params.comet_params, mzXMLiles)
    runPercolator(cometSearch.out)
    // Channel.fromPath(params.datadir + "/*/outputs/flnc.report.csv").collect().set { flncReports }
    // collectPolyATailLength(flncReports, isoseqCollapse.out.read_stat)
}