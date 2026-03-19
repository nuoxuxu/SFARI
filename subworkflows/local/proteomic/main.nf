process cometSearch {
    storeDir "${params.output_dir}/${params.orf_prediction}/pin"
    label "short_slurm_job"

    input:
    path comet_params
    path search_database
    path mzXMLfiles

    output:
    path "*.pin", emit: pin_files
    path "*.mzid", emit: mzid_files

    script:
    """
    ~/tools/comet.linux.exe -P$comet_params -D$search_database *.mzXML
    """
}

process runPercolator {
    storeDir "${params.output_dir}/${params.orf_prediction}"
    label "short_slurm_job"
    
    conda "$moduleDir/environment.yml"
    
    input:
    val mode
    path pin_files

    output:
    path "${mode}_percolator.tsv"

    script:
    """
    ls *.pin | xargs -I {} tail -n +2 {} > pooled.pin
    
    files=( tc*.pin )

    first_file="\${files[0]}"

    head -n 1 "\$first_file" | cat - pooled.pin > pooled2.pin

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

process peptideTrackUCSC {
    storeDir "${params.output_dir}/${params.orf_prediction}"
    conda "/scratch/nxu/SFARI/envs/python"
    
    input:
    val mode
    path annotation_gtf
    path final_sample_classification
    path predicted_cds_gtf
    path protein_search_database
    path peptides

    output:
    path "peptides_${mode}.gtf"

    script:
    """
    make_peptide_gtf_file.py \\
        --annotation_gtf $annotation_gtf \\
        --final_sample_classification $final_sample_classification \\
        --predicted_cds_gtf $predicted_cds_gtf \\
        --protein_search_database $protein_search_database \\
        --peptides $peptides \\
        --output "peptides_${mode}.gtf"
    """
}

process peptideMapping {
    storeDir "${params.output_dir}/${params.orf_prediction}"
    
    conda "/scratch/nxu/SFARI/envs/python"

    input:
    path annotation_gtf
    path classification
    path protein_search_database
    path percolator_res
    
    output:
    path "peptide_mapping.parquet", emit: peptide_mapping
    path "novel_peptides.csv", emit: novel_peptides
    
    script:
    """
    peptide_mapping.py \\
        --annotation_gtf $annotation_gtf \\
        --final_sample_classification $classification \\
        --protein_search_database $protein_search_database \\
        --percolator_res $percolator_res
    """
}

process filterScanNumber {
    storeDir "${params.output_dir}/${params.orf_prediction}/IPSA//"
    
    conda "/scratch/nxu/SFARI/envs/python"

    input:
    path mzXML_file
    path novel_peptides
    path percolator_res

    output:
    path "*_filtered.mzML", emit: filtered_mzXML
    path "*_identifications.csv", emit: filtered_identifications
    
    script:
    """
    filter_scan_number.py \\
        --mzXML_file $mzXML_file \\
        --novel_peptides $novel_peptides \\
        --percolator_res $percolator_res
    """
}

workflow proteomic {
    take:
    final_sample_classification
    predicted_cds_gtf
    protein_search_database
    main:
    Channel.fromPath(params.mzXMLfiles).collect().set { mzXMLfiles_collected }
    cometSearch(params.comet_params, protein_search_database, mzXMLfiles_collected)
    runPercolator(params.searchDB, cometSearch.out.pin_files)
    peptideTrackUCSC(params.searchDB, params.annotation_gtf, final_sample_classification, predicted_cds_gtf, protein_search_database, runPercolator.out)
    peptideMapping(params.annotation_gtf, final_sample_classification, protein_search_database, runPercolator.out)
    mzXMLfile = channel.fromPath(params.mzXMLfiles)
    filterScanNumber(mzXMLfile, peptideMapping.out.novel_peptides.first(), runPercolator.out.first())
}

workflow {
    final_sample_classification = Channel.fromPath("nextflow_results/V47/final_classification.parquet")
    predicted_cds_gtf = Channel.fromPath("nextflow_results/V47/orfanage/orfanage.gtf")
    protein_database = Channel.fromPath("nextflow_results/V47/orfanage/hybrid.fasta")
    proteomic(final_sample_classification, predicted_cds_gtf, protein_database)
}