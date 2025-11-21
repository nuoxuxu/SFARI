process cometSearch {
    publishDir "${params.output_dir}/${params.orf_prediction}/pin", mode: 'link'
    label "short_slurm_job"

    input:
    path comet_params
    path search_database
    path "file*.mzXML"

    output:
    path "file*.pin", emit: pin_files
    path "file*.mzid", emit: mzid_files

    script:
    """
    ~/tools/comet.linux.exe -P$comet_params -D$search_database file*.mzXML
    """
}

process runPercolator {
    publishDir "${params.output_dir}/${params.orf_prediction}", mode: 'link'
    label "short_slurm_job"
    
    conda "$moduleDir/environment.yml"
    
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

workflow proteomic {
    take:
    protein_database
    mzXMLfilesPath
    main:
    mzXMLfilesPath.collect().set { mzXMLfiles }
    cometSearch("/scratch/nxu/SFARI/data/comet.params", protein_database, mzXMLfiles)
    runPercolator(params.searchDB, cometSearch.out.pin_files)
    emit:
    peptides = runPercolator.out
}

workflow {
    protein_database = Channel.fromPath("nextflow_results/V47/orfanage/hybrid.fasta")
    Channel.fromPath("data/tc-1154/*.mzXML").collect().set { mzXMLfiles }
    proteomic(protein_database, mzXMLfiles)
}