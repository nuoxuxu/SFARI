from pathlib import Path

localrules: get_fofn, generate_filelist

tama_merge_out = ["_gene_report.txt", "_merge.txt", "_trans_report.txt", ".bed"]

rule all:
    input:
        expand("proc/{merge_priority}/merged_annos{output_file}", output_file=tama_merge_out, merge_priority=["no_cap"])

rule gtf_to_bed12:
    input:
        "data/long_read/LUO26876.20240514/{sample}/outputs/isoseq_transcripts.sorted.filtered_lite.gff"
    output:
        "proc/bed_files/{sample}.bed"
    conda: "tama"
    shell: 
        "python scripts/tama_format_gtf_to_bed12_ensembl.py {input} {output}"

SAMPLES = [path.stem for path in Path("data/long_read/LUO26876.20240514").iterdir()]

# capped or no_cap
rule generate_filelist:
    input:
        ["proc/bed_files/{sample}.bed".format(sample=sample) for sample in SAMPLES]
    output:
        "proc/{merge_priority}/filelist.txt"
    conda: "patch_seq_spl"
    script:
        "scripts/generate_filelist.py"

rule tama_merge:
    input:
        "proc/{merge_priority}/filelist.txt"
    params:
        outputdir="proc/{merge_priority}/merged_annos"
    output:
        expand("proc/{merge_priority}/merged_annos{output_file}", output_file=tama_merge_out, allow_missing=True)
    conda: "tama"
    resources:
        runtime=600
    shell:
        "python scripts/tama_merge.py -f {input} -p {params.outputdir}"

rule get_fofn:
    output:
        "proc/flnc.fofn"
    shell:
        "ls /scratch/s/shreejoy/nxu/SFARI/data/long_read/LUO26876.20240514/*/outputs/flnc.bam > {output}"
        
rule cluster2:
    input:
        "proc/flnc.fofn"
    output:
        "proc/clustered.bam"
    conda: "tama"
    shell:
        "isoseq cluster2 {input} {output}"

rule align:
    input:
        "proc/clustered.bam"
    output:
        "proc/mapped.bam"
    params:
        ref="/project/s/shreejoy/Genomic_references/Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    conda: "tama"
    shell:
        "pbmm2 align ${GENOMIC_DATA_DIR}/Ensembl/Human/Release_104/Raw/Homo_sapiens.GRCh38.dna.primary_assembly.fa {input} {output}"