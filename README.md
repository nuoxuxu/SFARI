# Recreate conda environment
```bash
mamba activate
mkdir envs
mamba env create -p ./envs/tama
ln -s ${PWD}/envs/tama ${CONDA_PREFIX}/envs/tama
```

# Install external datasets

```
wget https://github.com/gandallab/Dev_Brain_IsoSeq/raw/refs/heads/main/data/cp_vz_0.75_min_7_recovery_talon.gtf.gz -O data/patowary/data/cp_vz_0.75_min_7_recovery_talon.gtf.gz
gunzip data/patowary/data/cp_vz_0.75_min_7_recovery_talon.gtf.gz

wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz -O data/hg19ToHg38.over.chain.gz


```
# Long-read sequencing
## Recommended bulk Iso-Seq workflow

| Command         | Description                                                                                  | Output format     |
| --------------- | -------------------------------------------------------------------------------------------- | ----------------- |
| lima            | Remove cDNA primers                                                                          | fl.bam            |
| isoseq refine   | Remove polyA tail and artificial concatemers                                                 | flnc.bam          |
| isoseq cluster2 | De novo isoform-level clustering scalable to large number of reads (e.g. 40-100M FLNC reads) | clustered.bam     |
| pbmm2           | Align to the genome                                                                          | mapped.bam        |
| isoseq collapse | Collapse redundant transcripts based on exonic structures                                    | collapsed.gff     |
| pigeon classify | Classify transcripts against annotation                                                      | GFF and TXT files |
| pigeon filter   | Filter transcripts for potential artifacts                                                   | GFF and TXT files |
## How to harmonize transcript names across 15 samples?
### Merge mapped transcripts using `tama merge`
This approach was described in [@zhang2022a]. 

```shell
from pathlib import Path

localrules: generate_filelist, get_fofn

rule gtf_to_bed12:
    input:
        "data/long_read/LUO26876.20240514/{sample}/outputs/isoseq_transcripts.sorted.filtered_lite.gff"
    output:
        "proc/tama_merge_outputs/{sample}.bed"
    conda: "tama"
    shell: 
        "python scripts/tama_format_gtf_to_bed12_ensembl.py {input} {output}"

SAMPLES = [path.stem for path in Path("data/long_read/LUO26876.20240514").iterdir()]

rule generate_filelist:
    input:
        ["proc/tama_merge_outputs/{sample}.bed".format(sample=sample) for sample in SAMPLES]
    output:
        "proc/tama_merge_outputs/filelist.txt"
    conda: "patch_seq_spl"
    script:
        "scripts/generate_filelist.py"
        
rule tama_merge:
    input:
        "proc/tama_merge_outputs/filelist.txt"
    conda: "tama"
    resources:
        runtime=600
    shell:
        "python scripts/tama_merge.py -f {input} -p proc/tama_merge_outputs/merged_annos"
```
The problem with this approach is that the resulting annotation only contains 6457 genes and there are 586898 transcripts which means that there are around 90 transcripts per gene. This perhaps is a result of the way I make the file list. 
## isoseq cluster2
The fact that the Kinnex full-length RNA kit only has 12 cDNA primers that mark each sample and we have 15 samples is not a problem for clustering. Because the primers have already been removed in the step `isoseq refine`, which is the step that produces `flnc.bam`. `flnc.bam` is the only input needed for `isoseq cluster2`. The identifier is in the field QNAME of `flnc.bam`. For example, if QNAME is "m84090_240409_212646_s4/262737378/ccs/41_3294", then "m84090_240409_212646_s4" uniquely identifies a read from a particular sample, and this identifier also appears in the name of output from `lima` for each sample. In conclusion, this identifier is unique for each sample and it is easy to get a mapping from the sample names (e.g., CN_1_2) to the identifier.

Clustering one sample took around 4 hours on Niagara. See [slurm_log](slurm_log.md). Therefore, clustering 15 samples would probably take 4*15=60 hours, which exceeeds the maximum wall time on Niagara. To move `flnc.bam` files to other cluster.

```shell
rule cluster2:
    input:
        "proc/flnc.fofn"
    output:
        "proc/clustered.bam"
    conda: "tama"
    resources:
        runtime=7200
        
    shell:
        "isoseq cluster2 {input} {output}"
```
Clustering one sample took around 4 hours on Niagara. See [slurm_log](slurm_log.md). 

# Alignment
## pbmm2
[pbmm2 documentation](https://github.com/PacificBiosciences/pbmm2)

# Resources
[Elizabeth Tseng's LinkedIn page](https://www.linkedin.com/in/elizabeth-tseng-phd-95153696/recent-activity/all/)
[Isoseq Docs](https://isoseq.how/)
[Lima Docs](https://lima.how/)
[tama documentation](https://github.com/GenomeRIK/tama)

# Tests
```bash
sbatch -t 0-12:0 --wrap="~/miniforge3/envs/tama/bin/python scripts/tama_merge.py -f proc/no_cap/filelist.txt -d merge_dup -p proc/no_cap/merged_annos"
```