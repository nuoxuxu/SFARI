#!/bin/bash
export PATH="${PATH}:${HOME}/tools/annovar"

# Building the annovar database for Orfanage
gtfToGenePred -genePredExt nextflow_results/V47/orfanage/orfanage.gtf data/annovar_orfanage_db/hg38_refGene.txt

retrieve_seq_from_fasta.pl \
   --format refGene \
   --seqfile /project/s/shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa data/annovar_orfanage_db/hg38_refGene.txt \
   --out data/annovar_orfanage_db/hg38_refGeneMrna.fa

#Building the annovar database for Gencode
gtfToGenePred -genePredExt "${GENOMIC_DATA_DIR}"/GENCODE/gencode.v47.annotation.gtf data/annovar_gencode_db/hg38_refGene.txt

retrieve_seq_from_fasta.pl \
   --format refGene \
   --seqfile /project/s/shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa data/annovar_gencode_db/hg38_refGene.txt \
   --out data/annovar_gencode_db/hg38_refGeneMrna.fa

# ClinVar with Orfanage
table_annovar.pl data/clinvar_20250421.vcf \
    data/annovar_orfanage_db/ \
    --buildver hg38 \
    --out export/variant/annovar/ClinVar/ \
    --protocol refGene \
    --operation g \
    --nastring . \
    --polish \
    --thread 40 \
    --vcfinput

# De novo variants with Orfanage
table_annovar.pl export/variant/trost_2022_de_novo.avinput \
   data/annovar_orfanage_db/ \
   --buildver hg38 \
   --out export/variant/annovar/de_novo/orfanage/ \
   --protocol refGene \
   --operation g \
   --nastring . \
   --polish

table_annovar.pl export/variant/trost_2022_de_novo.avinput \
   data/annovar_gencode_db/ \
   --buildver hg38 \
   --out export/variant/annovar/de_novo/gencode/ \
   --protocol refGene \
   --operation g \
   --nastring . \
   --polish

# Testing

retrieve_seq_from_fasta.pl \
   --format refGene \
   --seqfile /project/s/shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa test.genePred \
   --out test/hg38_refGeneMrna.fa

table_annovar.pl data/clinvar_20250421.vcf \
    test/ \
    --buildver hg38 \
    --out test_out/ \
    --protocol refGene \
    --operation g \
    --nastring . \
    --polish \
    --thread 40 \
    --vcfinput