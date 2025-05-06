#!/bin/bash
export PATH="${PATH}:/home/s/shreejoy/nxu/tools/annovar"

# Building the annovar database for Orfanage
gtfToGenePred -genePredExt data/orfanage.gtf data/annovar_orfanage_db/hg38_refGene.txt

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

#Running annovar prediction for orfanage.gtf
table_annovar.pl export/variant/novel_splice_site_ClinVar.vcf \
    data/annovar_orfanage_db/ \
    --buildver hg38 \
    --out export/variant/annovar_with_orfanage_gtf/ \
    --protocol refGene \
    --operation g \
    --nastring . \
    --polish \
    --vcfinput

#Running annovar prediction for GENCODE reference
table_annovar.pl export/variant/novel_splice_site_ClinVar.vcf \
    data/annovar_gencode_db/ \
    --buildver hg38 \
    --out export/variant/annovar_with_gencode_gtf/ \
    --protocol refGene \
    --operation g \
    --nastring . \
    --polish \
    --vcfinput