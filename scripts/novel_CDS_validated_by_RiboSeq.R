library(GenomicRanges)
library(dplyr)
library(readr)
library(rtracklayer)
library(arrow)

detected_peptides <- import("nextflow_results/V47/orfanage/UCSC_tracks/detected_peptides.gtf")
tx_classification <- read_parquet("nextflow_results/V47/final_classification.parquet")
protein_classification <- read_tsv("nextflow_results/V47/orfanage/SFARI.protein_classification.tsv")

validated_ncORFs <- read_delim("validated_ncORFs.txt", delim = "\t") %>%
    left_join(
        tx_classification %>% select(c(isoform, associated_gene)),
        join_by(transcript_id == isoform)
    ) %>%
    left_join(
        protein_classification %>% select(c(pb, tx_cat, protein_classification_base)),
        join_by(transcript_id == pb)
    )

validated_ncORFs %>% distinct(associated_gene)

validated_ncORFs %>% write_csv("export/validated_ncORFs.csv")

validated_ncORFs <- GRanges(
    seqnames = validated_ncORFs$seqname,
    ranges = IRanges(start = validated_ncORFs$start, end = validated_ncORFs$end),
    strand = validated_ncORFs$strand,
    transcript_id = validated_ncORFs$transcript_id
)

validated_ncORFs[queryHits(findOverlaps(validated_ncORFs, detected_peptides))]

# Number of unique peptides identified
length(unique(detected_peptides$transcript_id))

# Number of unique genes with detected peptides
length(unique(detected_peptides$gene_name))