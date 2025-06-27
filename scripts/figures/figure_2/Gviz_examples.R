library(dplyr)
library(readr)
library(arrow)
library(rtracklayer)
source("src/utils.R")

# Load datasets
transcript_classification <- read_parquet("nextflow_results/V47/final_classification.parquet")
peptide_mapping <- read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet")

# Novel splice-junctions
peptides_gtf <- import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>%
    as_tibble() %>% 
    filter((novelty=="novel") & (detected == "True")) %>% 
    mutate(
        start = start - 1
    ) %>% 
    filter(type.1 == "splice-junction")
peptides_gr <- GRanges(
    seqnames = peptides_gtf$seqnames,
    ranges = IRanges(start = peptides_gtf$start, end = peptides_gtf$end),
    strand = peptides_gtf$strand,
    name = peptides_gtf$transcript_id
)
peptides_SJ_gr <- split(peptides_gr, mcols(peptides_gr)$name)
peptides_SJ_gr <- sapply(peptides_SJ_gr, function(x) gaps(x)[2]) %>% 
    unname() %>% 
    do.call(c, .)

riboseq_SJ_list <- lapply(Sys.glob("data/riboseq/*.out.tab"), read_SJ_out_tab)
riboseq_SJ <- bind_rows(riboseq_SJ_list) %>% 
    mutate(
        strand = case_when(strand == 1 ~ "+",
                           strand == 2 ~ "-")
    )
riboseq_SJ_gr <- GRanges(
    seqnames = riboseq_SJ$chrom,
    ranges = IRanges(start = riboseq_SJ$start, end = riboseq_SJ$end),
    strand = mcols(riboseq_SJ_gr)$strand
) %>% unique()

findOverlaps(peptides_SJ_gr, riboseq_SJ_gr, ignore.strand=TRUE, type="equal")

# Novel exonic regions
peptides_gtf <- import("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf") %>%
    as_tibble() %>% 
    filter((novelty=="novel") & (detected == "True")) %>% 
    mutate(
        start = start - 1
    ) %>% 
    filter(type.1 == "mono-exonic")
peptides_gr <- GRanges(
    seqnames = peptides_gtf$seqnames,
    ranges = IRanges(start = peptides_gtf$start, end = peptides_gtf$end),
    strand = peptides_gtf$strand,
    name = peptides_gtf$transcript_id
)


peptides_gr[unique(queryHits(hits))]

classification <- read_tsv("nextflow_results/V47/orfanage/SFARI.protein_classification.tsv")    

SFARI_gene <- read_csv("data/SFARI-Gene_genes_01-13-2025release_03-12-2025export.csv") %>% 
    filter(`gene-score` == 1)
    
peptide_mapping <- peptide_mapping %>% 
    left_join(
        peptides_gtf %>% 
            dplyr::select(transcript_id, detected),
        join_by(peptide == transcript_id)
    ) %>% 
    left_join(
        classification %>% dplyr::select(pb, protein_classification_base), 
        join_by(pb == pb)
    ) %>% 
    left_join(
        transcript_classification %>% dplyr::select(isoform, associated_gene),
        join_by(pb == isoform)
    )

peptide_mapping %>% 
    filter(detected == "True") %>%
    filter(!GENCODE) %>%
    group_by(peptide) %>%
    filter(all(protein_classification_base %in% c("pNIC", "pNNC"))) %>%
    ungroup() %>% 
    distinct(associated_gene, .keep_all = TRUE) %>% 
    filter(associated_gene %in% SFARI_gene$`gene-symbol`) %>% 
    View()
