library(ggtranscript)
library(rtracklayer)
library(readr)
library(tidyr)
library(stringr)
library(Biostrings)
library(arrow)
library(biomaRt)
library(ggtranscript)
library(ggplot2)
library(readxl)
library(GenomicFeatures)
library(txdbmaker)
library(dplyr)

# Create example figures that illustrate the usage of an alternative reference for interpreting variant effect

orfanage_gr <- import("nextflow_results/V47/orfanage/orfanage.gtf") %>% 
    as_tibble() %>% 
    filter(type %in% c("CDS", "exon")) %>% 
    select(c(seqnames, start, end, strand, transcript_id, type))

gencode_gr <- import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf")) %>% 
    as_tibble() %>% 
    filter(type %in% c("CDS", "exon")) %>% 
    select(c(seqnames, start, end, strand, transcript_id, type))

orfanage_exon_gr <- makeTxDbFromGFF("nextflow_results/V47/orfanage/orfanage.gtf", format="gtf") %>% 
    exonsBy(by=c("tx"), use.names=TRUE)

combined_gr <- bind_rows(orfanage_gr, gencode_gr)

ensembl_canonical_isoform <- read_tsv(
        paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf"),
        comment="#", 
        col_names = c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
    ) %>% 
    mutate(
        transcript_id = str_extract(attributes, "transcript_id \"[^\"]+\"") %>% str_remove_all("transcript_id |\""),
        gene_name = str_extract(attributes, "gene_name \"[^\"]+\"") %>% str_remove_all("gene_name |\""),
        canonical = str_detect(attributes, "Ensembl_canonical")
    ) %>% 
    filter(type=="transcript", canonical) %>%
    distinct(transcript_id, .keep_all=TRUE) %>% 
    select(gene_name, transcript_id)

novel_isoforms <- read_parquet("nextflow_results/V47/final_classification.parquet") %>% 
    filter(structural_category %in% c("novel_in_catalog", "novel_not_in_catalog")) %>% 
    pull(isoform)

peptides <- readAAStringSet("nextflow_results/V47/orfanage/orfanage_peptide.fasta")

isoforms_of_interest <- tibble(
        peptide_id = names(peptides),
        length = width(peptides)
    ) %>% 
    mutate(
        peptide_id = str_extract(peptide_id, "PB\\.\\d+\\.\\d+")
    ) %>% 
    filter(length >= 100) %>% 
    filter(peptide_id %in% novel_isoforms) %>% 
    pull(peptide_id)

all_variants <- read_tsv("data/All_variants_used_in_project.tsv") %>% 
    mutate(
        ExonicFunc.gencode.v47 = ifelse(ExonicFunc.gencode.v47==".", Func.gencode.v47, ExonicFunc.gencode.v47),
        ExonicFunc.orfanage_MAY_2025 = ifelse(ExonicFunc.orfanage_MAY_2025==".", Func.orfanage_MAY_2025, ExonicFunc.orfanage_MAY_2025)
    )

variants_of_interest <- all_variants %>% 
    filter(
        ExonicFunc.orfanage_MAY_2025 %in% c("stopgain", "nonsynonymous SNV", "splicing"),
        !(ExonicFunc.gencode.v47 %in% c("stopgain", "nonsynonymous SNV", "splicing"))
    ) %>% 
    separate_longer_delim(AAChange.orfanage_MAY_2025, delim = ";") %>% 
    mutate(
        isoform = str_extract(AAChange.orfanage_MAY_2025, "PB\\.\\d+\\.\\d+")
    ) %>% 
    filter(
        isoform %in% isoforms_of_interest,
        variant_type == "denovo",
        ExonicFunc.gencode.v47 != "nonsynonymous_SNV",
        Affection == 2
    ) %>%
    mutate(
        ensembl_id = str_extract(Gene.gencode.v47, "ENSG\\d+\\.\\d+")
    ) %>%
    distinct(start, ALT, .keep_all=TRUE)

mapping <- import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "GENCODE/gencode.v47.annotation.gtf")) %>% 
    as_tibble() %>% 
    filter(type == "gene") %>% 
    select(gene_id, gene_name)

variants_of_interest <- variants_of_interest %>% 
    left_join(mapping, by = c("ensembl_id" = "gene_id")) %>% 
    mutate(
        gene_name = ifelse(is.na(gene_name), Gene.gencode.v47, gene_name)
    )

DTE <- read_excel("data/DTE_DTU_Nuo.xlsx", sheet = "DTE")
DTU <- read_excel("data/DTE_DTU_Nuo.xlsx", sheet = "DTU")

SFARI_gene <- read_csv("data/SFARI-Gene_genes_07-08-2025release_09-09-2025export.csv")
ASiD_gene <- read_excel("data/PC_genes_output.xlsx", sheet="Predictions")

# Function starts here

plot_variant <- function(variant_coord) {
    print(sprintf("Plotting variant at coordinate: %d", variant_coord))
    variant_metadata <- variants_of_interest %>% filter(start == {{variant_coord}}) %>% distinct(start, .keep_all=TRUE)

    novel_isoform_id <- variant_metadata %>% 
        pull("isoform")

    novel_isoform_gr <- orfanage_exon_gr[[novel_isoform_id]]

    gene_name <- variant_metadata %>% 
        pull("gene_name")

    print(sprintf("Gene name: %s", gene_name))

    known_isoform_id <- ensembl_canonical_isoform %>% filter(gene_name == {{gene_name}}) %>% pull(transcript_id)

    variant_gr <- GRanges(
        seqnames = variant_metadata$CHROM,
        ranges = IRanges(start=variant_metadata$start, end=variant_metadata$start),
        strand = "*"
    )

    isoform_of_interest_exon <- combined_gr %>% 
        filter(transcript_id %in% c({{novel_isoform_id}}, {{known_isoform_id}}), type=="exon") %>% 
        mutate(
            transcript_id = factor(transcript_id, levels = c(novel_isoform_id, known_isoform_id))
        )

    isoform_of_interest_CDS <- combined_gr %>%
        filter(transcript_id %in% c(novel_isoform_id, known_isoform_id), type=="CDS") %>% 
        mutate(
            transcript_id = factor(transcript_id, levels = c(novel_isoform_id, known_isoform_id))
        )

    title <- sprintf("%s\n%s to %s", gene_name, pull(variant_metadata, "ExonicFunc.gencode.v47"), pull(variant_metadata, "ExonicFunc.orfanage_MAY_2025"))

    novel_isoform_gr <- orfanage_exon_gr[[novel_isoform_id]]
    
    if (as.vector(runValue(strand(novel_isoform_gr)[1])) == "-") {
        last_exon_exon_junction_tx_coord <- start(novel_isoform_gr[order(start(novel_isoform_gr))][2]) %>% 
            GRanges(seqnames=seqnames(novel_isoform_gr)[1], ranges = .,strand=strand(novel_isoform_gr)[1]) %>%
            mapToTranscripts(orfanage_exon_gr, intronJunctions=TRUE, ignore.strand=TRUE) %>% 
            subset(seqnames(.) == novel_isoform_id) %>% 
            end()
    } else if (as.vector(runValue(strand(novel_isoform_gr)[1])) == "+") {
        last_exon_exon_junction_tx_coord <- end(novel_isoform_gr[order(start(novel_isoform_gr), decreasing=TRUE)][2]) %>% 
            GRanges(seqnames=seqnames(novel_isoform_gr)[1], ranges = .,strand=strand(novel_isoform_gr)[1]) %>%
            mapToTranscripts(orfanage_exon_gr, intronJunctions=TRUE, ignore.strand=TRUE) %>% 
            subset(seqnames(.) == novel_isoform_id) %>% 
            end()
    }

    variant_tx_coord <- variant_gr %>%
        mapToTranscripts(orfanage_exon_gr, intronJunctions=TRUE, ignore.strand=TRUE) %>% 
        subset(seqnames(.) == novel_isoform_id) %>% start()

    distance_to_junction <- abs(variant_tx_coord - last_exon_exon_junction_tx_coord)
    
    print(sprintf("Distance to last exon-exon junction: %d bp", distance_to_junction))

    if (subjectHits((findOverlaps(variant_gr, novel_isoform_gr))) == length(novel_isoform_gr)) {
        title <- paste0(title, "\nVariant in terminal exon")
    } else {
        title <- paste0(title, sprintf("\nDistance to last exon-exon junction: %d bp", distance_to_junction))
    }
    
    if (novel_isoform_id %in% filter(DTE, (padj < 0.05)&(condition=="t00 vs t04"))$pb_id) {
        title <- paste0(title, "\nDTE: t00 vs t04")
    } else if (novel_isoform_id %in% filter(DTE, (padj < 0.05)&(condition=="t00 vs t30"))$pb_id) {
        title <- paste0(title, "\nDTE: t00 vs t30")
    } else if (novel_isoform_id %in% filter(DTE, (padj < 0.05)&(condition=="t04 vs t30"))$pb_id) {
        title <- paste0(title, "\nDTE: t04 vs t30")
    } else if (novel_isoform_id %in% DTE$pb_id) {
        title <- paste0(title, "\nno DTE")
    } else {
        title <- paste0(title, "\nnot passing prefilter for DTE")
    }

    if (novel_isoform_id %in% filter(DTU, (DTU_qval < 0.05)&(condition_1=="t00")&(condition_2=="t04"))$isoform_id) {
        title <- paste0(title, ", DTU: t00 vs t04")
    } else if (novel_isoform_id %in% filter(DTU, (DTU_qval < 0.05)&(condition_1=="t00")&(condition_2=="t30"))$isoform_id) {
        title <- paste0(title, ", DTU: t00 vs t30")
    } else if (novel_isoform_id %in% filter(DTU, (DTU_qval < 0.05)&(condition_1=="t04")&(condition_2=="t30"))$isoform_id) {
        title <- paste0(title, ", DTU: t04 vs t30")
    } else if (novel_isoform_id %in% DTU$isoform_id) {
        title <- paste0(title, ", no DTU")
    } else {
        title <- paste0(title, ", not passing prefilter for DTU")
    }

    if (gene_name %in% filter(SFARI_gene, `gene-score`==1)$`gene-symbol`) {
        title <- paste0(title, "\nSFARI gene score 1")
    } else if (gene_name %in% filter(SFARI_gene, `gene-score`==2)$`gene-symbol`) {
        title <- paste0(title, "\nSFARI gene score 2")
    } else if (gene_name %in% filter(SFARI_gene, `gene-score`==3)$`gene-symbol`) {
        title <- paste0(title, "\nSFARI gene score 3")
    }

    if (gene_name %in% filter(ASiD_gene, `Autism susceptibility`==1)$`Gene name`) {
        title <- paste0(title, "\nASiD gene")
    }

    isoform_of_interest_exon %>% 
        ggplot(aes(
            xstart = start,
            xend = end,
            y = transcript_id
        )) +
        geom_range(
            fill = "white",
            height = 0.25
        ) +
        geom_intron(
            data = to_intron(isoform_of_interest_exon, "transcript_id"),
            aes(strand = strand)
        ) +
        geom_range(
            data = isoform_of_interest_CDS
        ) +
        geom_intron(
            data = to_intron(isoform_of_interest_CDS, "transcript_id"),
            aes(strand = strand),
            arrow.min.intron.length = 500
        ) +
        geom_vline(xintercept = variant_coord, color="red") +
        ggtitle(title) +
        theme(
            axis.title.y = element_blank()
        )

    ggsave(str_glue("scripts/figures/figure_6/examples/{gene_name}.png"), width=9, height=3.5)
}

for (variant_coord in unique(variants_of_interest$start)) {
    tryCatch(
        {plot_variant(variant_coord)}
    )   
}

start_coord <- 50348873

variant_metadata <- variants_of_interest %>% filter(start == {{start_coord}}) %>% distinct(start, .keep_all=TRUE)

novel_isoform_id <- variant_metadata %>% 
    pull("isoform")

novel_isoform_gr <- orfanage_exon_gr[[novel_isoform_id]]

gene_name <- variant_metadata %>% 
    pull("gene_name")

known_isoform_id <- ensembl_canonical_isoform %>% filter(gene_name == {{gene_name}}) %>% pull(transcript_id)

variant_gr <- GRanges(
    seqnames = variant_metadata$CHROM,
    ranges = IRanges(start=variant_metadata$start, end=variant_metadata$start),
    strand = "*"
)

if (as.vector(runValue(strand(novel_isoform_gr)[1])) == "-") {
    last_exon_exon_junction_tx_coord <- start(novel_isoform_gr[order(start(novel_isoform_gr))][2]) %>% 
        GRanges(seqnames=seqnames(novel_isoform_gr)[1], ranges = .,strand=strand(novel_isoform_gr)[1]) %>%
        mapToTranscripts(orfanage_exon_gr, intronJunctions=TRUE, ignore.strand=TRUE) %>% 
        subset(seqnames(.) == novel_isoform_id) %>% 
        start()
} else if (as.vector(runValue(strand(novel_isoform_gr)[1])) == "+") {
    last_exon_exon_junction_tx_coord <- end(novel_isoform_gr[order(end(novel_isoform_gr), decreasing=TRUE)][2]) %>% 
        GRanges(seqnames=seqnames(novel_isoform_gr)[1], ranges = .,strand=strand(novel_isoform_gr)[1]) %>%
        mapToTranscripts(orfanage_exon_gr, intronJunctions=TRUE, ignore.strand=TRUE) %>% 
        subset(seqnames(.) == novel_isoform_id) %>% 
        end()
}

variant_tx_coord <- variant_gr %>%
    mapToTranscripts(orfanage_exon_gr, intronJunctions=TRUE, ignore.strand=TRUE) %>% 
    subset(seqnames(.) == novel_isoform_id) %>% start()

distance_to_junction <- abs(variant_tx_coord - last_exon_exon_junction_tx_coord)

subjectHits((findOverlaps(variant_gr, novel_isoform_gr))) == length(novel_isoform_gr)

DTU %>% filter(isoform_id=="PB.104608.73")