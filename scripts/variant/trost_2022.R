library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(GenomicRanges)

de_novo_variants <- read_excel("data/mmc2.xlsx", sheet="Table S2C", skip=1)
novel_exonic_regions <- readRDS("export/variant/novel_exonic_regions.rds")
# novel_splice_sites <- read_csv("export/variant/novel_splice_sites.csv")
novel_splice_sites <- GRanges(
    seqnames = novel_splice_sites$chrom,
    ranges   = IRanges(start = novel_splice_sites$start, end = novel_splice_sites$end),
    strand   = "*"
)

de_novo_variants <- de_novo_variants %>% 
    separate(
        col   = Variant,
        into  = c("chr", "pos", "ref", "alt"),
        sep   = ":",
        convert = TRUE
    ) %>% 
        filter(!grepl("_", chr)) %>%
        mutate(
            chr = paste0("chr", chr)
    ) %>% 
        filter(nchar(ref) == 1)

de_novo_variants_gr <- GRanges(
    seqnames = de_novo_variants$chr,
    ranges   = IRanges(start = de_novo_variants$pos, end = de_novo_variants$pos),
    strand   = "*"
)

de_novo_variants_in_novel_exonic_regions <- de_novo_variants[queryHits(findOverlaps(de_novo_variants_gr, novel_exonic_regions)), ] %>% 
    distinct(chr, pos, ref, alt, .keep_all = TRUE)

de_novo_variants_in_novel_exonic_regions %>% 
    # remove chr from chromosome names
    mutate(
        chr = gsub("chr", "", chr)
    ) %>%
    # insert a column identical to pos next to pos
    mutate(
        pos2 = pos
    ) %>%
    select(chr, pos, pos2, ref, alt) %>% 
    write_tsv(
        file = "export/variant/novel_exonic_regions_de_novo.avinput",
        col_names = FALSE
    )

read_tsv("export/variant/annovar/de_novo/orfanage/.hg38_multianno.txt") %>% 
    ggplot(aes(ExonicFunc.refGene)) +
    geom_bar() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

de_novo_variants_in_novel_exonic_regions %>% 
    ggplot(aes(`Variant type`)) +
    geom_bar() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
    )

# Format de novo variants for ANNOVAR

de_novo_variants %>% 
    separate(
        col   = Variant,
        into  = c("chr", "pos", "ref", "alt"),
        sep   = ":",
        convert = TRUE
    ) %>% 
    mutate(
        chr = gsub("chr", "", chr)
    ) %>%
    # insert a column identical to pos next to pos
    mutate(
        pos2 = pos
    ) %>%
    filter(nchar(ref) == 1) %>%
    select(chr, pos, pos2, ref, alt) %>% 
    write_tsv(
        file = "export/variant/trost_2022_de_novo.avinput",
        col_names = FALSE
    )

de_novo_orfanage <- read_tsv("export/variant/annovar/de_novo/orfanage/.hg38_multianno.txt")

de_novo_orfanage_gr <- GRanges(
    seqnames = de_novo_orfanage$Chr,
    ranges   = IRanges(start = de_novo_orfanage$Start, end = de_novo_orfanage$End),
    strand   = "*",
    Func.refGene = de_novo_orfanage$Func.refGene
)

seqlevelsStyle(de_novo_orfanage_gr) <- "UCSC"

de_novo_orfanage[queryHits(findOverlaps(de_novo_orfanage_gr, novel_exonic_regions)), ] %>%
    ggplot(aes(Func.refGene)) +
    geom_bar()