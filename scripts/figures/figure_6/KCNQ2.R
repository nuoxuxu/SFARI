library(ggtranscript)
library(rtracklayer)
library(readr)
library(tidyr)
library(stringr)
library(ggtranscript)
library(ggplot2)
library(readxl)
library(dplyr)

gene_of_interest <- "KCNQ2"

# Get variant info based on the gene of interest

variants_of_interest <- read_tsv("data/All_variants_used_in_project.tsv") %>% 
    mutate(
        ExonicFunc.gencode.v47 = ifelse(ExonicFunc.gencode.v47==".", Func.gencode.v47, ExonicFunc.gencode.v47),
        ExonicFunc.orfanage_MAY_2025 = ifelse(ExonicFunc.orfanage_MAY_2025==".", Func.orfanage_MAY_2025, ExonicFunc.orfanage_MAY_2025)
    ) %>% 
    filter(
        ExonicFunc.orfanage_MAY_2025 %in% c("stopgain", "nonsynonymous SNV", "splicing"),
        !(ExonicFunc.gencode.v47 %in% c("stopgain", "nonsynonymous SNV", "splicing"))
    ) %>% 
    separate_longer_delim(AAChange.orfanage_MAY_2025, delim = ";") %>% 
    mutate(
        isoform = str_extract(AAChange.orfanage_MAY_2025, "PB\\.\\d+\\.\\d+")
    ) %>% 
    filter(
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

variant_coord <- variants_of_interest %>% 
  filter(gene_name==gene_of_interest) %>% 
  pull(end)

variants_of_interest  %>% filter(gene_name==gene_of_interest) %>% pull(AAChange.gencode.v47)
variants_of_interest  %>% filter(gene_name==gene_of_interest) %>% pull(AAChange.orfanage_MAY_2025)

# Get known and novel annotation for plotting

gencode_gtf <- import(paste0(Sys.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf")) %>% 
  as_tibble()

gencode_tx <- gencode_gtf %>% 
  as_tibble() %>% 
  filter(gene_name==gene_of_interest) %>% 
  filter(type=="transcript") %>% 
  filter(tag=="GENCODE_Primary") %>%  
  slice(1) %>% 
  pull(transcript_id)  

known_annotation <- gencode_gtf %>% 
  filter(transcript_id == gencode_tx)  

novel_pb_id <- variants_of_interest %>% 
  filter(gene_name==gene_of_interest) %>% 
  pull(isoform)

novel_annotation <- import("nextflow_results/V47/orfanage/orfanage.gtf") %>% 
  as_tibble() %>% 
  filter(transcript_id == novel_pb_id)

known_annotation <- known_annotation %>% select(c(seqnames, start, end, strand, transcript_id, type))
novel_annotation <- novel_annotation %>% select(c(seqnames, start, end, strand, transcript_id, type))

combined_annotation <- bind_rows(known_annotation, novel_annotation) %>% 
  mutate(
    transcript_id = factor(transcript_id, levels = c(novel_pb_id, gencode_tx))
  )

# Plotting

my_theme <- theme_bw() +
    theme(
        axis.text.x = element_text(size = 20, angle = 45, color = "black", vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 20, color = "black"),
        strip.placement = "outside",
        legend.title = element_blank(),
        legend.position = c(0.08, 0.90),
        legend.text = element_text(size = 20)
    )

theme_set(my_theme)

exons <- combined_annotation %>%
  filter(type == "exon")

CDS <- combined_annotation %>%
  filter(type == "CDS")

Fig_6E <- exons %>% 
  ggplot(
    aes(
      xstart = start,
      xend = end,
      y = transcript_id,
      color = transcript_id,
      fill = transcript_id
    )
  ) +
  geom_range(
    height = 0.25,
    linewidth = 1
  ) +
  geom_intron(
      data = to_intron(exons, "transcript_id"),
      aes(strand = strand),
      color = "black"
  ) + 
  geom_range(
      data = CDS,
      linewidth = 1
  ) +
  geom_intron(
      data = to_intron(CDS, "transcript_id"),
      aes(strand = strand),
      arrow.min.intron.length = 500,
      color = "black"
  ) +
  geom_vline(xintercept = variant_coord, color = "red") +
  scale_color_manual(values = setNames(list("#009E73", "#E69F00"), c(gencode_tx, novel_pb_id))) +
  scale_fill_manual(values = setNames(list("#009E73", "#E69F00"), c(gencode_tx, novel_pb_id)))

ggsave("figures/figure_6/KCNQ2_full_transcript.pdf", width = 15, height = 5)

Fig_6E + coord_cartesian(xlim = c(variant_start-1200, variant_start+300))
ggsave("figures/figure_6/KCNQ2_full_transcript_zoomed.pdf", width = 8, height = 5)



# Figure 6 isoformSwitch barplot

isoformswitch <- readRDS("export/IsoformSwitchAnalyzeR/isoformswitch.rds")
isoformFeatures <- isoformswitch$isoformFeatures

isoformFeatures %>% 
  filter(
    isoform_id == novel_pb_id,
    condition_1 == "t04",
    condition_2 == "t04"
  ) %>% 
  pivot_longer(
    cols = c("iso_value_1", "iso_value_2"),
    names_to = "condition",
    values_to = "iso_value"
  ) %>%
  mutate(
    condition = recode(condition, "iso_value_1" = "t00", "iso_value_2" = "t30")
  ) %>% 
  ggplot(aes(x=condition, y=iso_value)) +
  geom_col() +
  labs(
    y = "Isoform Expression"
  ) +
  theme(
    axis.title.x = element_text(size = 25, color = "black")
  )

isoformFeatures %>% 
  filter(
    isoform_id == novel_pb_id
  )
