library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(scatterpie)

original <- read_excel("data/mmc2.xlsx", sheet = "Table S2C", skip = 1) %>%
    separate(
        col = Variant,
        into = c("chr", "pos", "ref", "alt"),
        sep = ":",
        convert = TRUE
    ) %>%
    mutate(
        chr = gsub("chr", "", chr)
    ) %>%
    # insert a column identical to pos next to pos
    mutate(
        pos2 = pos
    ) %>%
    filter(nchar(ref) == 1)

orfanage_predicted <- read_tsv("export/variant/annovar/de_novo/orfanage/.hg38_multianno.txt")
gencode_predicted <- read_tsv("export/variant/annovar/de_novo/gencode/.hg38_multianno.txt")

orfanage_func_to_be_replaced <- orfanage_predicted %>%
    distinct(Func.refGene, ExonicFunc.refGene) %>%
    filter(ExonicFunc.refGene == ".") %>%
    pull(Func.refGene)

orfanage_predicted <- orfanage_predicted %>%
    mutate(
        orfanage_variant_type = case_when(
            Func.refGene %in% orfanage_func_to_be_replaced ~ Func.refGene,
            .default = ExonicFunc.refGene
        )
    )

gencode_func_to_be_replaced <- gencode_predicted %>%
    distinct(Func.refGene, ExonicFunc.refGene) %>%
    filter(ExonicFunc.refGene == ".") %>%
    pull(Func.refGene)

gencode_predicted <- gencode_predicted %>%
    mutate(
        gencode_variant_type = case_when(
            Func.refGene %in% gencode_func_to_be_replaced ~ Func.refGene,
            .default = ExonicFunc.refGene
        )
    )

gencode_vt_ranked_by_severity <- c("splicing", "frameshift substitution", "stopgain", "stoploss", "startloss", "nonframeshift substitution", "nonsynonymous SNV", "synonymous SNV", "unknown", "ncRNA_splicing", "nc_RNA_exonic", "nc_RNA_exonic;splicing", "UTR5;UTR3", "UTR5", "UTR3", "intronic", "ncRNA_intronic", "upstream;downstream", "upstream", "downstream", "intergenic")
orfanage_vt_ranked_by_severity <- c("splicing", "frameshift substitution", "stopgain", "stoploss", "startloss", "nonframeshift substitution", "nonsynonymous SNV", "synonymous SNV", "unknown", "UTR5;UTR3", "UTR5", "UTR3", "upstream;downstream", "intronic", "upstream", "downstream", "intergenic")

gencode_predicted %>% distinct(gencode_variant_type) %>% View()
orfanage_predicted %>% distinct(orfanage_variant_type) %>% View()

df <- bind_cols(original[, "ASD status"], gencode_predicted[, "gencode_variant_type"], orfanage_predicted[, "orfanage_variant_type"]) %>%
    group_by(gencode_variant_type, orfanage_variant_type, `ASD status`) %>%
    summarise(
        n = n(),
        .groups = "drop"
    ) %>%
    # long to wide
    pivot_wider(
        names_from = `ASD status`,
        values_from = n,
        values_fill = 0
    ) %>%
    rename(
        "control" = `1`,
        "case" = `2`
    ) %>%
    mutate(
        gencode_variant_type = as.integer(factor(gencode_variant_type, levels = gencode_vt_ranked_by_severity)),
        orfanage_variant_type = as.integer(factor(orfanage_variant_type, levels = orfanage_vt_ranked_by_severity))
    )

df$region <- factor(seq_len(nrow(df)))
df$total <- log(df$case + df$control + 1, 2)^0.5 / 6.5

df %>%
    ggplot() +
    geom_scatterpie(aes(x = gencode_variant_type, y = orfanage_variant_type, r = total, group = region), cols = c("case", "control")) +
    scale_x_continuous(breaks = seq(1, 21), labels=gencode_vt_ranked_by_severity) +
    scale_y_continuous(breaks = seq(1, 17), labels=orfanage_vt_ranked_by_severity) +
    coord_fixed() +
    theme_bw(base_size = 16) +
    theme(
        legend.position = "bottom",
        text = element_text(family = "Arial"),
        panel.grid.major = element_line(size = 0.5),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    ) +
    ylab("Reassigned annovar Consequence") +
    xlab("Previous annovar Consequence")
