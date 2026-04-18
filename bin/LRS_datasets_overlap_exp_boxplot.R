library(ggplot2)
library(dplyr)
library(arrow)
library(edgeR)
library(ggpubr)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
expression_path  <- args[1]
lr_patowary_path <- args[2]
encode4_path     <- args[3]

my_theme <- theme_classic() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.x  = element_text(size = 14),
        axis.text.y  = element_text(size = 14),
        strip.text   = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)
    )
theme_set(my_theme)

colorVector <- c("TRUE" = "#00BFC4", "FALSE" = "#F8766D")

# ---- mean log2 CPM --------------------------------------------------------
expression <- read_parquet(expression_path)
dge <- DGEList(
    counts = as.matrix(expression[, grep("_", colnames(expression))]),
    genes  = expression[, "isoform"]
)
mean_log2_cpm <- data.frame(
    isoform       = expression$isoform,
    mean_log2_cpm = rowMeans(cpm(dge, log = TRUE))
)

# ---- Patowary et al. overlap -----------------------------------------------
patowary_df <- read_parquet(lr_patowary_path) %>%
    select(isoform, type, supported) %>%
    rename(in_comparison = supported) %>%
    mutate(comparison = "Patowary et al.") %>%
    left_join(mean_log2_cpm, by = "isoform")

# ---- ENCODE4 overlap -------------------------------------------------------
encode4_df <- read_parquet(encode4_path) %>%
    rename(in_comparison = in_encode4) %>%
    mutate(comparison = "ENCODE4") %>%
    left_join(mean_log2_cpm, by = "isoform")

# ---- Combined long-format dataframe ----------------------------------------
combined <- bind_rows(patowary_df, encode4_df) %>%
    mutate(
        type          = recode(type, "known" = "GENCODE", "novel" = "Novel"),
        in_comparison = factor(in_comparison, levels = c(TRUE, FALSE))
    )

# ---- Known transcripts figure ----------------------------------------------
known_overlap_plot <- combined %>%
    filter(type == "GENCODE") %>%
    ggplot(aes(x = comparison, y = mean_log2_cpm, fill = in_comparison)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual("Supported", values = colorVector) +
    stat_compare_means(
        aes(label = after_stat(p.signif)),
        method = "t.test"
    ) +
    labs(
        y     = expression("Mean" ~ log[2] * "(CPM + 1) across all samples"),
        title = "Known (GENCODE) transcripts"
    )

# ---- Novel transcripts figure ----------------------------------------------
novel_overlap_plot <- combined %>%
    filter(type == "Novel") %>%
    ggplot(aes(x = comparison, y = mean_log2_cpm, fill = in_comparison)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual("Supported", values = colorVector) +
    stat_compare_means(
        aes(label = after_stat(p.signif)),
        method = "t.test"
    ) +
    labs(
        y     = expression("Mean" ~ log[2] * "(CPM + 1) across all samples"),
        title = "Novel transcripts"
    )

# ---- Combined figure -------------------------------------------------------
combined_plot <- (known_overlap_plot + novel_overlap_plot) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 20, face = "bold"))

ggsave(
    "lrs_overlap_expression.pdf",
    combined_plot,
    width = 12, height = 5
)
