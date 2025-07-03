library(reticulate)
library(dplyr)
library(ggplot2)
library(patchwork)
library(arrow)

LR_SJ_novel <- read_parquet("export/LR_SJ_novel.parquet")

classification <- read_parquet("nextflow_results/V47/final_classification.parquet")

LR_SJ_novel %>%
    filter(LR, GENCODE) %>%
    filter(SR) %>%
    left_join(
        classification %>% select(c(isoform, structural_category)),
        by = join_by(transcript_id == isoform)
    ) %>%
    group_by(structural_category) %>%
    summarise(len = n()) %>%
    mutate(
        percent = len / sum(len) * 100
    ) %>%
    ggplot(
        aes(x = structural_category, y = len, fill = structural_category)
    ) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", percent)), vjust = -0.5, color = "black", size = 3.5) +
    labs(title = "Validated by short-read splice junctions") +
    theme_minimal()

LR_SJ_novel %>%
    filter(LR, !GENCODE) %>%
    filter(!SR) %>%
    left_join(
        classification %>% select(c(isoform, structural_category)),
        by = join_by(transcript_id == isoform)
    ) %>%
    group_by(structural_category) %>%
    summarise(len = n()) %>%
    mutate(
        percent = len / sum(len) * 100
    ) %>%
    ggplot(
        aes(x = structural_category, y = len, fill = structural_category)
    ) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f%%", percent)), vjust = -0.5, color = "black", size = 3.5) +
    labs(title = "Not validated by short-read splice junctions") +
    theme_minimal()

LR_SJ_novel %>%
    filter(LR, GENCODE) %>%
    filter(!SR) %>% 
    ggplot(aes(x=len)) +
    geom_histogram(bins = 100) +
    xlim(0, 20)

LR_SJ_novel %>%
    filter(LR, !GENCODE) %>%
    filter(!SR) %>%
    ggplot(aes(x=len)) +
    geom_histogram(bins = 100) +
    xlim(0, 20) +
    labs(title = "SJ not validated by short-read splice junctions")

LR_SJ_novel %>%
    filter(LR, !GENCODE) %>%
    filter(SR) %>%
    ggplot(aes(x=len)) +
    geom_histogram(bins = 100) +
    xlim(0, 20) +
    labs(title = "SJ validated by short-read splice junctions")

df <- LR_SJ_novel %>%
    filter(LR, GENCODE) %>%
    count(SR, name = "len") %>%
    mutate(
        percent = len / sum(len) * 100
    ) %>%
    mutate(
        SR = if_else(SR, "Validated", "Not validated")
    ) %>%
    arrange(desc(SR)) %>%
    mutate(ypos = cumsum(percent) - 0.5 * percent)

p1 <- df %>%
    ggplot(aes(x = "", y = percent, fill = SR)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    labs(fill = "Validation by\n short-read\n splice junctions") +
    theme_void() +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 5)

df <- LR_SJ_novel %>%
    filter(LR, !GENCODE) %>%
    count(SR, name = "len") %>%
    mutate(
        percent = len / sum(len) * 100
    ) %>%
    mutate(
        SR = if_else(SR, "Validated", "Not validated")
    ) %>%
    arrange(desc(SR)) %>%
    mutate(ypos = cumsum(percent) - 0.5 * percent)

p2 <- df %>%
    ggplot(aes(x = "", y = percent, fill = SR)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    labs(fill = "Validation by\n short-read\n splice junctions") +
    theme_void() +
    geom_text(aes(y = ypos, label = sprintf("%.1f%%", percent)), color = "white", size = 5)

p1 + p2
ggsave("figures/figure_1/sr_validation_pie.pdf", width = 200, height = 100, units = "mm")
