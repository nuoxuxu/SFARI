library(reticulate)
library(dplyr)
library(ggplot2)
library(patchwork)

py_run_string("
from src.utils import read_gtf, read_SJ, gtf_to_SJ
import polars as pl
from pathlib import Path
import os

SR_SJ = (
    pl.concat([read_SJ(file) for file in Path('STAR_results').rglob('*_SJ.out.tab')], how='vertical')
    .unique(['chrom', 'start', 'end', 'strand'])
    .select(['chrom', 'start', 'end', 'strand'])
    .with_columns(
        strand = pl.col('strand').map_elements(lambda s: '+' if s == 1 else '-', return_dtype=pl.String)
    )
    .with_columns(
        SR = pl.lit(True)
    )
)

LR_SJ = (
    read_gtf('nextflow_results/V47/final_transcripts.gtf')
    .pipe(gtf_to_SJ)
)

gencode_V47_SJ = (
    read_gtf(''.join([os.getenv('GENOMIC_DATA_DIR'), '/GENCODE/gencode.v47.annotation.gtf']))
    .filter(pl.col('feature') == 'exon')
    .pipe(gtf_to_SJ)
)

LR_SJ_novel = (
    LR_SJ
    .filter(pl.col('start').is_null().not_())
    .with_columns(
        LR=pl.lit(True)
    )
    .join(
        gencode_V47_SJ['chrom', 'start', 'end', 'strand'].with_columns(GENCODE=pl.lit(True)),
        on=['chrom', 'start', 'end', 'strand'],
        how='full',
        coalesce=True
    )
    .join(
        SR_SJ,
        on=['chrom', 'start', 'end', 'strand'],
        how='full',
        coalesce=True
    )
)

LR_SJ_novel = (
    LR_SJ_novel
    .with_columns(
        pl.col('LR').fill_null(False),
        pl.col('GENCODE').fill_null(False),
        pl.col('SR').fill_null(False)
    )
)
")

LR_SJ_novel <- as_tibble(py$LR_SJ_novel$to_pandas())

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
ggsave("test.pdf", width = 200, height = 100, units ="mm")