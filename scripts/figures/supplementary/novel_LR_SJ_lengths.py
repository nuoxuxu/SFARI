from src.utils import read_gtf, read_SJ, gtf_to_SJ
import polars as pl
from pathlib import Path
import os
from src.ryp import r, to_r

SR_SJ = (
    pl.concat([read_SJ(file) for file in Path('nextflow_results/STAR_results').rglob('*_SJ.out.tab')], how='vertical')
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

novel_LR_SJ = LR_SJ_novel\
    .drop_nulls()\
    .filter(
        pl.col("LR"), 
        pl.col("GENCODE").not_()
    )\
    .with_columns(
        len = pl.col("end") - pl.col("start") +1
    )

to_r(novel_LR_SJ, "novel_LR_SJ")

r(
"""
library(ggplot2)
library(dplyr)

novel_LR_SJ %>%
    ggplot(aes(len, fill=SR)) +
    geom_histogram(alpha=0.5) +
    theme_minimal() + 
    ggtitle("Novel long-read SJ lengths distribution")
ggsave("figures/supplementary/novel_LR_SJ_lengths.pdf")
"""
)
