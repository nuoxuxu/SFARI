from src.utils import read_gtf, read_SJ, gtf_to_SJ
import polars as pl
from pathlib import Path
import os

SR_SJ = (
    pl.concat([read_SJ(file) for file in Path('data/riboseq').rglob('*_SJ.out.tab')], how='vertical')
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
    read_gtf('nextflow_results/V47/orfanage/orfanage.gtf')
    .filter(pl.col('feature') == 'CDS')
    .pipe(gtf_to_SJ)
)

gencode_V47_SJ = (
    read_gtf(''.join([os.getenv('GENOMIC_DATA_DIR'), '/GENCODE/gencode.v47.annotation.gtf']))
    .filter(pl.col('feature') == 'CDS')
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

LR_SJ_novel.write_parquet("data/riboseq/riboseq_SJ.parquet")