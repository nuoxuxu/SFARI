from src.utils import read_gtf, read_SJ, gtf_to_SJ
import polars as pl
from pathlib import Path
import os

SR_SJ = (
    pl.concat([read_SJ(file) for file in Path('export/STAR_results').rglob('*_SJ.out.tab')], how='vertical')
    .unique(['chrom', 'start', 'end', 'strand', 'motif'])
    .select(['chrom', 'start', 'end', 'strand', 'motif'])
    .with_columns(
        strand = pl.col('strand').map_elements(lambda s: '+' if s == 1 else '-', return_dtype=pl.String)
    )
    .with_columns(
        SR = pl.lit(True)
    )
)

LR_SJ = read_gtf('nextflow_results/V47/final_transcripts.gtf')\
    .filter(pl.col('feature') == 'CDS')\
    .pipe(gtf_to_SJ)\
    .group_by(['chrom', 'start', 'end', 'strand'])\
    .len()

gtf = read_gtf('nextflow_results/V47/final_transcripts.gtf')\
    .filter(pl.col('feature') == 'exon')

gtf\
    .filter(
        pl.col("feature")!="transcript"
    )\
    .filter(
        pl.col("transcript_id").is_in(gtf.group_by("transcript_id").len().filter(pl.col("len")!=1).select("transcript_id")["transcript_id"])
    )\
    .group_by("transcript_id", maintain_order=True)\
    .agg(
        pl.col("strand").unique().map_elements(lambda l: l[0], return_dtype=pl.String),
        pl.col("seqname").unique().map_elements(lambda l: l[0], return_dtype=pl.String),
        pl.col("start").map_elements(lambda l: np.sort(np.array(l)-1)[1:].tolist(), return_dtype=pl.List(pl.Int64)).alias("end"),
        pl.col("end").map_elements(lambda l: np.sort(np.array(l)+1)[:-1].tolist(), return_dtype=pl.List(pl.Int64)).alias("start")
    )\
    .explode("start", "end")\
    .rename({"seqname": "chrom"})

gencode_V47_SJ = (
    read_gtf(os.path.join(os.getenv('GENOMIC_DATA_DIR', ''), 'GENCODE', 'gencode.v47.annotation.gtf'))
    .filter(pl.col('feature') == 'exon')
    .pipe(gtf_to_SJ)\
    .unique(['strand', 'chrom', 'start', 'end'])
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

LR_SJ_novel.write_parquet('export/LR_SJ_novel.parquet')