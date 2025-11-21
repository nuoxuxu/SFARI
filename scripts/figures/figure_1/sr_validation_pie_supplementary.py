from src.utils import read_gtf, read_SJ, gtf_to_SJ
import polars as pl
from pathlib import Path
import os
from SingleCell.ryp import to_r, to_py, r

def read_refmap(file):
    out = pl.read_csv(file, separator="\t")\
    .filter(pl.col("class_code")=="=")\
    .with_columns(
        pl.col("qry_id_list").str.split(",").map_elements(lambda s: [e.split("|")[1] for e in s], return_dtype=pl.List(pl.String)).alias("qry_id_list"))\
    .explode("qry_id_list")\
    .select(["qry_id_list", "ref_id"])\
    .rename({"qry_id_list": "transcript_id"})
    return out

classification = pl.read_parquet("nextflow_results/V47/final_classification.parquet")

TALON_GENCODE = read_refmap("nextflow_results/V47/compare/TALON_GENCODE_V47.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap") 
TALON_SFARI = read_refmap("nextflow_results/V47/compare/TALON_SFARI.cp_vz_0.75_min_7_recovery_talon_hg38.gtf.refmap")

r(
    """
    library(edgeR)
    library(arrow)
    library(dplyr)

    expression <- read_parquet("nextflow_results/V47/final_expression.parquet")
    dge <- DGEList(counts=as.matrix(expression[, grep("_", colnames(expression))]), genes=expression[, "isoform"])
    log2_cpm <- cpm(dge, log=TRUE)
    mean_log2_cpm <- rowMeans(log2_cpm)
    out <- as_tibble(data.frame(isoform=expression$isoform, mean_log2_cpm=mean_log2_cpm))
    """
)
lr_log2_cpm_edgeR = to_py("out").drop("index")

lr_log2_cpm = pl.read_parquet("nextflow_results/V47/final_expression.parquet")\
    .rename(
        {
            "NPC_1_3": "NPC_1_2",
            "NPC_3_3": "NPC_3_2",
            "CN_1_2": "CN_1_1",
            "CN_1_3": "CN_1_2",
        }
    )\
    .with_columns(
        (pl.selectors.numeric() / pl.sum_horizontal(pl.selectors.numeric()) * 1e6 + 1).log(2)
    )\
    .with_columns(
        mean_log2_cpm = pl.mean_horizontal(pl.selectors.numeric())
    )\
    .select(['isoform', 'mean_log2_cpm'])

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

# Without incorporating GENCODE annotation

LR_SJ = read_gtf('nextflow_results/V47/final_transcripts.gtf')\
    .filter(pl.col('feature') == 'exon')\
    .pipe(gtf_to_SJ)\
    .join(classification["isoform", "structural_category"], left_on = "transcript_id", right_on='isoform', how="left")\
    .with_columns(
        type = pl.when(pl.col('structural_category') == 'full-splice_match').then(pl.lit('known'))
               .when(pl.col('structural_category')!='full-splice_match').then(pl.lit('novel'))
    )\
    .join(lr_log2_cpm_edgeR, left_on="transcript_id", right_on="isoform", how="left")\
    .group_by(["chrom", "start", "end", "strand"])\
    .agg(
        pl.col("type").unique(),
        pl.col("mean_log2_cpm").mean()
    )\
    .with_columns(
        type = pl.when((pl.col("type").list.len()==1) & (pl.col('type').list.first()=="novel")).then(pl.lit('novel'))\
                .when((pl.col("type").list.len()==1) & (pl.col('type').list.first()=="known")).then(pl.lit('known'))\
                .otherwise(pl.lit("shared"))
    )

LR_SJ.write_parquet("nextflow_results/V47/LR_SJ.parquet")

# Another method to define the set of SJ

gencode_V47_SJ = read_gtf(os.path.join(os.getenv('GENOMIC_DATA_DIR', ''), 'GENCODE', 'gencode.v47.annotation.gtf'))\
    .filter(pl.col('feature') == 'exon')\
    .pipe(gtf_to_SJ)\
    .unique(['strand', 'chrom', 'start', 'end'])

LR_SJ = read_gtf('nextflow_results/V47/final_transcripts.gtf')\
    .filter(pl.col('feature') == 'exon')\
    .pipe(gtf_to_SJ)\
    .unique(['strand', 'chrom', 'start', 'end'])

full_SJ_label = pl.concat([gencode_V47_SJ, LR_SJ], how='vertical')\
    .join(
        classification["isoform", "structural_category"], left_on = "transcript_id", right_on='isoform', how="left"
    )\
    .with_columns(
        type = pl.when(pl.col("structural_category").is_null()).then(pl.lit("known"))\
        .when(pl.col("structural_category") == "full-splice_match").then(pl.lit("known"))\
        .when(pl.col("structural_category") != "full-splice_match").then(pl.lit("novel"))
    )

out = read_gtf('nextflow_results/V47/final_transcripts.gtf')\
    .pipe(gtf_to_SJ)\
    .join(full_SJ_label["chrom", "strand", "start", "end", "type"], on = ["chrom", "strand", "start", "end"], how="left")\
    .join(lr_log2_cpm_edgeR, left_on="transcript_id", right_on="isoform", how="left")\
    .join(SR_SJ["chrom", "start", "end", "strand", "SR"], on=["chrom", "start", "end", "strand"], how="left")\
    .group_by(["chrom", "start", "end", "strand", "SR"])\
    .agg(
        pl.col("type").unique(),
        pl.col("mean_log2_cpm").mean()
    )\
    .with_columns(
        type = pl.when((pl.col("type").list.len()==1) & (pl.col('type').list.first()=="novel")).then(pl.lit('novel'))\
                .when((pl.col("type").list.len()==1) & (pl.col('type').list.first()=="known")).then(pl.lit('known'))\
                .otherwise(pl.lit("shared"))
    )

out.write_parquet("nextflow_results/V47/LR_SJ_2.parquet")

to_r(out, "LR_SJ")
r(
    """
    library(ggplot2)
    LR_SJ %>% 
        ggplot(aes(x=type, y=mean_log2_cpm, fill=SR)) +
        geom_boxplot() +
        labs(x = "Splice junctions", y = "Long-read RNA-seq expression\n(log2(CPM + 1))")
    ggsave("figures/figure_1/sj_expression_by_whether_validated_3.png", width=4, height=4)
    """
)

# Use Patowaey dataset as a validation

classification = pl.read_parquet("nextflow_results/V47/final_classification.parquet")

lr_log2_cpm = pl.read_parquet("nextflow_results/V47/final_expression.parquet")\
    .rename(
        {
            "NPC_1_3": "NPC_1_2",
            "NPC_3_3": "NPC_3_2",
            "CN_1_2": "CN_1_1",
            "CN_1_3": "CN_1_2",
        }
    )\
    .with_columns(
        (pl.selectors.numeric() / pl.sum_horizontal(pl.selectors.numeric()) * 1e6 + 1).log(2)
    )\
    .with_columns(
        mean_log2_cpm = pl.mean_horizontal(pl.selectors.numeric())
    )\
    .select(['isoform', 'mean_log2_cpm'])

classification\
    .select(["isoform", "structural_category"])\
    .with_columns(
        type = pl.when(pl.col("structural_category") == "full-splice_match")\
            .then(pl.lit('known'))\
            .when(pl.col("structural_category") != "full-splice_match")\
            .then(pl.lit('novel'))
    )\
    .join(lr_log2_cpm, left_on="isoform", right_on="isoform", how="left")\
    .drop("structural_category")\
    .join(TALON_SFARI.unique("ref_id"), left_on="isoform", right_on="ref_id", how="left")\
    .with_columns(
        supported = pl.when(pl.col("transcript_id").is_null()).then(pl.lit(False)).otherwise(pl.lit(True))
    )\
    .write_parquet("nextflow_results/V47/LR_patowary.parquet")

gtf = read_gtf(os.path.join(os.getenv('GENOMIC_DATA_DIR', ''), 'GENCODE', 'gencode.v47.annotation.gtf')).filter(pl.col("feature")=="exon")