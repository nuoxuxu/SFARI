from src.utils import read_gtf
import polars as pl

annot_peptides_hybrid = read_gtf("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf", attributes=["gene_name", "transcript_id", "detected", "type", "novelty"])
gencode = read_gtf("/Users/xunuo/Genomic_references/GENCODE/gencode.v47.annotation.gtf", attributes=["gene_name", "gene_id", "transcript_id"])
final_classification = pl.read_parquet("nextflow_results/V47/final_classification.parquet")
peptide_mapping = pl.read_parquet("nextflow_results/V47/orfanage/peptide_mapping.parquet")
genes_to_plot = peptide_mapping\
    .join(
        final_classification.rename({"isoform": "pb"})["pb", "structural_category", "associated_gene", "chrom", "strand"], 
        on = ["pb"],
        how = "left"
    )\
    .join(
        annot_peptides_hybrid.rename({"transcript_id": "peptide"}),
        on = ["peptide"],
        how = "left"
    )\
    .filter(
        pl.col("structural_category") == "novel_not_in_catalog",
        pl.col("detected") == "True",
        pl.col("novelty") == "novel"
    )\
    .group_by("pb")\
    .agg(
        pl.col("associated_gene").first(),
        pl.col("peptide").unique().map_elements(lambda x: len(x)).alias("n_peptides")
    )\
    .filter(
        pl.col("n_peptides") == 2
    )
gencode_genes_to_plot = gencode\
    .filter(
        pl.col("feature") == "transcript",
        pl.col("gene_name").is_in(genes_to_plot["associated_gene"])
    )\
    .select("seqname", "gene_name", "transcript_id", "strand")\
    .with_columns(
        pl.lit("gencode")
    )
isoseq_genes_to_plot = final_classification\
    .filter(
        pl.col("associated_gene").is_in(genes_to_plot["associated_gene"])
    )\
    .select("chrom", "associated_gene", "isoform", "strand")\
    .with_columns(
        pl.lit("isoseq")
    )\
    .rename(
        {"chrom": "seqname", "associated_gene": "gene_name", "isoform": "transcript_id"}
    )
plotgene_transcripts_codingOnly = pl.concat([gencode_genes_to_plot, isoseq_genes_to_plot], how = "vertical")
plotgene_transcripts_codingOnly\
    .join(
        final_classification.rename({"isoform": "transcript_id"})["transcript_id", "associated_transcript"],
        on = ["transcript_id"],
        how = "left"
    )\
    .with_columns(
        pl.when(pl.col("associated_transcript").is_null().or_(pl.col("associated_transcript")=="novel"))\
            .then(pl.col("transcript_id"))\
            .otherwise(pl.col("associated_transcript"))
    ).drop("associated_transcript")\
    .write_csv("export/plotgene_transcripts_codingOnly.txt", separator="\t", include_header=False)