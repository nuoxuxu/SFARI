import polars as pl

pl.read_csv("nextflow_results/orfanage/novel_peptides.csv")\
    .unique("pep")\
    .select("pep", "transcript_id", "gene_name")\
    .write_csv("export/TableS2_novel_peptides.csv")