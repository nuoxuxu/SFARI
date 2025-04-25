import polars as pl
from src.utils import read_gtf

signalp_res = pl.read_csv("export/results_20250325-202643.csv")

pbids = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")\
    .filter(
        pl.col("feature") == "transcript"
    )\
    .unique("transcript_id")["transcript_id"]\
    .to_list()

signalp_res\
    .filter(
        pl.col("Protein_ID").is_in(pbids)
    )\
    .write_csv("export/deeploc2.csv", include_header=True)