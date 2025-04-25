import polars as pl
from src.utils import read_gtf

pfam = pl.read_csv("export/full_pfam_scan_results.txt")

pbids = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")\
    .filter(
        pl.col("feature") == "transcript"
    )\
    .unique("transcript_id")["transcript_id"]\
    .to_list()

pfam\
    .filter(
        pl.col("seq_id").is_in(pbids)
    )\
    .write_csv("export/processed_pfam.txt", include_header=True)