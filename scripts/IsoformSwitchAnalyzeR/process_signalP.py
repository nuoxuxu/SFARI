import polars as pl
from src.utils import read_gtf

signalp5 = pl.read_csv("export/orfanage_peptide_summary.signalp5", separator = "\t", comment_prefix="#", new_columns=["ID","Prediction","SP(Sec/SPI)","OTHER","CS Position"])

pbids = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")\
    .filter(
        pl.col("feature") == "transcript"
    )\
    .unique("transcript_id")["transcript_id"]\
    .to_list()

signalp5 = signalp5\
    .filter(
        pl.col("ID").is_in(pbids)
    )

signalp5\
    .write_csv("export/processed.signalp5", separator="\t", include_header=False)