from pathlib import Path
import re
import polars as pl
from src.utils import fastq_dir_to_samplesheet

fastq_dir_to_samplesheet("data/SFARI_data", "assets/samplesheet.csv", strandedness="unstranded")
pattern = re.compile(r'([^_/]+_[^_/]+_[^_/]+)_')

pl.read_csv("assets/samplesheet.csv")\
    .with_columns(
        pl.col("sample").map_elements(lambda s: pattern.search(s).group(1), return_dtype = pl.String)
    )\
    .with_columns(
        pl.when(pl.col("sample").str.starts_with("iPSC"))\
            .then(pl.col("sample").map_elements(lambda s: s.rsplit("_", 1)[0], return_dtype = pl.String))\
            .otherwise(pl.col("sample"))
            
    )\
    .write_csv("assets/samplesheet.csv", include_header = True)