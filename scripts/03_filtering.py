import polars as pl
from src.single_cell import SingleCell
import numpy as np

lr_bulk = SingleCell("results/long_read/pbid.h5ad")
classification = pl.read_csv("proc/merged_collapsed_classification.filtered_lite_classification.txt", separator="\t", null_values="NA")

# Filter pacbio id count matrix to keep only transcripts in `classification`
# Filter Settings:
#   -a,--polya-percent     FLOAT  Adenine percentage at genomic 3' end to flag an isoform as intra-priming. [0.6]
#   -r,--polya-run-length  INT    Continuous run-A length at genomic 3' end to flag an isoform as intra-priming. [6]
#   -m,--max-distance      INT    Maximum distance to an annotated 3' end to preserve as a valid 3' end. [50]
#   -c,--min-cov           INT    Minimum junction coverage for each isoform. Only used if min_cov field is not 'NA'. [3]
#   --mono-exon                   Filter out all mono-exonic transcripts.
#   --skip-junctions              Skip junctions.txt filtering.

lr_bulk = lr_bulk.filter_var( pl.col("pbid").is_in(classification["isoform"]))
lr_bulk = lr_bulk[:, np.where(np.sum(lr_bulk.X > 5, axis = 0) > 2)[0]]

lr_bulk.save("results/long_read/pbid_filtered.h5ad")