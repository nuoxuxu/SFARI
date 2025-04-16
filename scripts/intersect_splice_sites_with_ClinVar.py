from src.utils import read_gtf
import polars as pl
from src.utils import gtf_to_SJ
import os



#--------------------------------Read the CSV with the defined column names--------------------------------
# ClinVar coordinates are 1-based

# Define the column names in a clear list structure

variant_columns = [
    "Variant ID",
    "Location",
    "Variant type",
    "Gene",
    "Molecular consequences",
    "Most severe clinical significance",
    "1000G MAF",
    "GO-ESP MAF",
    "ExAC MAF",
    "Publications (PMIDs)"
]

vv_table = pl.read_csv(
    "data/vv_table_download-2025314155934.txt",
    separator="\t",
    comment_prefix="#",
    new_columns=variant_columns
)

vv_table = vv_table.filter(pl.col("Variant type") == "single nucleotide variant")\
    .with_columns(
        pl.col("Location").cast(pl.Int64)
    )

#--------------------------------Read in final_transcripts and GENCODE SJ--------------------------------
final_transcripts_SJ = read_gtf("nextflow_results/V47/final_transcripts.gtf")\
    .filter(pl.col("feature") == "exon")\
    .pipe(gtf_to_SJ)\
    .unique(["chrom", "start", "end"])
GENCODE_SJ = read_gtf("".join([os.getenv("GENOMIC_DATA_DIR"), "/GENCODE/gencode.v47.annotation.gtf"]))\
    .filter(pl.col("feature") == "exon")\
    .pipe(gtf_to_SJ)\
    .unique(["chrom", "start", "end"])
#--------------------------------Get novel splice sites--------------------------------
novel_splice_sites = final_transcripts_SJ\
    .join(
        GENCODE_SJ["chrom", "start", "end"], 
        on=["chrom", "start", "end"], 
        how="left",
        coalesce=False
        )\
    .filter(pl.col("start_right").is_null())\
    .select(["chrom", "start", "end"])\
    .with_columns(
        end_1 = pl.col("end") + 1,
        start_1 = pl.col("start") - 1
    )\
    .unpivot(
        index="chrom",
        on=["start", "end", "start_1", "end_1"]
    )

novel_splice_sites\
    .join(
        vv_table["Variant ID", "Location", "Gene"],
        left_on="value",
        right_on="Location",
        how="left",
        coalesce=False
    )