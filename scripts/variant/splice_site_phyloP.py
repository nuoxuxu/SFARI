import pyBigWig
import polars as pl
from src.utils import gtf_to_SJ, read_gtf
import os

bw = "data/hg38.phyloP100way.bw"
pbw = pyBigWig.open(bw)

def phylop(row):
    """ Get the phyloP score for a given variant.
    """
    return pbw.values(row["chrom"], row["pos"]-1, row["pos"])[0]

def add_phylop_to_df(df):
    Acceptor = df\
        .map_rows(lambda row: tuple(pbw.values(row[1], row[3]-1-25, row[3]+10)))\
        .drop_nans()\
        .mean()\
        .transpose(header_name="mean")\
        .hstack([pl.Series("pos", range(-25, 11))])\
        .with_columns(
            region = pl.lit("Acceptor")
        )

    Donor = df\
        .map_rows(lambda row: tuple(pbw.values(row[1], row[2]-1-10, row[2]+10)))\
        .drop_nans()\
        .mean()\
        .transpose(header_name="mean")\
        .hstack([pl.Series("pos", range(-10, 11))])\
        .with_columns(
            region = pl.lit("Donor")
        )

    export = pl.concat([Acceptor, Donor], how="vertical")
    return export

def export_phyloP(feature, out):
    """ Export the phyloP scores to a CSV file.
    """
    final_transcripts_SJ = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")\
        .filter(pl.col("feature") == feature)\
        .pipe(gtf_to_SJ)\
        .select("transcript_id", "chrom", "start", "end")

    gencode_gtf = os.getenv("GENOMIC_DATA_DIR") + "/GENCODE/gencode.v47.annotation.gtf" # type: ignore
    GENCODE_SJ = read_gtf(gencode_gtf)\
        .filter(pl.col("feature") == feature)\
        .pipe(gtf_to_SJ)\
        .select("transcript_id", "chrom", "start", "end")
    
    known_SJ_ss = GENCODE_SJ\
        .unique(["chrom", "start", "end"])\
        .pipe(add_phylop_to_df)\
        .with_columns(
            spl_type = pl.lit("known")
        )

    novel_3prime_ss = final_transcripts_SJ\
        .join(
                GENCODE_SJ,
                on=["chrom", "start"]
            )\
        .filter(
            pl.col("end").is_in(GENCODE_SJ["end"]).not_()
        ).unique(["chrom", "start", "end"])\
        .pipe(add_phylop_to_df)\
        .with_columns(
            spl_type = pl.lit("novel_3prime")
        )

    novel_5prime_ss = final_transcripts_SJ\
        .join(
                GENCODE_SJ,
                on=["chrom", "end"]
            )\
        .filter(
            pl.col("start").is_in(GENCODE_SJ["start"]).not_()
        ).unique(["chrom", "start", "end"])\
        .pipe(add_phylop_to_df)\
        .with_columns(
            spl_type = pl.lit("novel_5prime")
        )

    novel_both_ss = final_transcripts_SJ\
        .filter(
            pl.col("start").is_in(GENCODE_SJ["start"]).not_() &
            pl.col("end").is_in(GENCODE_SJ["end"]).not_()
        ).unique(["chrom", "start", "end"])\
        .pipe(add_phylop_to_df)\
        .with_columns(
            spl_type = pl.lit("novel_both")
        )

    export = pl.concat([known_SJ_ss, novel_3prime_ss, novel_5prime_ss, novel_both_ss], how="vertical")
    export.write_csv(out)

export_phyloP("exon", "export/variant/exon_ss_phyloP.csv")
export_phyloP("CDS", "export/variant/CDS_ss_phyloP.csv")