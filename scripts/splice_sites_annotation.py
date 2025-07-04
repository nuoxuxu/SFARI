from src.utils import read_gtf, read_SJ
import polars as pl
from Bio import SeqIO
import os
import pyBigWig
from pathlib import Path

def gtf_to_SJ(gtf):
    import numpy as np

    assert gtf["feature"].unique().is_in(["exon", "CDS"]).all(), "GTF must only contain exon and CDS features"

    return gtf\
        .filter(
            pl.col("feature")!="transcript"
        )\
        .group_by(["transcript_id", "feature"], maintain_order=True)\
        .agg(
            pl.col("seqname").unique().map_elements(lambda l: l[0], return_dtype=pl.String),
            pl.col("start").map_elements(lambda l: np.sort(np.array(l)-1)[1:].tolist(), return_dtype=pl.List(pl.Int64)).alias("end"),
            pl.col("end").map_elements(lambda l: np.sort(np.array(l)+1)[:-1].tolist(), return_dtype=pl.List(pl.Int64)).alias("start"),
            pl.col("strand").unique().map_elements(lambda l: l[0], return_dtype=pl.String)
        )\
        .explode("start", "end")\
        .rename({"seqname": "chrom"})\
        .filter(pl.col("start").is_null().not_())\
        .select("chrom", "feature", "start", "end", "strand", "transcript_id")

def add_phyloP(ss_df, bw, column_prefix):
    pbw = pyBigWig.open(bw)
    
    chrom_idx = ss_df.columns.index("chrom")
    coord_idx = ss_df.columns.index("coord")

    end_phyloP = ss_df\
        .map_rows(lambda row: tuple(pbw.values(row[chrom_idx], row[coord_idx]-2, row[coord_idx])))\
        .rename({"column_0": "end_phyloP_2", "column_1" : "end_phyloP_1"})

    start_phyloP = ss_df\
        .map_rows(lambda row: tuple(pbw.values(row[chrom_idx], row[coord_idx]-1, row[coord_idx]+1)))\
        .rename({"column_0": "start_phyloP_1", "column_1" : "start_phyloP_2"})

    ss_df = pl.concat([ss_df, end_phyloP, start_phyloP], how="horizontal")\
        .with_columns(
            pl.when(pl.col("start_or_end")=="start")\
                .then(pl.col("start_phyloP_1"))\
                .when(pl.col("start_or_end")=="end")\
                .then(pl.col("end_phyloP_1"))\
                .alias(f"{column_prefix}_phyloP_1"),
            pl.when(pl.col("start_or_end")=="start")\
                .then(pl.col("start_phyloP_2"))\
                .when(pl.col("start_or_end")=="end")\
                .then(pl.col("end_phyloP_2"))\
                .alias(f"{column_prefix}_phyloP_2")
        )\
        .drop(
            "end_phyloP_1", "end_phyloP_2", "start_phyloP_1", "start_phyloP_2"
        )
    return ss_df

def unique_and_add_is_CDS(ss_df):
    return ss_df.group_by(
            ["chrom", "strand", "start_or_end", "coord"],
        )\
        .agg(
            pl.col("feature"),
            pl.col("transcript_id").unique().map_elements(lambda l: l[0], return_dtype=pl.String)
        )\
        .with_columns(
            is_CDS = pl.col("feature").list.unique()
        )\
        .with_columns(
            is_CDS = pl.when(pl.col("is_CDS").list.len()==2)\
            .then(pl.lit("CDS"))\
            .otherwise(pl.col("is_CDS").list.first())
        )\
        .with_columns(
            pl.col("is_CDS").replace_strict({"CDS": True, "exon": False})
        ).drop("feature")

def add_canonical(ss_df):
    genome = list(SeqIO.parse(os.getenv("GENOMIC_DATA_DIR") + "/GENCODE/GRCh38.primary_assembly.genome.fa", "fasta"))
    genome = {record.id: record for record in genome}

    chrom_idx = ss_df.columns.index("chrom")
    coord_idx = ss_df.columns.index("coord")

    acceptor_seq = ss_df\
        .map_rows(lambda row: str(genome[row[chrom_idx]].seq[row[coord_idx]-2:row[coord_idx]]))\
        .rename({"map": "acceptor_seq"})

    donor_seq = ss_df\
        .map_rows(lambda row: str(genome[row[chrom_idx]].seq[row[coord_idx]-1:row[coord_idx]+1]))\
        .rename({"map": "donor_seq"})

    return pl.concat([ss_df, acceptor_seq, donor_seq], how="horizontal")\
        .with_columns(
            pl.when(pl.col("start_or_end") == "end")
                .then(pl.col("acceptor_seq"))
                .otherwise(pl.col("donor_seq"))
                .alias("splice_site_seq")
        )\
        .drop("acceptor_seq", "donor_seq")\
        .with_columns(
            canonical = pl.when((pl.col("strand") == "+") & (pl.col("start_or_end") == "start"))\
                .then(pl.col("splice_site_seq") == "GT")\
                .when((pl.col("strand") == "-") & (pl.col("start_or_end") == "end"))\
                .then(pl.col("splice_site_seq") == "AC")\
                .when((pl.col("strand") == "+") & (pl.col("start_or_end") == "end"))\
                .then(pl.col("splice_site_seq") == "AG")\
                .when((pl.col("strand") == "-") & (pl.col("start_or_end") == "start"))\
                .then(pl.col("splice_site_seq") == "CT")\
        )\
        .drop("splice_site_seq")

def ad_SR(ss_df):
    SR_ss = pl.concat([read_SJ(file) for file in Path('export/STAR_results').rglob('*_SJ.out.tab')], how='vertical')\
        .unique(['chrom', 'start', 'end', 'strand'])\
        .select(['chrom', 'start', 'end', 'strand'])\
        .with_columns(
            strand = pl.col('strand').map_elements(lambda s: '+' if s == 1 else '-', return_dtype=pl.String)
        )\
        .unpivot(
            on = ["start", "end"],
            index=["chrom", "strand"],
            variable_name="start_or_end",
            value_name="coord"
        )\
        .with_columns(
            SR = pl.lit(True)
        )\
        .unique(["chrom", "strand", "start_or_end", "coord"])
     
    return ss_df\
        .join(
            SR_ss, 
            on=["chrom", "strand", "start_or_end", "coord"],
            how="left"
        )\
        .with_columns(
            SR = pl.col("SR").fill_null(False)
        )

def add_is_novel(ss_df):
    GENCODE_ss = read_gtf(os.getenv("GENOMIC_DATA_DIR") + "/GENCODE/gencode.v47.annotation.gtf")\
        .filter(pl.col("feature") == "exon")\
        .pipe(gtf_to_SJ)\
        .unpivot(
            on=["start", "end"],
            index=["chrom", "strand", "transcript_id"],
            variable_name="start_or_end",
            value_name="coord"
        )\
        .unique(["chrom", "start_or_end", "coord"])
    
    return ss_df\
        .join(
            GENCODE_ss,
            on=["chrom", "strand", "start_or_end", "coord"],
            how="left"
        )\
        .with_columns(
            pl.when(pl.col("transcript_id_right").is_null())\
                .then(pl.lit(True))\
                .otherwise(pl.lit(False)).alias("is_novel")
        ).drop("transcript_id_right")

def main():
    read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")\
        .filter(pl.col("feature").is_in(["CDS", "exon"]))\
        .filter(pl.col("seqname")!= "chrM")\
        .pipe(gtf_to_SJ)\
        .unpivot(
            on=["start", "end"],
            index=["chrom", "feature", "strand", "transcript_id"],
            variable_name="start_or_end",
            value_name="coord"
        )\
        .pipe(unique_and_add_is_CDS)\
        .pipe(add_canonical)\
        .pipe(ad_SR)\
        .pipe(add_is_novel)\
        .pipe(add_phyloP, "data/cactus241way.phyloP.bw", "mammal")\
        .write_csv("export/ORFanage_splice_sites.csv")

if __name__ == "__main__":
    main()