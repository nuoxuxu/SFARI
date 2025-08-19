#!/usr/bin/env python3
import polars as pl
import pyBigWig
import argparse

def add_phyloP(df, pbw_path):
    """
    Add phyloP scores to the splice sites DataFrame.
    
    Args:
        df (pl.DataFrame): DataFrame containing splice sites with 'chrom' and 'coord' columns.
        
    Returns:
        pl.DataFrame: DataFrame with an additional 'phyloP' column.
    """
    df = df\
        .filter(
            pl.col("chrom")!="chrM"
        )
    pbw = pyBigWig.open(pbw_path)
    
    chrom_idx = df.columns.index("chrom")
    coord_idx = df.columns.index("coord")
    
    phyloP_column = df.map_rows(lambda row: tuple(pbw.values(row[chrom_idx], row[coord_idx]-1, row[coord_idx]))).rename({"column_0": "phyloP"})
    
    return pl.concat([df, phyloP_column], how="horizontal")

def add_coord(df):
    return df\
    .with_columns([
        pl.when(pl.col("start_or_end") == "start").then(pl.col("start"))
          .when(pl.col("start_or_end") == "end").then(pl.col("end"))
          .alias("coord")
    ])

def main():
    argparser = argparse.ArgumentParser(description="Add phyloP scores to splice sites.")
    argparser.add_argument("--known_splice_sites_cds", type=str, help="Path to known splice sites CDS CSV file.")
    argparser.add_argument("--novel_splice_sites_cds", type=str, help="Path to novel splice sites CDS CSV file.")
    argparser.add_argument("--known_splice_sites_all", type=str, help="Path to known splice sites all CSV file.")
    argparser.add_argument("--novel_splice_sites_all", type=str, help="Path to novel splice sites all CSV file.")
    argparser.add_argument("--phyloP_bigwig", type=str, required=True, help="Path to the phyloP bigWig file.")
    
    
    params = argparser.parse_args()

    known_splice_sites_cds = pl.read_csv(params.known_splice_sites_cds)\
        .with_columns(type=pl.lit("known"))\
        .pipe(add_coord)\
        .pipe(add_phyloP, pbw_path=params.phyloP_bigwig)\
        .drop_nans()

    novel_splice_sites_cds = pl.read_csv(params.novel_splice_sites_cds)\
        .with_columns(type=pl.lit("novel"))\
        .pipe(add_coord)\
        .pipe(add_phyloP, pbw_path=params.phyloP_bigwig)\
        .drop_nans()

    known_splice_sites_cds.write_csv("known_splice_sites_cds_phyloP.csv")
    novel_splice_sites_cds.write_csv("novel_splice_sites_cds_phyloP.csv")

    # Not just CDS, splice junctions from all exons

    known_splice_sites_all = pl.read_csv(params.known_splice_sites_all)\
        .with_columns(type=pl.lit("known"))\
        .pipe(add_coord)\
        .pipe(add_phyloP, pbw_path=params.phyloP_bigwig)\
        .drop_nans()

    novel_splice_sites_all = pl.read_csv(params.novel_splice_sites_all)\
        .with_columns(type=pl.lit("novel"))\
        .pipe(add_coord)\
        .pipe(add_phyloP, pbw_path=params.phyloP_bigwig)\
        .drop_nans()

    known_splice_sites_all.write_csv("known_splice_sites_all_phyloP.csv")
    novel_splice_sites_all.write_csv("novel_splice_sites_all_phyloP.csv")

if __name__ == "__main__":
    main()