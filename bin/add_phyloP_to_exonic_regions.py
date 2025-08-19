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
            pl.col("seqnames")!="chrM"
        )
    pbw = pyBigWig.open(pbw_path)
    
    chrom_idx = df.columns.index("seqnames")
    coord_idx = df.columns.index("coord")
    
    phyloP_column = df.map_rows(lambda row: tuple(pbw.values(row[chrom_idx], row[coord_idx]-1, row[coord_idx]))).rename({"column_0": "phyloP"})
    
    return pl.concat([df, phyloP_column], how="horizontal")

def reformat(exonic_regions_df):
    """
    Reformat the exonic regions DataFrame to match the expected format.
    
    Args:
        exonic_regions_df (pl.DataFrame): DataFrame containing exonic regions.
        
    Returns:
        pl.DataFrame: Reformatted DataFrame with 'seqnames', 'strand', 'start', 'end', and 'type' columns.
    """
    return exonic_regions_df\
        .unique(["seqnames", "strand", "start", "end"])\
        .filter(pl.col("end")-pl.col("start") > 3)\
        .with_columns(
            coord = pl.int_ranges(pl.col("start") - 1, pl.col("end")).list.sample(2, with_replacement=False)
        )\
        .explode("coord")

def main():
    parser = argparse.ArgumentParser(description="Add phyloP scores to exonic regions.")
    parser.add_argument("--cds_regions_csv", type=str, required=True, help="Path to the CDS regions CSV file.")
    parser.add_argument("--utr_regions_csv", type=str, required=True, help="Path to the UTR regions CSV file.")
    parser.add_argument("--phyloP_bigwig", type=str, required=True, help="Path to the phyloP bigWig file.")

    params = parser.parse_args()

    pl.read_csv(params.utr_regions_csv)\
        .pipe(reformat)\
        .pipe(add_phyloP, pbw_path=params.phyloP_bigwig)\
        .select(["seqnames", "strand", "coord", "type", "phyloP"])\
        .write_csv("UTR_regions_with_phyloP.csv")

    pl.read_csv(params.cds_regions_csv)\
        .pipe(reformat)\
        .pipe(add_phyloP, pbw_path=params.phyloP_bigwig)\
        .select(["seqnames", "strand", "coord", "type", "phyloP"])\
        .write_csv("CDS_regions_with_phyloP.csv")

if __name__ == "__main__":
    main()