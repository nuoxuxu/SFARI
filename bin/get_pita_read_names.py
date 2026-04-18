#!/usr/bin/env python3
import argparse
import polars as pl


def main():
    parser = argparse.ArgumentParser(description="Get read names for PITA genes from spearman TSV + BED.")
    parser.add_argument("--spearman", required=True, help="Per-sample spearman TSV")
    parser.add_argument("--bed", required=True, help="7-column BED file with gene_name in last column")
    parser.add_argument("--out", required=True, help="Output file with one read name per line")
    args = parser.parse_args()

    spearman = pl.read_csv(args.spearman, separator="\t")
    pita_genes = spearman.filter(
        (pl.col("exon_type") == "dualT") &
        (pl.col("pvalue") < 0.001) &
        (pl.col("spearman_r") > 0)
    )["gene_name"]

    bed = pl.read_csv(
        args.bed, separator="\t", has_header=False,
        new_columns=["chrom", "chromStart", "chromEnd", "name", "score", "strand", "gene_name"]
    )
    reads = bed.filter(pl.col("gene_name").is_in(pita_genes))["name"].unique()

    with open(args.out, "w") as f:
        f.write("\n".join(reads.to_list()))


if __name__ == "__main__":
    main()
