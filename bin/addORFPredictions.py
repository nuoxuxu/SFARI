#!/home/s/shreejoy/nxu/miniforge3/envs/genomekit/bin/python

import genome_kit as gk
import polars as pl
from src.utils import read_gff
from src.single_cell import SingleCell
import argparse

def main():
    parser = argparse.ArgumentParser(description='Add ORF predictions to the high confidence set')
    parser.add_argument('--h5ad_file', action='store', type=str, required=True)
    parser.add_argument('--longest_orfs_gff3', action='store', type=str, required=True)
    parser.add_argument('--pred_orfs_gff3', action='store', type=str, required=True)
    parser.add_argument('--output', action='store', type=str, required=True)
    params = parser.parse_args()

    genome = gk.Genome("SFARI")

    def get_CDS_coords(pbid):
        return [cds.interval.start for cds in genome.transcripts[pbid].cdss] + [cds.interval.end for cds in genome.transcripts[pbid].cdss]

    def get_start_coord(pbid):
        return genome.transcripts[pbid].cdss[0].start

    def get_end_coord( pbid):
        return genome.transcripts[pbid].cdss[-1].end
    
    lr_bulk = SingleCell(params.h5ad_file)

    longest_orfs_header = read_gff(params.longest_orfs_gff3, attributes=["ID", "Parent", "Name"]).drop_nulls("seqname")\
        .with_columns(
            pl.col("Name").str.strip_chars('"')
        )\
        .filter(pl.col("feature")=="mRNA")

    pred_orfs_header = read_gff(params.pred_orfs_gff3, attributes=["ID", "Parent", "Name"]).drop_nulls("seqname")\
        .with_columns(
            pl.col("Name").str.strip_chars('"')
        )\
        .filter(pl.col("feature")=="mRNA")
    
    ORF_annot = lr_bulk.var["isoform"]\
        .to_frame()\
        .join(
            longest_orfs_header["seqname"].unique().to_frame().with_columns(pl.lit(True).alias("at_least_one_orf")).rename({"seqname": "isoform"}),
            how = "left",
            on = "isoform"
            )\
        .with_columns(
            pl.when(pl.col("at_least_one_orf").is_null()).then(pl.lit(False)).otherwise(pl.lit(True)).alias("at_least_one_orf")
        )

    pred_orfs_annot = pred_orfs_header["seqname", "Name"]\
        .select(
            pl.col("seqname"),
            pl.lit(True).alias("predicted_orf"),
            ORF_type = pl.col("Name").str.extract(r"ORF type:(\w+)")
        )\
        .rename({"seqname": "isoform"})
    
    ORF_annot = ORF_annot\
        .join(
            pred_orfs_annot,
            how = "left",
            on = "isoform"
            )\
        .with_columns(
            pl.when(pl.col("predicted_orf").is_null()).then(pl.lit(False)).otherwise(pl.lit(True)).alias("predicted_orf")
        )
    
    lr_bulk.var = lr_bulk.var\
        .join(
            ORF_annot,
            how = "left",
            on = "isoform"
        )

    CDS_genomic_coords = lr_bulk.var\
        .filter(pl.col("predicted_orf"))\
        .select(
            pl.col("isoform"),
            CDS_genomic_start = pl.col("isoform").map_elements(get_start_coord),
            CDS_genomic_end = pl.col("isoform").map_elements(get_end_coord),
            CDS_coord_list = pl.col("isoform").map_elements(get_CDS_coords)
        )

    lr_bulk.var = lr_bulk.var\
        .join(
            CDS_genomic_coords,
            on = "isoform",
            how = "left"
        )\
        .with_columns(
            pl.when(pl.col("CDS_genomic_start_right")\
                .is_null())\
                .then(pl.lit(None))\
                .otherwise(pl.col("CDS_genomic_start_right"))\
                .alias("CDS_genomic_start"),
            pl.when(pl.col("CDS_genomic_end_right")\
                .is_null())\
                .then(pl.lit(None))\
                .otherwise(pl.col("CDS_genomic_end_right"))\
                .alias("CDS_genomic_end")
        ).drop("CDS_genomic_start_right", "CDS_genomic_end_right")

    coordList_to_baseIsoform = lr_bulk.var\
        .group_by("CDS_coord_list")\
        .agg(
            pl.col("isoform")
        )\
        .with_columns(
            base_isoform = pl.when(pl.col("CDS_coord_list").is_null())\
                .then(pl.lit(None))\
                .otherwise(pl.col("isoform").map_elements(lambda x: x[0], return_dtype=pl.String))
        )\
        .explode("isoform")\
        .drop("CDS_coord_list")

    lr_bulk.var = lr_bulk.var\
        .drop("CDS_coord_list")\
        .join(
            coordList_to_baseIsoform,
            on = "isoform",
            how = "left"
        )

    lr_bulk.save(params.output, overwrite=True)

if __name__ == "__main__":
    main()