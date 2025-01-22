#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
import polars as pl
from src.utils import read_gtf
import argparse

def main():
    parser = argparse.ArgumentParser(description='Add base protein classification and gene name to TransDecoder-predicted transcripts')
    parser.add_argument("--genome_gff3_gtf", action="store", type=str, required=True)
    parser.add_argument("--unfiltered_protein_classification", action="store", type=str, required=True)
    parser.add_argument("--reference_gtf", action="store", type=str, required=True)
    parser.add_argument("--output", action="store", type=str, required=True)
    
    params = parser.parse_args()
    CDS_gtf  = read_gtf(
        params.genome_gff3_gtf,
        attributes = ["transcript_id", "gene_id"]
        )\
        .filter(
            pl.col("feature")!="transcript",
            pl.col("feature") == "CDS"
            )\
        ["seqname", "feature", "start", "end", "strand", "gene_id", "transcript_id"]


    protein_classification = pl.read_csv(params.unfiltered_protein_classification, separator = "\t").rename({"pb": "transcript_id"})["transcript_id", "pr_gene", "protein_classification_base"]

    gene_id_to_gene_name = read_gtf(params.reference_gtf, attributes = ["gene_name", "gene_id"]).unique("gene_id")["gene_name", "gene_id"]

    CDS_gtf = CDS_gtf\
        .join(
            protein_classification,
            on = "transcript_id",
            how = "left"
        )\
        .join(
            gene_id_to_gene_name.rename({"gene_id": "pr_gene"}),
            on = "pr_gene",
            how = "left"
        )["seqname", "feature", "start", "end", "strand", "gene_name", "transcript_id", "protein_classification_base"]

    CDS_gtf.write_csv(params.output)