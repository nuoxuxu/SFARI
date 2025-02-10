#!/usr/bin/env python3
from src.utils import read_fasta
import polars as pl
import argparse

def main():
    parser = argparse.ArgumentParser("Get best_orf.tsv")
    parser.add_argument("--final_sample_classification", action='store', type=str, required=True)
    parser.add_argument("--final_sample_fasta", action='store', type=str, required=True)
    parser.add_argument("--cds_fasta", action='store', type=str, required=True)
    parser.add_argument("--output", action='store', type=str, required=True)

    params = parser.parse_args()
    
    def get_CDS_start_end_coord(transcript_seq, CDS_seq):
        
        # Find start position
        start_pos = transcript_seq.find(CDS_seq)

        # Compute end position
        if start_pos != -1:  # Ensure substring exists
            end_pos = start_pos + len(CDS_seq) - 1
            return start_pos, end_pos
        else:
            print("Substring not found.")

    cds = read_fasta(params.cds_fasta)
    full_tx = read_fasta(params.final_sample_fasta)

    df = cds\
        .rename({"seq": "cds_seq"})\
        .join(
            full_tx.rename({"seq": "full_seq"}),
            on="transcript_id",
            how="left"
        )

    df = pl.concat([df.select("transcript_id"), df.map_rows(lambda row: get_CDS_start_end_coord(row[2], row[1]))], how="horizontal")\
        .rename(
            {"column_0": "start", "column_1": "end"}
        )\
        .with_columns(
            pl.col("start") + 1,
            pl.col("end") + 1
        )

    best_orf = pl.read_parquet(params.final_sample_classification)\
        ["isoform", "length"]\
        .with_columns(
            orf_frame = 1
        )\
        .rename({"isoform": "pb_acc", "length": "len"})

    df = df\
        .rename({"transcript_id": "pb_acc", "start": "orf_start", "end": "orf_end"})

    best_orf = best_orf\
        .join(
            df, on = "pb_acc", how = "left"
        )\
        .drop_nulls("orf_start")\
        .with_columns(
            orf_len = pl.col("orf_end") - pl.col("orf_start") + 1
        )

    best_orf.write_csv(params.output, separator="\t")

if __name__ == "__main__":
    main()