#!/bin/env python
import polars as pl
from tqdm import tqdm
from src.utils import read_gtf
import argparse

def make_IL_regex(seq: str) -> str:
    pattern = ''.join('[IL]' if c in 'IL' else c for c in seq)
    return f'{pattern}'

def read_fasta(fasta_file, gencode = False):
    """
    Reads a FASTA file and converts it into a Polars DataFrame.
    Removes '*' from sequences.
    
    Args:
        fasta_file (str): Path to the FASTA file.
    
    Returns:
        polars.DataFrame: A DataFrame with 'transcript_id' and 'seq' columns.
    """
    sequences = []
    transcript_id = None
    seq = []
    
    with open(fasta_file, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Header line
                if transcript_id is not None:  # Save previous entry
                    sequences.append((transcript_id, "".join(seq).replace("*", "")))  # Strip '*'
                # Extract transcript_id (substring before the first space)
                if gencode is False:
                    transcript_id = line[1:].split(" ", 1)[0]
                else:
                    transcript_id = line[1:].split("|")[1]
                seq = []  # Reset sequence
            else:
                seq.append(line)  # Collect sequence lines

    # Add the last sequence
    if transcript_id is not None:
        sequences.append((transcript_id, "".join(seq).replace("*", "")))  # Strip '*'
    
    # Convert to Polars DataFrame
    df = pl.DataFrame(sequences, schema=["transcript_id", "seq"], orient="row")
    return df

def remove_matches_not_following_RK(peptide_mapping):
    def get_regex_patern(txt):
        return "".join(["[RK]", txt])

    return peptide_mapping\
        .with_columns(
            regex = pl.col("pep").map_elements(lambda x: get_regex_patern(x))
        )\
        .with_columns(
            pl.col("seq").str.extract_all(pl.col("regex"))
        )\
        .explode("seq")\
        .filter(pl.col("seq").is_not_null())\
        .drop("pep", "regex", "seq")

def get_novel_peptide_list(peptide_mapping, gencode_gtf, classification):
    novel_peptides = peptide_mapping\
        .with_columns(
            GENCODE = pl.col("transcript_id").str.starts_with("ENST")
        )\
        .group_by("original_pep")\
        .agg(
            pl.col("GENCODE")
        )\
        .filter(
            pl.col("GENCODE").list.contains(True).not_()
        ).unique("original_pep")["original_pep"].to_list()
    return peptide_mapping\
        .with_columns(
            novel_peptide = pl.col("original_pep").is_in(novel_peptides)
        )\
        .join(
            classification["isoform", "associated_gene"],
            left_on = "transcript_id",
            right_on = "isoform",
            how = "left"
        )\
        .join(
            gencode_gtf.select(["transcript_id", "gene_name"]),
            on = "transcript_id",
            how = "left"
        )\
        .with_columns(
            gene_name = pl.coalesce(pl.col("gene_name"), pl.col("associated_gene"))
        )\
        .drop("associated_gene")\
        .filter(pl.col("novel_peptide"))

def main():
    parser = argparse.ArgumentParser(description='Map peptides to transcripts considering I/L ambiguity')
    parser.add_argument('--annotation_gtf', action='store', type=str, required=True)
    parser.add_argument('--final_sample_classification', action='store', type=str, required=True)
    parser.add_argument('--protein_search_database', action='store', type=str, required=True)
    parser.add_argument('--percolator_res', action='store', type=str, required=True)
    params = parser.parse_args()

    peptide_seq = read_fasta(params.protein_search_database)
    classification = pl.read_parquet(params.final_sample_classification)
    gencode_gtf = read_gtf(params.annotation_gtf, attributes=["gene_name", "transcript_id"])\
        .filter(pl.col("feature")=="transcript")

    percolator_res = pl.read_csv(params.percolator_res, has_header=True, separator="\t")\
            .with_columns(
                proteinIds = pl.col("proteinIds").map_elements(lambda s: s.split(",")[0], return_dtype=pl.String)
            )\
            .with_columns(
                pl.col("peptide").str.replace_all(r"M\[15.9949\]", "M")
            )\
            .with_columns(
                pep = pl.col("peptide").str.split(".").map_elements(lambda x: x[1], return_dtype=pl.String),
                prev_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[0], return_dtype=pl.String),
                next_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[2], return_dtype=pl.String)
            )\
            .unique("pep")\
            .filter(
                pl.col("q-value") < 0.05
            )\
            ["pep"].to_list()

    my_list = []
    for pep in tqdm(percolator_res):
        pep_regex = make_IL_regex(pep)
        df = peptide_seq\
            .with_columns(
                pl.col("seq").str.extract_all(pep_regex).alias("pep")
            )\
            .filter(
                pl.col("pep").list.len() > 0
            )\
            .explode("pep")\
            .with_columns(
                pl.lit(pep).alias("original_pep")
            )\
            .unique(["transcript_id", "original_pep"])
        my_list.append(df)

    peptide_mapping = pl.concat(my_list, how="vertical")

    peptide_mapping_RK = remove_matches_not_following_RK(peptide_mapping)
    novel_peptides_RK = get_novel_peptide_list(peptide_mapping_RK, gencode_gtf, classification)
    
    peptide_mapping.write_parquet("peptide_mapping.parquet")
    novel_peptides_RK.write_csv("novel_peptides.csv")

if __name__ == "__main__":
    main()

# def print_novel_peptide_info(novel_peptides):
#     print(f"""There are {novel_peptides.unique('original_pep').shape[0]} novel peptides mapped to \
# {novel_peptides.unique('transcript_id').shape[0]} novel isoforms \
# that correspond to {novel_peptides.unique('gene_name').shape[0]} genes.""")

# print_novel_peptide_info(novel_peptides_RK)