from src.utils import read_gtf, read_fasta
import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
import random

def get_nucleotide_pos(pep_dict, CDS_coords, strand, query_seq, pb_id):
    """
    Get the nucleotide positions of a peptide sequence in a protein coding gene.
    """
    
    pep_seq = pep_dict[pb_id]
    
    start = pep_seq.find(query_seq)
    if start == -1:
        print("Peptide not found in the sequence")
        return None
    else:
        nucleotide_pos = CDS_coords.filter(pl.col("transcript_id")==pb_id)["coords"][0].to_list()
        if strand == "-":
            end_idx = len(nucleotide_pos) - start*3
            start_idx = len(nucleotide_pos) - (start + len(query_seq)) * 3
            # return Seq("".join([chromosome_seq[i] for i in nucleotide_pos[start_idx:end_idx]])).reverse_complement().translate()
            return nucleotide_pos[start_idx:end_idx]
        else:
            start_idx = start*3
            end_idx = (start + len(query_seq)) * 3
            # return Seq("".join([chromosome_seq[i] for i in nucleotide_pos[start_idx:end_idx]])).translate()
            return nucleotide_pos[start_idx:end_idx]

def sample_substring(s, min_len=1, max_len=None):
    if not s:
        return ''
    if max_len is None:
        max_len = len(s)
    start = random.randint(0, len(s) - min_len)
    end = random.randint(start + min_len, min(start + max_len, len(s)))
    return s[start:end]