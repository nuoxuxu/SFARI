from src.utils import read_gtf, read_fasta
import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
import random

records = list(SeqIO.parse("/project/s/shreejoy/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa", "fasta"))

orfanage_CDS = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")\
    .filter(feature="CDS")

pep_dict = read_fasta("nextflow_results/V47/orfanage/orfanage_peptide_formatted.fasta")\
    .to_pandas()\
    .set_index("transcript_id")\
    .to_dict()["seq"]

CDS_coords = orfanage_CDS\
    .sort("transcript_id", "start")\
    .with_columns(
        pl.struct(["start", "end"]).
        map_elements(lambda s: range(s["start"]-1, s["end"]), return_dtype=pl.List(pl.Int64)).alias("coords")
    )\
    .group_by("transcript_id")\
    .agg(
        pl.col("coords").flatten()
    )

def get_nucleotide_pos(query_seq, pb_id):
    """
    Get the nucleotide positions of a peptide sequence in a protein coding gene.
    """
    # Find the start and end indices of the peptide sequence in the protein coding gene
    pep_seq = pep_dict[pb_id]
    chromosome = orfanage_CDS\
        .filter(pl.col("transcript_id")==pb_id)\
        .select("seqname").unique()["seqname"][0]
    strand = orfanage_CDS\
        .filter(pl.col("transcript_id")==pb_id)\
        .select("strand").unique()["strand"][0]
    
    start = pep_seq.find(query_seq)
    if start == -1:
        print("Peptide not found in the sequence")
        return None
    else:
        nucleotide_pos = CDS_coords.filter(pl.col("transcript_id")==pb_id)["coords"][0].to_list()
        chromosome_seq = [record for record in records if record.id == chromosome][0].seq
        if strand == "-":
            end_idx = len(nucleotide_pos) - start*3
            start_idx = len(nucleotide_pos) - (start + len(query_seq)) * 3
            # return Seq("".join([chromosome_seq[i] for i in nucleotide_pos[start_idx:end_idx]])).reverse_complement().translate()
            return start_idx, end_idx
        else:
            start_idx = start*3
            end_idx = (start + len(query_seq)) * 3
            # return Seq("".join([chromosome_seq[i] for i in nucleotide_pos[start_idx:end_idx]])).translate()
            return start_idx, end_idx

def sample_substring(s, min_len=1, max_len=None):
    if not s:
        return ''
    if max_len is None:
        max_len = len(s)
    start = random.randint(0, len(s) - min_len)
    end = random.randint(start + min_len, min(start + max_len, len(s)))
    return s[start:end]

def main():
    tx_to_be_tested = orfanage_CDS.select("transcript_id").unique().to_series().to_list()
    seq_to_be_tested = [sample_substring(pep_dict[tx], 10) for tx in tx_to_be_tested]

    [print(get_nucleotide_pos(seq_to_be_tested[i], tx_to_be_tested[i]), tx_to_be_tested[i], seq_to_be_tested[i]) for i in range(10)]

if __name__ == "__main__":
    main()