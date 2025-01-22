#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python
import polars as pl
from collections import defaultdict, namedtuple
from Bio import SeqIO
from src.single_cell import SingleCell
import pandas as pd
import copy
import argparse
from src.utils import read_gtf

def read_fasta(fasta_file):
    seqs = defaultdict()
    for rec in SeqIO.parse(fasta_file, 'fasta'):
        if rec.id.startswith("ENSP"):            
            seqs[rec.id.split("|")[1]] = str(rec.seq)
        elif rec.id.count("p")==1:
            seqs[rec.id.rsplit(".", 1)[0]] = str(rec.seq)
        else:
            seqs[rec.id] = str(rec.seq)
    return seqs

def find_start_pep_index(row):
    pep = row[0]
    pb_acc = row[1]
    try:
        seq = seqs[pb_acc]
    except:
        return 0
    start = seq.find(pep)
    if start == -1:
        return 0
    start_idx = start + 1
    return start_idx

def find_end_pep_index(row):
    pep = row[0]
    pb_acc = row[1]
    try:
        seq = seqs[pb_acc]
    except:
        return 0
    start = seq.find(pep)
    if start == -1:
        return 0
    start_idx = start + 1
    end_idx = start_idx + len(pep) - 1
    return end_idx

def make_cumulative_blens(blocks):
    cblocks = []
    cbl = 0 # cumulative block length
    for b in blocks:
        cbl += b
        cblocks.append(cbl)
    return cblocks

def process_gtf(gtf):
    # CDS coords into dict
    pbs = defaultdict(lambda: ['chr', 'strand', [], [], []]) # pb -> [chr, strand, [start, end], [block lengths], [cum. block lengths]]
    # PB.1.1 -> ['chr1', '+', [[100,150], [200,270]], [50, 70], [50, 120], [150-200]]
    for i, row in gtf.iterrows():
        chr, feat, start, end, strand, acc = row
        pbs[acc][0] = chr
        pbs[acc][1] = strand
        pbs[acc][2].append([int(start), int(end)])
    
    for acc, infos in pbs.items():
        strand = infos[1]
        if strand == '+':
            infos[2] = sorted(infos[2])
        elif strand == '-':
            infos[2] = sorted(infos[2], reverse=True)
        infos[3] = [end-start+1 for [start, end] in infos[2]]
        infos[4] = make_cumulative_blens(infos[3])
    return pbs

def read_reference_gtf(gtf_file):
    gtf = pd.read_table(gtf_file, skiprows=5, header=None)
    gtf = gtf[[0, 2, 3, 4, 6, 8]]
    gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']
    gtf = gtf.loc[gtf['feat']=='CDS']

    gtf['acc'] = gtf['acc'].str.split('transcript_id "').str[1]
    gtf['acc'] = gtf['acc'].str.split('";').str[0]
    return gtf

def read_sample_gtf(gtf_file):
    gtf = pd.read_table(gtf_file, skiprows=1, header=None)
    gtf = gtf[[0, 2, 3, 4, 6, 8]]
    gtf.columns = ['chr', 'feat', 'start', 'end', 'strand', 'acc']
    gtf = gtf.loc[gtf['feat']=='CDS']
    gtf['acc'] = gtf['acc'].str.extract(r'transcript_id "([^;]*)";')
    return gtf

def get_first_block_index(orf_coord, cblens):
    # get the index corresponding to the first block containing the orf start
    # return index, and the dela (spacing upstream of end)
    for i, cblen in enumerate(cblens):
        if orf_coord <= cblen:
            delta = cblen - orf_coord
            return i, delta
    return i, 0

def most_frequent(List):
	counter = 0
	element = List[0]
	for i in List:
		curr_frequency = List.count(i)
		if(curr_frequency > counter):
			counter = curr_frequency
			element = i
	return element

def make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1: i2+1]
    # trim ends to orf start/end
    orf_coords[0][0] = orf_coords[0][1] - delta1
    orf_coords[-1][1] = orf_coords[-1][1] - delta2
    return orf_coords

def make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords):
    orf_coords = copy.deepcopy(coords)
    orf_coords = orf_coords[i1: i2+1]
    # trim ends to orf start/end
    orf_coords[0][1] = orf_coords[0][0] + delta1
    orf_coords[-1][0] = orf_coords[-1][0] + delta2
    return orf_coords

def write_peptide_gtf(output_name, pep_ranges, pbs):
    # Note - code below was modified from script that found CDS orf range wtihin
    # a full-length RNA isoforms, so orf-range is now peptide range 
    with open(output_name, 'w') as ofile:
        # write UCSC track header
        # remove for conversion to bed12 (genePred complains)
        # ofile.write('track name=peptide color=0,0,0\n')
        for i, row in pep_ranges.iterrows():
            pep_seq, pb_acc, prev_aa, next_aa, gene, PSMId, pep_start, pep_end = row
            # convert from protein (AA) to CDS (nt) coords
            pep_start = pep_start * 3 - 2
            pep_end = pep_end * 3
            if pb_acc in pbs:
                infos = pbs[pb_acc]
                chr, strand, coords, blens, cblens = infos
                i1, delta1 = get_first_block_index(pep_start, cblens)
                i2, delta2 = get_first_block_index(pep_end, cblens)
                if strand == '+':
                    orf_coords = make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords)
                elif strand == '-':
                    orf_coords = make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords)
                # write out the coordinates
                prev_aa = most_frequent(prev_aa.split('|'))
                next_aa = most_frequent(next_aa.split('|'))
                if chr in ['chrX','chrY']:
                    gene = f"{gene}_{chr}"
                acc_id= f"{prev_aa}.{pep_seq}.{next_aa}({gene})"
                pep_acc = f'gene_id "{PSMId}"; gene_name "{gene}"; transcript_id "{acc_id}";'
                for [start, end] in orf_coords:
                    ofile.write('\t'.join([chr, 'hg38_canon', 'exon', str(start), str(end), '.', strand,
                                '.', pep_acc]) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter genome_gffs from TransDecoder with novel transripts')
    parser.add_argument('--filter', action='store_true')
    parser.add_argument('--sample_gtf', action='store', type=str, required=True)
    parser.add_argument('--reference_gtf', action='store', type=str, required=True)
    parser.add_argument('--peptides', action='store', type=str, required=True)
    parser.add_argument('--h5ad_file', action='store', type=str, required=True)
    parser.add_argument('--sample_fasta', action='store', type=str, required=True)
    parser.add_argument('--output', action='store', type=str, required=True)
    params = parser.parse_args()
    
    seqs = read_fasta(params.sample_fasta)
    pb_gene = SingleCell(params.h5ad_file).var["isoform", "associated_gene"].with_columns(pl.col("associated_gene").cast(pl.String))
    gencode_gene = read_gtf("/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.annotation.gtf", attributes = ["gene_name", "transcript_id"])\
        .filter(
            pl.col("feature") == "transcript",
            pl.col("transcript_id").is_in(seqs.keys())
        )["transcript_id", "gene_name"].rename({"transcript_id": "isoform", "gene_name": "associated_gene"})
    pb_gene = pl.concat([pb_gene, gencode_gene])

    if params.filter:
        percolator_res = pl.read_csv(params.peptides, has_header=True, separator="\t")\
            .with_columns(
            proteinIds = pl.col("proteinIds").map_elements(lambda s: s.split(","), return_dtype=pl.List(pl.String))
            )\
            .explode("proteinIds")\
            .filter(
                pl.col("q-value") < 0.05,
                pl.col("proteinIds").str.contains("DECOY").not_()
            )\
            .with_columns(
                pl.when(pl.col("proteinIds").str.starts_with("ENSP"))\
                    .then(pl.col("proteinIds").str.extract(r"^[^|]*\|([^|]*)"))\
                    .otherwise(pl.col("proteinIds"))
            )\
            .with_columns(
                pl.col("peptide").str.replace_all(r"M\[15.9949\]", "M")
            )\
            .with_columns(
                pep = pl.col("peptide").str.split(".").map_elements(lambda x: x[1], return_dtype=pl.String),
                prev_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[0], return_dtype=pl.String),
                next_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[2], return_dtype=pl.String)
            )\
            .rename(
                {"proteinIds": "pb_acc"}
            )\
            .join(pb_gene.rename({"isoform": "pb_acc"}), on = "pb_acc", how = "left")\
            .rename({"associated_gene": "gene"})[['pep', 'pb_acc', 'prev_aa','next_aa', 'gene', 'PSMId']]\
            .unique("pep")
    else:
        percolator_res = pl.read_csv(params.peptides, has_header=True, separator="\t")\
            .with_columns(
            proteinIds = pl.col("proteinIds").map_elements(lambda s: s.split(","), return_dtype=pl.List(pl.String))
            )\
            .explode("proteinIds")\
            .with_columns(
                pl.when(pl.col("proteinIds").str.starts_with("ENSP"))\
                    .then(pl.col("proteinIds").str.extract(r"^[^|]*\|([^|]*)"))\
                    .otherwise(pl.col("proteinIds"))
            )\
            .with_columns(
                pl.col("peptide").str.replace_all(r"M\[15.9949\]", "M")
            )\
            .with_columns(
                pep = pl.col("peptide").str.split(".").map_elements(lambda x: x[1], return_dtype=pl.String),
                prev_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[0], return_dtype=pl.String),
                next_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[2], return_dtype=pl.String)
            )\
            .rename(
                {"proteinIds": "pb_acc"}
            )\
            .join(pb_gene.rename({"isoform": "pb_acc"}), on = "pb_acc", how = "left")\
            .rename({"associated_gene": "gene"})[['pep', 'pb_acc', 'prev_aa','next_aa', 'gene', 'PSMId']]\
            .unique("pep")        

    start_idx = percolator_res.map_rows(find_start_pep_index).rename({"map": "pep_start"})
    end_idx = percolator_res.map_rows(find_end_pep_index).rename({"map": "pep_end"})

    pep_ranges = pl.concat([pl.concat([percolator_res, start_idx], how="horizontal"), end_idx], how="horizontal")\
        .filter(
            pl.col("pb_acc").is_in(seqs.keys())
        )\
        .filter(
            ((pl.col("pep_start")<0).or_(pl.col("pep_end")<0)).not_()
        )
    
    sample_gtf = read_sample_gtf(params.sample_gtf)
    reference_gtf = read_reference_gtf(params.reference_gtf)
    pbs = process_gtf(pd.concat([sample_gtf, reference_gtf], ignore_index=True))

    write_peptide_gtf(params.output, pep_ranges.to_pandas(), pbs)
        
else:
    Args = namedtuple("Args", "sample_gtf reference_gtf peptides h5ad_file sample_fasta output")
    params = Args(
        sample_gtf = "nextflow_results/transcripts_filtered.fasta.transdecoder.genome.gff3.gtf",
        reference_gtf = "/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.annotation.gtf",
        peptides = "nextflow_results/hybrid_percolator.tsv",
        h5ad_file = "nextflow_results/pbid_orf.h5ad",
        sample_fasta = "nextflow_results/hybrid.fasta",
        output = "test_peptides.gtf"
    )
    seqs = read_fasta(params.sample_fasta)
    pb_gene = SingleCell(params.h5ad_file).var["isoform", "associated_gene"].with_columns(pl.col("associated_gene").cast(pl.String))
    gencode_gene = read_gtf("/project/s/shreejoy/Genomic_references/GENCODE/gencode.v39.annotation.gtf", attributes = ["gene_name", "transcript_id"])\
        .filter(
            pl.col("feature") == "transcript",
            pl.col("transcript_id").is_in(seqs.keys())
        )["transcript_id", "gene_name"].rename({"transcript_id": "isoform", "gene_name": "associated_gene"})
    pb_gene = pl.concat([pb_gene, gencode_gene])
    percolator_res = pl.read_csv(params.peptides, has_header=True, separator="\t")\
        .with_columns(
        proteinIds = pl.col("proteinIds").map_elements(lambda s: s.split(","), return_dtype=pl.List(pl.String))
        )\
        .explode("proteinIds")\
        .filter(
            pl.col("q-value") < 0.05,
            pl.col("proteinIds").str.contains("DECOY").not_()
        )\
        .with_columns(
            pl.when(pl.col("proteinIds").str.starts_with("ENSP"))\
                .then(pl.col("proteinIds").str.extract(r"^[^|]*\|([^|]*)"))\
                .otherwise(pl.col("proteinIds"))
        )\
        .with_columns(
            pl.col("peptide").str.replace_all(r"M\[15.9949\]", "M")
        )\
        .with_columns(
            pep = pl.col("peptide").str.split(".").map_elements(lambda x: x[1], return_dtype=pl.String),
            prev_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[0], return_dtype=pl.String),
            next_aa = pl.col("peptide").str.split(".").map_elements(lambda x: x[2], return_dtype=pl.String)
        )\
        .rename(
            {"proteinIds": "pb_acc"}
        )\
        .join(pb_gene.rename({"isoform": "pb_acc"}), on = "pb_acc", how = "left")\
        .rename({"associated_gene": "gene"})[['pep', 'pb_acc', 'prev_aa','next_aa', 'gene']]\
        .unique("pep")

    start_idx = percolator_res.map_rows(find_start_pep_index).rename({"map": "pep_start"})
    end_idx = percolator_res.map_rows(find_end_pep_index).rename({"map": "pep_end"})

    pep_ranges = pl.concat([pl.concat([percolator_res, start_idx], how="horizontal"), end_idx], how="horizontal")\
        .with_columns(
            pep_start = pl.col("pep_start") * 3 - 2,
            pep_end = pl.col("pep_end") * 3
        )\
        .filter(
            pl.col("pb_acc").is_in(seqs.keys())
        )\
        .filter(
            ((pl.col("pep_start")<0).or_(pl.col("pep_end")<0)).not_()
        )
        
    cds = read_fasta("nextflow_results/transcripts_filtered.fasta.transdecoder.cds")
    transcript = read_fasta("nextflow_results/transcripts_filtered.fasta")

    sample_gtf = read_sample_gtf(params.sample_gtf)
    reference_gtf = read_reference_gtf(params.reference_gtf)
    pbs = process_gtf(pd.concat([sample_gtf, reference_gtf], ignore_index=True))

    chr, strand, coords, blens, cblens = pbs["PB.100609.85"]
    pep_seq, pb_acc, prev_aa, next_aa,gene, pep_start, pep_end = pep_ranges.filter(pl.col("pep").str.contains("TLTAEEAEEEWERR"))

    i1, delta1 = get_first_block_index(pep_start[0], cblens)
    i2, delta2 = get_first_block_index(pep_end[0], cblens)

    orf_coords = make_coords_trimmed_to_orf_range_neg_strand(i1, delta1, i2, delta2, coords)
    make_coords_trimmed_to_orf_range(i1, delta1, i2, delta2, coords)

    seqs["PB.100609.85"].find("TLTAEEAEEEWERR")