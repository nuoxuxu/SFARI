import pyBigWig
from Bio import SeqIO
import polars as pl
from src.utils import gtf_to_SJ, read_gtf, read_SJ
import os
import numpy as np
import polars.selectors as cs
from pathlib import Path

bw = "data/hg38.phyloP100way.bw"
pbw = pyBigWig.open(bw)

genome = list(SeqIO.parse(os.getenv("GENOMIC_DATA_DIR") + "/GENCODE/GRCh38.primary_assembly.genome.fa", "fasta")) # type: ignore

def phylop(row):
    """ Get the phyloP score for a given variant.
    """
    return pbw.values(row["chrom"], row["pos"]-1, row["pos"])[0]

def add_phylop_to_df(df):
    acceptor_raw_phyloP= df\
        .map_rows(lambda row: tuple(pbw.values(row[1], row[3]-1-25, row[3]+10)))\
        .drop_nans()
    
    acceptor_sem = acceptor_raw_phyloP\
        .select(
            cs.all().std()/np.sqrt(df.shape[0])
        )\
        .transpose(column_names=["sem"]).to_series()
    
    Acceptor = acceptor_raw_phyloP\
        .mean()\
        .transpose(column_names=["phyloP"])\
        .hstack([acceptor_sem, pl.Series("pos", range(-25, 11))])\
        .with_columns(
            region = pl.lit("Acceptor")
        )

    donor_raw_phyloP = df\
        .map_rows(lambda row: tuple(pbw.values(row[1], row[2]-1-10, row[2]+10)))\
        .drop_nans()
    
    donor_sem = donor_raw_phyloP\
        .select(
            cs.all().std()/np.sqrt(df.shape[0])
        )\
        .transpose(column_names=["sem"]).to_series()

    Donor = donor_raw_phyloP\
        .mean()\
        .transpose(column_names=["phyloP"])\
        .hstack([donor_sem, pl.Series("pos", range(-10, 11))])\
        .with_columns(
            region = pl.lit("Donor")
        )

    export = pl.concat([Acceptor, Donor], how="vertical")
    return export

def get_acceptor_seq(df):
    """ Get the acceptor sequence for given splice junctions.
    """
    return df\
        .map_rows(lambda row: str(genome[np.nonzero([record.id == row[1] for record in genome])[0].tolist()[0]].seq[row[3]-2:row[3]]))\
        .rename({"map": "acceptor_seq"})

def get_donor_seq(df):
    """ Get the donor sequence for given splice junctions.
    """
    return df\
        .map_rows(lambda row: str(genome[np.nonzero([record.id == row[1] for record in genome])[0].tolist()[0]].seq[row[2]-1:row[2]+1]))\
        .rename({"map": "donor_seq"})

def filter_for_canonical(df):
    """ Filter for canonical splice sites.
    """
    return pl.concat([df, get_acceptor_seq(df), get_donor_seq(df)], how="horizontal")\
        .filter(
            (pl.col("acceptor_seq")=="AG") & (pl.col("donor_seq")=="GT")
        )

def filter_for_riboseq_evidence(df):
    riboseq_SJ = pl.concat([read_SJ(file) for file in Path('data/riboseq').rglob('*_SJ.out.tab')], how='vertical')\
        .unique(['chrom', 'start', 'end', 'strand'])\
        .select(['chrom', 'start', 'end', 'strand'])\
        .with_columns(
            strand = pl.col('strand').map_elements(lambda s: '+' if s == 1 else '-', return_dtype=pl.String)
        )

    return df\
        .join(
            riboseq_SJ.drop("strand"),
            on=["chrom", "start", "end"],
            how="inner"
        )

def filter_for_peptide_evidence(df):
    annot_peptides_hybrid = read_gtf("nextflow_results/V47/orfanage/annot_peptides_hybrid.gtf", attributes=["detected", "type", "transcript_id"])
    annot_peptides_hybrid = annot_peptides_hybrid\
        .filter(
            (pl.col("detected") == "True") & (pl.col("type") == "splice-junction")
        )
    peptide_SJ = gtf_to_SJ(annot_peptides_hybrid)

    return df\
        .join(
            peptide_SJ.drop("transcript_id", "strand"),
            on=["chrom", "start", "end"],
            how="inner"
        )

def filter_combined(df, canonical_ss=False, riboseq_evidence = False, translational_evidence=False):
    """ Filter the splice junctions for canonical splice sites and/or translational evidence.
    """
    if canonical_ss and riboseq_evidence and translational_evidence:
        return filter_for_peptide_evidence(filter_for_riboseq_evidence(filter_for_canonical(df)))
    if canonical_ss and riboseq_evidence:
        return filter_for_riboseq_evidence(filter_for_canonical(df))
    if canonical_ss and translational_evidence:
        return filter_for_peptide_evidence(filter_for_canonical(df))    
    elif riboseq_evidence and translational_evidence:
        return filter_for_peptide_evidence(filter_for_riboseq_evidence(df))
    elif canonical_ss:
        return filter_for_canonical(df)
    elif translational_evidence:
        return filter_for_peptide_evidence(df)
    elif riboseq_evidence:
        return filter_for_riboseq_evidence(df)
    else:
        return lambda df: df

def export_phyloP(feature, out, canonical_ss=False, riboseq_evidence=False, translational_evidence=False):
    """ Export the phyloP scores to a CSV file.
    """
    final_transcripts_SJ = read_gtf("nextflow_results/V47/orfanage/orfanage.gtf")\
        .filter(pl.col("feature") == feature)\
        .pipe(gtf_to_SJ)\
        .select("transcript_id", "chrom", "start", "end")

    gencode_gtf = os.getenv("GENOMIC_DATA_DIR") + "/GENCODE/gencode.v47.annotation.gtf" # type: ignore
    GENCODE_SJ = read_gtf(gencode_gtf)\
        .filter(pl.col("feature") == feature)\
        .pipe(gtf_to_SJ)\
        .select("transcript_id", "chrom", "start", "end")
    
    
    known_SJ_ss = GENCODE_SJ\
        .unique(["chrom", "start", "end"])\
        .pipe(filter_combined, canonical_ss=canonical_ss, translational_evidence=translational_evidence, riboseq_evidence= riboseq_evidence)\
        .pipe(add_phylop_to_df)\
        .with_columns(
            spl_type = pl.lit("known")
        )

    novel_3prime_ss = final_transcripts_SJ\
        .join(
                GENCODE_SJ,
                on=["chrom", "start"]
            )\
        .filter(
            pl.col("end").is_in(GENCODE_SJ["end"]).not_()
        ).\
        unique(["chrom", "start", "end"])\
        .pipe(filter_combined, canonical_ss=canonical_ss, translational_evidence=translational_evidence, riboseq_evidence= riboseq_evidence)\
        .pipe(add_phylop_to_df)\
        .with_columns(
            spl_type = pl.lit("novel_3prime")
        )

    novel_5prime_ss = final_transcripts_SJ\
        .join(
                GENCODE_SJ,
                on=["chrom", "end"]
            )\
        .filter(
            pl.col("start").is_in(GENCODE_SJ["start"]).not_()
        )\
        .unique(["chrom", "start", "end"])\
        .pipe(filter_combined, canonical_ss=canonical_ss, translational_evidence=translational_evidence, riboseq_evidence= riboseq_evidence)\
        .pipe(add_phylop_to_df)\
        .with_columns(
            spl_type = pl.lit("novel_5prime")
        )

    novel_both_ss = final_transcripts_SJ\
        .filter(
            pl.col("start").is_in(GENCODE_SJ["start"]).not_() &
            pl.col("end").is_in(GENCODE_SJ["end"]).not_()
        )\
        .unique(["chrom", "start", "end"])\
        .pipe(filter_combined, canonical_ss=canonical_ss, translational_evidence=translational_evidence, riboseq_evidence= riboseq_evidence)\
        .pipe(add_phylop_to_df)\
        .with_columns(
            spl_type = pl.lit("novel_both")
        )
        
    export = pl.concat([known_SJ_ss, novel_3prime_ss, novel_5prime_ss, novel_both_ss], how="vertical")
    export.write_csv(out)

export_phyloP("exon", "export/variant/exon_ss_phyloP.csv")
export_phyloP("CDS", "export/variant/CDS_ss_phyloP.csv")
export_phyloP("exon", "export/variant/exon_ss_phyloP_canonical.csv", canonical_ss=True)
export_phyloP("CDS", "export/variant/CDS_ss_phyloP_canonical.csv", canonical_ss=True)
export_phyloP("exon", "export/variant/exon_ss_phyloP_canonical_translational.csv", canonical_ss=True, translational_evidence=True)
export_phyloP("CDS", "export/variant/CDS_ss_phyloP_canonical_translational.csv", canonical_ss=True, translational_evidence=True)
export_phyloP("exon", "export/variant/exon_ss_phyloP_canonical_riboseq.csv", canonical_ss=True, riboseq_evidence=True)
export_phyloP("CDS", "export/variant/CDS_ss_phyloP_canonical_riboseq.csv", canonical_ss=True, riboseq_evidence=True)