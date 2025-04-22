from src.utils import read_gtf
import polars as pl
from tqdm import tqdm
import numpy as np
import multiprocessing as mp
from src.utils import read_CAGE_peaks, read_polyA_peaks

#--------Read in datasets--------
cage_peak_obj = read_CAGE_peaks("data/refTSS_v3.3_human_coordinate.hg38.sorted.bed")
polya_peak_obj = read_polyA_peaks("data/atlas.clusters.2.0.GRCh38.96.bed")

gtf = read_gtf("nextflow_results/V47/final_transcripts.gtf")\
    .filter(pl.col("feature")=="transcript")\
    .rename({"seqname": "chrom", "start": "txStart", "end": "txEnd"})\
    .select(["transcript_id", "chrom", "txStart", "txEnd", "strand"]).to_pandas()

classification = pl.read_parquet("nextflow_results/V47/final_classification.parquet")

#--------CAGE peaks--------
def process_cage_peak_info(pbid):
    rec = gtf.loc[gtf["transcript_id"] == pbid, :].transpose().iloc[:, 0].to_dict()
    if rec["strand"] == '+':
        within_CAGE, dist_CAGE = cage_peak_obj.find(rec["chrom"], rec["strand"], rec["txStart"])
    else:
        within_CAGE, dist_CAGE = cage_peak_obj.find(rec["chrom"], rec["strand"], rec["txEnd"])
    return within_CAGE, dist_CAGE

with mp.Pool(processes=30) as pool:
    within_CAGE_peak = []
    dist_to_CAGE_peak = []
    for result in tqdm(pool.imap(process_cage_peak_info, classification["isoform"]), total=len(classification["isoform"])):
        within_CAGE, dist_CAGE = result
        within_CAGE_peak.append(within_CAGE)
        dist_to_CAGE_peak.append(dist_CAGE)

dist_to_CAGE_peak = np.array([float(x) if x != 'NA' else np.nan for x in dist_to_CAGE_peak], dtype=np.float64)
within_CAGE_peak = np.array([x == "TRUE" for x in within_CAGE_peak], dtype=bool)

classification = classification\
    .with_columns(
        dist_to_CAGE_peak = pl.Series(dist_to_CAGE_peak),
        within_CAGE_peak = pl.Series(within_CAGE_peak, dtype=pl.Boolean)
    )

#--------polyA peaks--------
def process_polya_peak_info(pbid):
    rec = gtf.loc[gtf["transcript_id"] == pbid, :].transpose().iloc[:, 0].to_dict()
    if rec.strand == '+':
        within_polyA_site, dist_polyA_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txEnd)
    else:
        within_polyA_site, dist_polyA_site = polya_peak_obj.find(rec.chrom, rec.strand, rec.txStart)
    return within_polyA_site, dist_polyA_site

with mp.Pool(processes=30) as pool:
    within_polyA_site = []
    dist_to_polyA_site = []
    for result in tqdm(pool.imap(process_polya_peak_info, classification["isoform"]), total=len(classification["isoform"])):
        within_polyA_site, dist_polyA_site = result
        within_polyA_site.append(within_polyA_site)
        dist_to_polyA_site.append(dist_polyA_site)

dist_to_polyA_site = np.array([float(x) if x != 'NA' else np.nan for x in dist_to_polyA_site], dtype=np.float64)
within_polyA_site = np.array([x == "TRUE" for x in within_polyA_site], dtype=bool)

classification = classification\
    .with_columns(
        dist_to_polyA_site = pl.Series(dist_to_polyA_site),
        within_polyA_site = pl.Series(within_polyA_site, dtype=pl.Boolean)
    )

#--------SQANTI3_qc--------
SQANTI3_qc = pl.read_csv("SQANTI3_qc_classification.txt", separator="\t", null_values=["NA"])\
    .filter(
        pl.col("isoform").is_in(classification["isoform"])
    )

SQANTI3_qc.filter(
        pl.col("dist_to_CAGE_peak").abs() < 100
    )

classification.filter(
        pl.col("dist_to_CAGE_peak").abs() < 100
    )

classification["dist_to_CAGE_peak"]