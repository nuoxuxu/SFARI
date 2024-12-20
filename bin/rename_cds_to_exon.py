#!/home/s/shreejoy/nxu/miniforge3/envs/patch_seq_spl/bin/python

import polars as pl
from src.utils import read_gtf
from src.single_cell import SingleCell

sample_gtf = read_gtf("nextflow_results/transcripts_filtered.fasta.transdecoder.genome.gff3.gtf")