from gtfparse import read_gtf
from transcript_transformer.data import parse_transcriptome, save_transcriptome_to_h5
import h5py
import time
data_dict = parse_transcriptome("export/final_transcripts_RiboTIE.gtf", "/Users/xunuo/Genomic_references/GENCODE/GRCh38.primary_assembly.genome.fa")

no_handle = True
max_wait = 900
waited = 0
h5_path = "test.h5"
while no_handle and (waited < max_wait):
    try:
        f = h5py.File(h5_path, "a")
        no_handle = False
    except Exception as e:
        if waited < max_wait:
            time.sleep(120)
            waited += 120
if not no_handle:
    try:
        f = save_transcriptome_to_h5(f, data_dict)
        f.close()
