#!/usr/bin/env python
import pysam, numpy as np, sys
from pathlib import Path

bam_in     = sys.argv[1]          # collapsed Iso-Seq BAM
timepoint = Path(bam_in).name.removesuffix(".bam")  # e.g. "iPSC_merged"
bedgraph = f"data/katherine/{timepoint}.coverage.weighted.bedgraph"

# Pre-allocate one array per contig
sam = pysam.AlignmentFile(bam_in)
cov = {r: np.zeros(l, dtype=np.int32)
       for r,l in zip(sam.references, sam.lengths)}

for rec in sam.fetch(until_eof=True):
    # 1) how many raw HiFi reads built this isoform?
    try:
        n = rec.get_tag("is")
    except KeyError:
        n = len(rec.get_tag("im").split(","))

    # 2) add that count to every reference block in the CIGAR
    for start, end in rec.get_blocks():      # gets exon blocks, skips N introns
        cov[rec.reference_name][start:end] += n

sam.close()

# 3) emit BedGraph (only non-zero stretches)
with open(bedgraph, "w") as out:
    for chrom, arr in cov.items():
        run_start, current = None, 0
        for pos, depth in enumerate(arr):
            if depth != current:
                if current:                   # close previous run
                    out.write(f"{chrom}\t{run_start}\t{pos}\t{current}\n")
                run_start, current = (pos, depth) if depth else (None, 0)
        if current:                           # flush last run
            out.write(f"{chrom}\t{run_start}\t{len(arr)}\t{current}\n")