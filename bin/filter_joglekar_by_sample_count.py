#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict


def intron_chain_key(fields):
    """
    Group reads by splice junctions only, ignoring variable TSS/TTS.
    Multi-exon: 'chrom:strand:e1end-e2start,e2end-e3start,...'
    Mono-exon:  'chrom:strand:start-end:mono'
    """
    chrom       = fields[0]
    chrom_start = int(fields[1])
    strand      = fields[5]
    block_count = int(fields[9])
    sizes  = [int(x) for x in fields[10].rstrip(',').split(',') if x]
    starts = [int(x) for x in fields[11].rstrip(',').split(',') if x]
    exons = []
    for i in range(block_count):
        s = chrom_start + starts[i]
        exons.append((s, s + sizes[i]))
    if len(exons) == 1:
        return f"{chrom}:{strand}:{exons[0][0]}-{exons[0][1]}:mono"
    introns = ','.join(f"{exons[i][1]}-{exons[i+1][0]}" for i in range(len(exons) - 1))
    return f"{chrom}:{strand}:{introns}"


def main():
    parser = argparse.ArgumentParser(
        description="Keep reads whose intron chain appears in more than N samples."
    )
    parser.add_argument("min_samples", type=int,
                        help="Minimum number of samples a transcript must appear in")
    parser.add_argument("output", help="Output BED file path")
    args = parser.parse_args()

    bed_files = sorted(f for f in os.listdir('.') if f.endswith('_corrected_reads.bed'))

    # Count how many samples contain each intron chain
    chain_samples = defaultdict(set)
    for bed_file in bed_files:
        sample = bed_file.replace('_corrected_reads.bed', '')
        with open(bed_file) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) < 12:
                    continue
                chain_samples[intron_chain_key(fields)].add(sample)

    keep = {k for k, s in chain_samples.items() if len(s) > args.min_samples}

    # Write a single combined BED — one representative read per passing intron chain
    written_keys = set()
    with open(args.output, 'w') as out:
        for bed_file in bed_files:
            with open(bed_file) as fh:
                for line in fh:
                    if line.startswith('#'):
                        continue
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 12:
                        continue
                    k = intron_chain_key(fields)
                    if k in keep and k not in written_keys:
                        out.write(line)
                        written_keys.add(k)


if __name__ == "__main__":
    main()
