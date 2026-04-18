#!/usr/bin/env python3
"""
Compute overlap between Xu and Rynard et al. transcripts and ENCODE4
brain-expressed transcripts by intron chain identity.

Outputs encode4_overlap.parquet with columns:
  isoform     - SFARI transcript ID
  type        - "known" (FSM) or "novel" (all other categories)
  in_encode4  - True if the intron chain is present in ENCODE4 brain-expressed transcripts
"""
import argparse
import re
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from collections import defaultdict


def get_intron_chain_keys(gtf_path):
    tx_exons = defaultdict(list)
    tid_re = re.compile(r'transcript_id "([^"]+)"')
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'exon':
                continue
            m = tid_re.search(fields[8])
            if not m:
                continue
            tx_exons[m.group(1)].append((fields[0], fields[6], int(fields[3]), int(fields[4])))
    result = {}
    for tid, exons in tx_exons.items():
        chrom, strand = exons[0][0], exons[0][1]
        coords = sorted(exons, key=lambda x: x[2])
        if len(coords) == 1:
            result[tid] = f"{chrom}:{strand}:{coords[0][2]}-{coords[0][3]}:mono"
        else:
            result[tid] = ','.join(
                f"{coords[i][3]}-{coords[i+1][2]}" for i in range(len(coords) - 1)
            )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sfari_classification', required=True)
    parser.add_argument('--sfari_gtf',            required=True)
    parser.add_argument('--encode4_gtf',           required=True)
    parser.add_argument('--encode4_abundance',     required=True)
    args = parser.parse_args()

    # ENCODE4: apply expression filter, build intron chain set
    enc = pd.read_csv(args.encode4_abundance, sep='\t', low_memory=False)
    brodmann_cols = [c for c in enc.columns if c.startswith('brodmann')]
    enc_filtered  = enc[(enc[brodmann_cols] > 5).sum(axis=1) > 2]
    print(f"ENCODE4: {len(enc_filtered):,} transcripts pass expression filter")

    enc_keys   = get_intron_chain_keys(args.encode4_gtf)
    enc_chains = {enc_keys[t] for t in enc_filtered['annot_transcript_id'] if t in enc_keys}
    print(f"ENCODE4: {len(enc_chains):,} unique intron chains")

    # SFARI: intron chains and classification
    sfari_keys = get_intron_chain_keys(args.sfari_gtf)
    sfari_cls  = pq.read_table(args.sfari_classification).to_pandas()[['isoform', 'structural_category']]

    sfari_cls['type']      = sfari_cls['structural_category'].map(
        lambda c: 'known' if c == 'full-splice_match' else 'novel'
    )
    sfari_cls['in_encode4'] = sfari_cls['isoform'].map(
        lambda t: sfari_keys.get(t) in enc_chains if sfari_keys.get(t) is not None else False
    )

    out = sfari_cls[['isoform', 'type', 'in_encode4']]
    pq.write_table(
        pa.Table.from_pandas(out),
        "encode4_overlap.parquet",
        compression='snappy'
    )
    print(f"\nWrote encode4_overlap.parquet ({len(out):,} rows)")
    print(out.groupby(['type', 'in_encode4']).size().to_string())


if __name__ == '__main__':
    main()
