#!/usr/bin/env python3
"""
Make a 2-way Venn diagram comparing fusion transcripts between
SFARI (Xu and Rynard et al.) and Patowary et al. datasets.

Transcripts are compared by intron chain:
  multi-exon : 'chrom:strand:e1end-e2start,e2end-e3start,...'
  mono-exon  : 'chrom:strand:start-end:mono'

Only transcripts with structural_category == 'fusion' are included.
"""
import argparse
import re
from collections import defaultdict

import pandas as pd
import pyarrow.parquet as pq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2


# ---------------------------------------------------------------------------
# GTF helpers
# ---------------------------------------------------------------------------

def get_intron_chain_keys(gtf_path):
    """
    Return {transcript_id: key} where key is:
      multi-exon : 'chrom:strand:e1end-e2start,e2end-e3start,...'
      mono-exon  : 'chrom:strand:start-end:mono'
    """
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
            tid = m.group(1)
            tx_exons[tid].append((fields[0], fields[6], int(fields[3]), int(fields[4])))

    result = {}
    for tid, exons in tx_exons.items():
        chrom, strand = exons[0][0], exons[0][1]
        coords = sorted(exons, key=lambda x: x[2])
        if len(coords) == 1:
            result[tid] = f"{chrom}:{strand}:{coords[0][2]}-{coords[0][3]}:mono"
        else:
            introns = ','.join(
                f"{coords[i][3]}-{coords[i+1][2]}"
                for i in range(len(coords) - 1)
            )
            result[tid] = f"{chrom}:{strand}:{introns}"
    return result


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def read_sqanti_classification(path):
    return pd.read_csv(path, sep='\t', low_memory=False)


def build_fusion_set(df, keys):
    """Return a set of intron-chain keys for fusion transcripts."""
    id_col = 'isoform' if 'isoform' in df.columns else df.columns[0]
    fusion_ids = df.loc[df['structural_category'] == 'fusion', id_col]
    return {keys[t] for t in fusion_ids if t in keys}


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def save_venn2(set_a, set_b, labels, title, filename):
    """Save a 2-way Venn diagram."""
    fig, ax = plt.subplots(figsize=(6, 5))
    venn2([set_a, set_b], set_labels=labels, ax=ax)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sfari_classification',    required=True,
                        help='SFARI final classification parquet file')
    parser.add_argument('--sfari_gtf',               required=True,
                        help='SFARI transcripts GTF')
    parser.add_argument('--patowary_classification', required=True,
                        help='Patowary SQANTI classification TSV')
    parser.add_argument('--patowary_gtf',            required=True,
                        help='Patowary transcripts GTF (hg38)')
    args = parser.parse_args()

    # --- SFARI (parquet) ---
    sfari_df   = pq.read_table(args.sfari_classification).to_pandas()
    sfari_keys = get_intron_chain_keys(args.sfari_gtf)
    sfari_set  = build_fusion_set(sfari_df, sfari_keys)

    # --- Patowary ---
    pat_df   = read_sqanti_classification(args.patowary_classification)
    pat_keys = get_intron_chain_keys(args.patowary_gtf)
    pat_set  = build_fusion_set(pat_df, pat_keys)

    print(f"SFARI fusion transcripts:    {len(sfari_set)}")
    print(f"Patowary fusion transcripts: {len(pat_set)}")
    print(f"Overlap:                     {len(sfari_set & pat_set)}")

    save_venn2(
        sfari_set, pat_set,
        labels=('Xu and Rynard et al.', 'Patowary et al.'),
        title='Fusion transcripts (SFARI / Patowary)',
        filename='venn_fusion_sfari_patowary.pdf',
    )
    print("Saved: venn_fusion_sfari_patowary.pdf")


if __name__ == '__main__':
    main()
