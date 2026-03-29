#!/usr/bin/env python3
"""
Make three-way Venn diagrams comparing transcripts across SFARI, Patowary, and Joglekar studies.

Transcripts are compared by intron chain:
  multi-exon : 'chrom:strand:e1end-e2start,e2end-e3start,...'
  mono-exon  : 'chrom:strand:start-end:mono'

This is GENCODE-version-agnostic — two transcripts match if and only if they
share the same internal splice junctions, regardless of TSS/TTS variation or
which reference version was used for alignment.

ISM filtering (Patowary and Joglekar only — SFARI is pre-filtered):
  Keep 5prime_fragment with within_polyA_site == True
  Keep 3prime_fragment with within_CAGE_peak  == True
  Drop all other ISM subcategories
"""
import argparse
import re
from collections import defaultdict

import pandas as pd
import pyarrow.parquet as pq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3


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
    """Read a SQANTI3 TSV classification file into a DataFrame."""
    return pd.read_csv(path, sep='\t', low_memory=False)


def apply_ism_filter(df):
    """
    For ISM transcripts keep only:
      5prime_fragment  with within_polyA_site == True
      3prime_fragment  with within_CAGE_peak  == True
    All other ISM rows are dropped.
    """
    is_ism = df['structural_category'] == 'incomplete-splice_match'
    keep_5p = (df['subcategory'] == '5prime_fragment') & (df['within_polyA_site'] == True)
    keep_3p = (df['subcategory'] == '3prime_fragment') & (df['within_CAGE_peak'] == True)
    return df[~is_ism | keep_5p | keep_3p].copy()


def build_sets(df, keys):
    """
    Return (fsm_set, non_fsm_set, all_set) for one study using intron chain keys.
    """
    id_col = 'isoform' if 'isoform' in df.columns else df.columns[0]
    fsm_ids     = df.loc[df['structural_category'] == 'full-splice_match',     id_col]
    non_fsm_ids = df.loc[df['structural_category'] != 'full-splice_match', id_col]

    fsm_set     = {keys[tid] for tid in fsm_ids     if tid in keys}
    non_fsm_set = {keys[tid] for tid in non_fsm_ids if tid in keys}

    return fsm_set, non_fsm_set, fsm_set | non_fsm_set


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def save_venn(sets, labels, title, filename):
    fig, ax = plt.subplots(figsize=(6, 5))
    venn3(sets, set_labels=labels, ax=ax)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(filename)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sfari_classification',    required=True)
    parser.add_argument('--sfari_gtf',               required=True)
    parser.add_argument('--patowary_classification', required=True)
    parser.add_argument('--patowary_gtf',            required=True)
    parser.add_argument('--joglekar_classification', required=True)
    parser.add_argument('--joglekar_gtf',            required=True)
    args = parser.parse_args()

    # --- SFARI (parquet, already ISM-filtered) ---
    sfari_df   = pq.read_table(args.sfari_classification).to_pandas()
    sfari_keys = get_intron_chain_keys(args.sfari_gtf)
    sfari_fsm, sfari_non_fsm, sfari_all = build_sets(sfari_df, sfari_keys)

    # --- Patowary ---
    pat_df   = apply_ism_filter(read_sqanti_classification(args.patowary_classification))
    pat_keys = get_intron_chain_keys(args.patowary_gtf)
    pat_fsm, pat_non_fsm, pat_all = build_sets(pat_df, pat_keys)

    # --- Joglekar ---
    jog_df   = apply_ism_filter(read_sqanti_classification(args.joglekar_classification))
    jog_keys = get_intron_chain_keys(args.joglekar_gtf)
    jog_fsm, jog_non_fsm, jog_all = build_sets(jog_df, jog_keys)

    labels = ('SFARI', 'Patowary', 'Joglekar')

    save_venn(
        (sfari_fsm, pat_fsm, jog_fsm),
        labels,
        'Full-splice match transcripts',
        'venn_fsm.pdf',
    )
    save_venn(
        (sfari_non_fsm, pat_non_fsm, jog_non_fsm),
        labels,
        'Non-FSM transcripts',
        'venn_non_fsm.pdf',
    )
    save_venn(
        (sfari_all, pat_all, jog_all),
        labels,
        'All transcripts',
        'venn_all.pdf',
    )


if __name__ == '__main__':
    main()
