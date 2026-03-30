#!/usr/bin/env python3
"""
Make UpSet plots and Venn diagrams comparing transcripts across SFARI, Patowary,
Joglekar, and ENCODE4 studies.

Transcripts are compared by intron chain:
  multi-exon : 'chrom:strand:e1end-e2start,e2end-e3start,...'
  mono-exon  : 'chrom:strand:start-end:mono'

ISM filtering (Patowary and Joglekar only — SFARI is pre-filtered, ENCODE4 uses expression filter):
  Keep 5prime_fragment with within_polyA_site == True
  Keep 3prime_fragment with within_CAGE_peak  == True
  Drop all other ISM subcategories

ENCODE4 expression filter:
  Keep transcripts with >5 reads in >2 brodmann brain sample columns.
  Known/novel determined by transcript_novelty column ('Known' vs all others).
"""
import argparse
import re
from collections import defaultdict

import pandas as pd
import pyarrow.parquet as pq
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
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
    return pd.read_csv(path, sep='\t', low_memory=False)


def apply_ism_filter(df):
    """
    Keep ISM transcripts only if:
      5prime_fragment  with within_polyA_site == True
      3prime_fragment  with within_CAGE_peak  == True
    """
    is_ism  = df['structural_category'] == 'incomplete-splice_match'
    keep_5p = (df['subcategory'] == '5prime_fragment') & (df['within_polyA_site'] == True)
    keep_3p = (df['subcategory'] == '3prime_fragment') & (df['within_CAGE_peak']  == True)
    return df[~is_ism | keep_5p | keep_3p].copy()


def build_sets(df, keys):
    """Return (known_set, novel_set, all_set) using intron chain keys.
    'known' = full-splice_match; 'novel' = everything else.
    """
    id_col      = 'isoform' if 'isoform' in df.columns else df.columns[0]
    known_ids   = df.loc[df['structural_category'] == 'full-splice_match', id_col]
    novel_ids   = df.loc[df['structural_category'] != 'full-splice_match', id_col]
    known_set   = {keys[t] for t in known_ids if t in keys}
    novel_set   = {keys[t] for t in novel_ids if t in keys}
    return known_set, novel_set, known_set | novel_set


def read_encode4(tsv_path):
    """Read ENCODE4 abundance TSV and apply brain expression filter.
    Returns DataFrame filtered to transcripts with >5 reads in >2 brodmann columns.
    """
    df = pd.read_csv(tsv_path, sep='\t', low_memory=False)
    brodmann_cols = [c for c in df.columns if c.startswith('brodmann')]
    passing = (df[brodmann_cols] > 5).sum(axis=1) > 2
    return df[passing].copy()


def build_encode4_sets(df, keys):
    """Return (known_set, novel_set, all_set) for ENCODE4 using intron chain keys.
    'known' = transcript_novelty == 'Known'; 'novel' = all other novelty categories.
    """
    id_col    = 'annot_transcript_id'
    known_ids = df.loc[df['transcript_novelty'] == 'Known',  id_col]
    novel_ids = df.loc[df['transcript_novelty'] != 'Known',  id_col]
    known_set = {keys[t] for t in known_ids if t in keys}
    novel_set = {keys[t] for t in novel_ids if t in keys}
    return known_set, novel_set, known_set | novel_set


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def save_upset(named_sets, title, filename):
    """Save a 4-way UpSet plot using plain matplotlib."""
    labels = list(named_sets.keys())
    n      = len(labels)

    # Enumerate all non-empty intersections and their sizes
    intersections = {}
    for r in range(1, n + 1):
        for combo in itertools.combinations(range(n), r):
            members = set(named_sets[labels[combo[0]]])
            for i in combo[1:]:
                members &= named_sets[labels[i]]
            # subtract elements that belong to larger combos (exclusive intersection)
            for other in itertools.combinations(range(n), r + 1):
                if set(combo).issubset(other):
                    pass  # we want inclusive intersections for UpSet style
            intersections[combo] = len(members)

    # Keep only non-zero intersections, sort descending by size
    combos  = [c for c, s in intersections.items() if s > 0]
    sizes   = [intersections[c] for c in combos]
    order   = np.argsort(sizes)[::-1]
    combos  = [combos[i] for i in order]
    sizes   = [sizes[i]  for i in order]
    n_bars  = len(combos)

    fig = plt.figure(figsize=(max(10, n_bars * 0.6), 6))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax_bar = fig.add_subplot(gs[0])
    ax_dot = fig.add_subplot(gs[1], sharex=ax_bar)

    x = np.arange(n_bars)

    # Bar chart
    ax_bar.bar(x, sizes, color='steelblue', width=0.6)
    for xi, s in zip(x, sizes):
        ax_bar.text(float(xi), s + max(sizes) * 0.01, str(s),
                    ha='center', va='bottom', fontsize=7)
    ax_bar.set_ylabel('Intersection size')
    ax_bar.set_title(title)
    ax_bar.set_xlim(-0.5, n_bars - 0.5)
    ax_bar.tick_params(bottom=False, labelbottom=False)
    ax_bar.spines['bottom'].set_visible(False)

    # Dot matrix
    ax_dot.set_ylim(-0.5, n - 0.5)
    ax_dot.set_yticks(range(n))
    ax_dot.set_yticklabels(labels, fontsize=9)
    ax_dot.tick_params(bottom=False, labelbottom=False)
    ax_dot.spines['bottom'].set_visible(False)
    ax_dot.spines['top'].set_visible(False)

    for xi, combo in enumerate(combos):
        # grey background dots for all sets
        ax_dot.scatter([xi] * n, range(n), color='lightgrey', s=60, zorder=2)
        # filled dots for sets in this intersection
        ax_dot.scatter([xi] * len(combo), list(combo),
                       color='steelblue', s=60, zorder=3)
        # vertical line connecting filled dots
        if len(combo) > 1:
            ax_dot.plot([xi, xi], [min(combo), max(combo)],
                        color='steelblue', linewidth=2, zorder=2)

    fig.savefig(filename, bbox_inches='tight')
    plt.close(fig)


def save_venn(sets, labels, title, filename):
    """Save a 3-way Venn diagram."""
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
    parser.add_argument('--encode4_gtf',             required=True)
    parser.add_argument('--encode4_abundance',       required=True)
    args = parser.parse_args()

    # --- SFARI (parquet, already ISM-filtered) ---
    sfari_df              = pq.read_table(args.sfari_classification).to_pandas()
    sfari_keys            = get_intron_chain_keys(args.sfari_gtf)
    sfari_known, sfari_novel, sfari_all = build_sets(sfari_df, sfari_keys)

    # --- Patowary ---
    pat_df                = apply_ism_filter(read_sqanti_classification(args.patowary_classification))
    pat_keys              = get_intron_chain_keys(args.patowary_gtf)
    pat_known, pat_novel, pat_all = build_sets(pat_df, pat_keys)

    # --- Joglekar ---
    jog_df                = apply_ism_filter(read_sqanti_classification(args.joglekar_classification))
    jog_keys              = get_intron_chain_keys(args.joglekar_gtf)
    jog_known, jog_novel, jog_all = build_sets(jog_df, jog_keys)

    # --- ENCODE4 ---
    enc_df                = read_encode4(args.encode4_abundance)
    enc_keys              = get_intron_chain_keys(args.encode4_gtf)
    enc_known, enc_novel, enc_all = build_encode4_sets(enc_df, enc_keys)

    # --- 4-way UpSet plots ---
    for suffix, sets in [
        ('known', {'Xu and Rynard et al.': sfari_known, 'Patowary et al.': pat_known, 'Joglekar': jog_known, 'ENCODE4': enc_known}),
        ('novel', {'Xu and Rynard et al.': sfari_novel, 'Patowary et al.': pat_novel, 'Joglekar': jog_novel, 'ENCODE4': enc_novel}),
        ('all',   {'Xu and Rynard et al.': sfari_all,   'Patowary et al.': pat_all,   'Joglekar': jog_all,   'ENCODE4': enc_all}),
    ]:
        title = {'known': 'Known transcripts', 'novel': 'Novel transcripts', 'all': 'All transcripts'}[suffix]
        save_upset(sets, title, f'upset_{suffix}.pdf')

    # --- 3-way Venn: SFARI × ENCODE4 × Patowary ---
    venn_labels = ('Xu and Rynard et al.', 'ENCODE4', 'Patowary et al.')
    for suffix, triple, title in [
        ('known', (sfari_known, enc_known, pat_known), 'Known transcripts (SFARI / ENCODE4 / Patowary)'),
        ('novel', (sfari_novel, enc_novel, pat_novel), 'Novel transcripts (SFARI / ENCODE4 / Patowary)'),
        ('all',   (sfari_all,   enc_all,   pat_all),   'All transcripts (SFARI / ENCODE4 / Patowary)'),
    ]:
        save_venn(triple, venn_labels, title, f'venn_{suffix}_sfari_encode4_patowary.pdf')


if __name__ == '__main__':
    main()
