#!/usr/bin/env python3
"""
Venn diagram of unique CDS sequence overlaps across ORFanage, GeneMark, and CPAT.
Sequences are compared by their nucleotide content (uppercased, gaps stripped).
"""

import argparse
from matplotlib_venn import venn3
import matplotlib.pyplot as plt


def load_sequences(fasta_path):
    """Return set of unique uppercased sequences from a FASTA file."""
    seqs = set()
    current = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current:
                    seqs.add(''.join(current).upper())
                current = []
            else:
                current.append(line)
    if current:
        seqs.add(''.join(current).upper())
    return seqs


def main():
    parser = argparse.ArgumentParser(description='Venn diagram of CDS sequence overlaps')
    parser.add_argument('--orfanage',  required=True, help='ORFanage CDS FASTA')
    parser.add_argument('--genemark',  required=True, help='GeneMarkS-T CDS FASTA')
    parser.add_argument('--cpat',      required=True, help='CPAT best-ORF FASTA')
    parser.add_argument('--output',    required=True, help='Output figure path (e.g. venn_orf_overlap.pdf)')
    args = parser.parse_args()

    print('Loading sequences...', flush=True)
    orfanage = load_sequences(args.orfanage)
    genemark = load_sequences(args.genemark)
    cpat     = load_sequences(args.cpat)
    print(f'  ORFanage : {len(orfanage):,} unique sequences')
    print(f'  GeneMark : {len(genemark):,} unique sequences')
    print(f'  CPAT     : {len(cpat):,} unique sequences')

    fig, ax = plt.subplots(figsize=(7, 6))
    venn3(
        [orfanage, genemark, cpat],
        set_labels=('ORFanage', 'GeneMarkS-T', 'CPAT'),
        ax=ax,
    )
    ax.set_title('Overlap of unique CDS sequences\nacross ORF prediction methods')
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f'Figure saved to {args.output}')


if __name__ == '__main__':
    main()
