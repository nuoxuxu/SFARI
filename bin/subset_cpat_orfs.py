#!/usr/bin/env python3
"""
Subset a CPAT ORF FASTA file to sequences listed in the CPAT best-ORF TSV (column 'ID').
"""

import argparse
import sys


def load_best_ids(tsv_path):
    ids = set()
    with open(tsv_path) as fh:
        header = fh.readline().split('\t')
        id_col = header.index('ID')
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) > id_col:
                ids.add(parts[id_col])
    return ids


def subset_fasta(fasta_path, ids, output_path):
    written = 0
    skipped = 0
    keep = False
    with open(fasta_path) as fh, open(output_path, 'w') as out:
        for line in fh:
            if line.startswith('>'):
                name = line[1:].split()[0]
                keep = name in ids
                if keep:
                    written += 1
                else:
                    skipped += 1
            if keep:
                out.write(line)
    return written, skipped


def main():
    parser = argparse.ArgumentParser(description='Subset CPAT ORF FASTA to best-ORF predictions')
    parser.add_argument('--orfs_fasta', required=True, help='CPAT ORF sequences FASTA (e.g. SFARI.ORF_seqs.fa)')
    parser.add_argument('--best_tsv',   required=True, help='CPAT best-ORF TSV with an ID column (e.g. SFARI.ORF_prob.best.tsv)')
    parser.add_argument('--output',     required=True, help='Output subsetted FASTA file')
    args = parser.parse_args()

    print('Loading best-ORF IDs...', flush=True)
    ids = load_best_ids(args.best_tsv)
    print(f'  {len(ids)} IDs loaded', flush=True)

    print('Subsetting FASTA...', flush=True)
    written, skipped = subset_fasta(args.orfs_fasta, ids, args.output)
    print(f'  Written: {written}, Skipped: {skipped}')

    missing = ids - set()  # check coverage
    # report IDs not found
    found_ids = set()
    with open(args.output) as fh:
        for line in fh:
            if line.startswith('>'):
                found_ids.add(line[1:].split()[0])
    not_found = ids - found_ids
    if not_found:
        print(f'Warning: {len(not_found)} IDs from TSV not found in FASTA', file=sys.stderr)


if __name__ == '__main__':
    main()
