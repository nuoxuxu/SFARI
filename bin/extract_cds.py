#!/usr/bin/env python3
"""
Parse GeneMark .lst file and extract protein-coding nucleotide sequences
from the corresponding FASTA file.
"""

import re
import sys
from pathlib import Path

def reverse_complement(seq):
    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    return seq.translate(comp)[::-1]

def parse_lst(lst_path):
    """Parse GeneMark .lst file, return dict of transcript -> list of (strand, start, end)."""
    content = Path(lst_path).read_text()
    blocks = re.split(r'Model information:.*?\n', content)[1:]

    transcript_genes = {}
    for block in blocks:
        name_match = re.search(r'FASTA definition line: (\S+)', block)
        if not name_match:
            continue
        name = name_match.group(1)
        # Coordinates are 1-based; strip leading < or > (partial gene markers)
        genes = re.findall(
            r'^\s+\d+\s+([+-])\s+<?>?(\d+)\s+<?>?(\d+)\s+\d+\s+\d+',
            block, re.MULTILINE
        )
        if genes:
            transcript_genes[name] = [(strand, int(left), int(right))
                                       for strand, left, right in genes]
    return transcript_genes

def load_fasta(fasta_path):
    """Load FASTA into dict {name: sequence}."""
    seqs = {}
    name = None
    parts = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs[name] = ''.join(parts)
    return seqs

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Extract CDS sequences using GeneMark predictions')
    parser.add_argument('--lst',    required=True, help='GeneMark .lst prediction file')
    parser.add_argument('--fasta',  required=True, help='Full transcript FASTA file')
    parser.add_argument('--output', required=True, help='Output CDS FASTA file')
    args = parser.parse_args()

    print('Parsing GeneMark predictions...', flush=True)
    transcript_genes = parse_lst(args.lst)
    print(f'  {len(transcript_genes)} transcripts with predictions', flush=True)

    print('Loading FASTA sequences...', flush=True)
    seqs = load_fasta(args.fasta)
    print(f'  {len(seqs)} sequences loaded', flush=True)

    missing = 0
    written = 0
    with open(args.output, 'w') as out:
        for transcript, genes in transcript_genes.items():
            if transcript not in seqs:
                missing += 1
                continue
            full_seq = seqs[transcript]
            for i, (strand, start, end) in enumerate(genes):
                # Convert to 0-based slice: start-1 .. end
                cds = full_seq[start - 1 : end]
                if strand == '-':
                    cds = reverse_complement(cds)
                suffix = f'.gene{i+1}' if len(genes) > 1 else ''
                header = f'>{transcript}{suffix} strand={strand} start={start} end={end} len={len(cds)}'
                out.write(header + '\n')
                # Write sequence in 60-char lines
                for j in range(0, len(cds), 60):
                    out.write(cds[j:j+60] + '\n')
                written += 1

    print(f'Written {written} CDS sequences to {args.output}')
    if missing:
        print(f'Warning: {missing} transcripts in .lst not found in FASTA', file=sys.stderr)

if __name__ == '__main__':
    main()
