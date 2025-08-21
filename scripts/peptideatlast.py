#!/usr/bin/env python3
import sys, requests, time, re
from collections import Counter

API = "https://www.peptideatlas.org/api/promast/v1/map"

def read_fasta_peptides(path):
    pep = []
    with open(path) as fh:
        seq = []
        for line in fh:
            if line.startswith(">"):
                if seq:
                    pep.append("".join(seq))
                    seq = []
            else:
                seq.append(line.strip())
        if seq: pep.append("".join(seq))
    # keep only valid AA letters; drop empty
    pep = [re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', p.upper()) for p in pep]
    return [p for p in pep if p]

def pa_present(peptide, proteome="human"):
    # Returns True if peptide is found in the requested proteome
    r = requests.get(API, params={"peptide": peptide, "proteome": proteome}, timeout=20)
    r.raise_for_status()
    js = r.json()
    # ProMaST returns mappings when found; treat any non-empty mapping as a hit
    return bool(js.get("mappings"))

def main(fasta, proteome):
    peptides = list(dict.fromkeys(read_fasta_peptides(fasta)))  # unique, order-preserving
    hits = 0
    for i, p in enumerate(peptides, 1):
        try:
            if pa_present(p, proteome=proteome):
                hits += 1
        except Exception as e:
            sys.stderr.write(f"[warn] {p}: {e}\n")
        if i % 10 == 0:
            time.sleep(0.2)  # be polite
    pct = 100.0 * hits / len(peptides) if peptides else 0.0
    print(f"Proteome/build: {proteome}")
    print(f"Unique peptides checked: {len(peptides)}")
    print(f"In PeptideAtlas: {hits}")
    print(f"Percentage: {pct:.2f}%")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(f"Usage: {sys.argv[0]} peptides.fasta <proteome_name>  # e.g., human, mouse, plasma")
    main(sys.argv[1], sys.argv[2])
