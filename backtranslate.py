#!/usr/bin/env python3
"""Codon-aware alignment: translate -> MAFFT on protein -> back-translate to nt.

Keeps every gap a multiple of three so codon numbering stays meaningful.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile

CODONS = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(_B):
    for j, b2 in enumerate(_B):
        for k, b3 in enumerate(_B):
            CODONS[b1 + b2 + b3] = _AA[i * 16 + j * 4 + k]


def read_fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def translate(nt):
    nt = nt.upper().replace("-", "")
    aa = []
    for i in range(0, len(nt) - len(nt) % 3, 3):
        aa.append(CODONS.get(nt[i:i + 3], "X"))
    return "".join(aa)


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nt", required=True)
    ap.add_argument("--out-nt", required=True)
    ap.add_argument("--out-aa", required=True)
    ap.add_argument("--threads", default="4")
    a = ap.parse_args()

    nt = read_fasta(a.nt)
    aa = [(n, translate(s)) for n, s in nt]
    for n, s in aa:
        if "*" in s[:-1]:
            print(f"[warn] internal stop codon in {n} - frame or sequence problem", file=sys.stderr)

    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as tf:
        aa_in = tf.name
    write_fasta(aa_in, aa)

    with open(a.out_aa, "w") as out:
        subprocess.run(["mafft", "--auto", "--thread", a.threads, "--preservecase", aa_in],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)

    aligned_aa = dict(read_fasta(a.out_aa))
    ntd = {n: s.upper().replace("-", "") for n, s in nt}

    back = []
    for n, _ in nt:
        prot, src, i, outs = aligned_aa[n], ntd[n], 0, []
        for ch in prot:
            if ch == "-":
                outs.append("---")
            else:
                outs.append(src[i:i + 3].ljust(3, "N"))
                i += 3
        back.append((n, "".join(outs)))
    write_fasta(a.out_nt, back)


if __name__ == "__main__":
    main()
