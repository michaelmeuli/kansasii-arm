#!/usr/bin/env python3
"""
05_annotate.py - map published resistance markers onto your MAFFT alignment and
write a Jalview features file plus a machine-readable call table.

The hard part this solves: almost no marker in markers.tsv is expressed in
M. kansasii coordinates. rrl/rrs positions are E. coli numbering; rpoB, gyrA,
rpsL, gid codons are M. tuberculosis or M. avium numbering. Writing those
numbers straight onto an M. kansasii alignment is wrong, sometimes by tens of
residues. So for each marker this script:

  1. aligns the M. kansasii reference gene to the external numbering reference
     (protein-level for CDS, nucleotide-level for rRNA),
  2. walks the pairwise alignment to convert external position -> M. kansasii
     position,
  3. converts M. kansasii position -> column in your multi-sample alignment,
  4. reads off what each sample actually has at that column.

Markers with numbering_ref = NONE are reported but not placed; they need manual
curation against the original paper's amplicon.

Usage:
  python3 05_annotate.py --aln-dir work/aln --refs refs --markers markers.tsv \
      --resolved work/targets.resolved.tsv --outdir work/report
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict

CODONS = {}
_B, _AA = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(_B):
    for j, b2 in enumerate(_B):
        for k, b3 in enumerate(_B):
            CODONS[b1 + b2 + b3] = _AA[i * 16 + j * 4 + k]

# Jalview feature colours, keyed by evidence level
COLOURS = {
    "VALIDATED": "e41a1c",
    "VALIDATED_OTHER_NTM": "ff7f00",
    "VALIDATED_OTHER_SPECIES": "ffbf00",
    "VALIDATED_N1": "ff7f00",
    "SUGGESTIVE": "999999",
    "REGION": "80b1d3",
    "DO_NOT_CALL": "4daf4a",
    "NEGATIVE_FINDING": "cccccc",
    "GAP": "cccccc",
    "HYPERMUTATOR": "984ea3",
    "SCREENED_NEGATIVE": "cccccc",
}


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
    return "".join(CODONS.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def mafft_pair(a_name, a_seq, b_name, b_seq):
    """Align two sequences, return the aligned pair."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        tf.write(f">{a_name}\n{a_seq}\n>{b_name}\n{b_seq}\n")
        p = tf.name
    res = subprocess.run(["mafft", "--auto", "--quiet", "--preservecase", p],
                         capture_output=True, text=True, check=True)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
        tf.write(res.stdout)
        q = tf.name
    recs = dict(read_fasta(q))
    os.unlink(p)
    os.unlink(q)
    return recs[a_name], recs[b_name]


def transfer(ext_seq, mkn_seq, ext_pos):
    """External 1-based position -> M. kansasii 1-based position, or None."""
    ea, ma = mafft_pair("EXT", ext_seq, "MKN", mkn_seq)
    ei = mi = 0
    for ec, mc in zip(ea, ma):
        if ec != "-":
            ei += 1
        if mc != "-":
            mi += 1
        if ei == ext_pos:
            return (mi if mc != "-" else None), (ec, mc)
    return None, (None, None)


def ungapped_to_column(aligned_seq):
    """Map 1-based ungapped position -> 0-based alignment column."""
    m, u = {}, 0
    for col, ch in enumerate(aligned_seq):
        if ch != "-":
            u += 1
            m[u] = col
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln-dir", default="work/aln")
    ap.add_argument("--refs", default="refs")
    ap.add_argument("--markers", default="markers.tsv")
    ap.add_argument("--resolved", default="work/targets.resolved.tsv")
    ap.add_argument("--outdir", default="work/report")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    numbering = {}
    with open(args.resolved) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            numbering[r["gene"]] = r["numbering_ref"]

    markers = defaultdict(list)
    with open(args.markers) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            markers[f[1]].append(dict(zip(
                ["drug", "gene", "kind", "marker", "position", "numbering_ref",
                 "evidence_species", "evidence", "notes", "doi"], f)))

    calls_path = os.path.join(args.outdir, "marker_calls.tsv")
    unplaced = []
    with open(calls_path, "w", newline="") as cf:
        cw = csv.writer(cf, delimiter="\t", lineterminator="\n")
        cw.writerow(["gene", "drug", "marker", "evidence", "numbering_ref",
                     "ext_position", "mkn_position", "aln_column",
                     "ref_state", "sample", "sample_state", "differs", "doi"])

        for aln_file in sorted(os.listdir(args.aln_dir)):
            if not aln_file.endswith((".codon.aln.fasta", ".nt.aln.fasta")):
                continue
            tag = aln_file.split(".")[0]
            gene = re.sub(r"_c\d+$", "", tag)
            if gene not in markers:
                continue
            # prefer the codon-aware alignment when both exist
            codon = os.path.join(args.aln_dir, f"{tag}.codon.aln.fasta")
            path = codon if os.path.exists(codon) else os.path.join(args.aln_dir, aln_file)
            if path != os.path.join(args.aln_dir, aln_file) and aln_file.endswith(".nt.aln.fasta"):
                continue

            recs = read_fasta(path)
            ref_name, ref_aln = recs[0]
            ref_ungapped = ref_aln.replace("-", "")
            colmap = ungapped_to_column(ref_aln)
            is_cds = os.path.exists(codon)

            feats = []
            for m in markers[gene]:
                nref, pos = m["numbering_ref"], m["position"]
                if nref == "NONE" or pos == "NA":
                    unplaced.append(m)
                    continue

                ext_path = os.path.join(args.refs, f"{nref}.fa")
                if not os.path.exists(ext_path):
                    unplaced.append({**m, "notes": m["notes"] + f" [missing {ext_path}]"})
                    continue
                ext_seq = read_fasta(ext_path)[0][1]
                mkn_cmp = translate(ref_ungapped) if nref.endswith("_PROT") else ref_ungapped

                lo, hi = (pos.split("-") + [pos])[:2] if "-" in pos else (pos, pos)
                cols, mkn_lo, mkn_hi = [], None, None
                for p in (int(lo), int(hi)):
                    mkn_p, (ec, mc) = transfer(ext_seq, mkn_cmp, p)
                    if mkn_p is None:
                        continue
                    if mkn_lo is None:
                        mkn_lo, ref_state_ext = mkn_p, (ec, mc)
                    mkn_hi = mkn_p
                if mkn_lo is None:
                    unplaced.append({**m, "notes": m["notes"] + " [no alignment anchor]"})
                    continue

                # protein position -> first nt of that codon
                def to_nt(p):
                    return (p - 1) * 3 + 1 if nref.endswith("_PROT") else p

                nt_lo, nt_hi = to_nt(mkn_lo), to_nt(mkn_hi) + (2 if nref.endswith("_PROT") else 0)
                c_lo, c_hi = colmap.get(nt_lo), colmap.get(nt_hi)
                if c_lo is None:
                    unplaced.append({**m, "notes": m["notes"] + " [outside alignment]"})
                    continue
                c_hi = c_hi if c_hi is not None else c_lo

                ref_state = ref_aln[c_lo:c_hi + 1]
                feats.append((m, c_lo, c_hi, ref_state, mkn_lo))

                for sname, saln in recs[1:]:
                    s_state = saln[c_lo:c_hi + 1]
                    differs = "YES" if s_state.upper() != ref_state.upper() and "N" not in s_state.upper() else \
                              ("NOCALL" if "N" in s_state.upper() or "-" in s_state else "no")
                    cw.writerow([tag, m["drug"], m["marker"], m["evidence"], nref,
                                 pos, mkn_lo, c_lo + 1, ref_state, sname, s_state,
                                 differs, m["doi"]])

            # ---- Jalview features file -------------------------------------
            jf = os.path.join(args.outdir, f"{tag}.jalview_features.txt")
            with open(jf, "w") as fh:
                used = {m["evidence"] for m, *_ in feats}
                for ev in sorted(used):
                    fh.write(f"{ev}\t{COLOURS.get(ev, '777777')}\n")
                fh.write("\nSTARTGROUP\tAMR\n")
                for m, c_lo, c_hi, ref_state, mkn_p in feats:
                    desc = f"{m['drug']} | {m['marker']} | {m['evidence']} | {m['evidence_species']} | doi:{m['doi']}"
                    # Jalview feature coords are 1-based on the named sequence,
                    # in UNGAPPED residue numbering - so emit per sequence.
                    for sname, saln in recs:
                        smap = {}
                        u = 0
                        for col, ch in enumerate(saln):
                            if ch != "-":
                                u += 1
                                smap[col] = u
                        s_lo, s_hi = smap.get(c_lo), smap.get(c_hi)
                        if s_lo is None:
                            continue
                        fh.write(f"{desc}\t{sname}\t-1\t{s_lo}\t{s_hi or s_lo}\t{m['evidence']}\n")
                fh.write("ENDGROUP\tAMR\n")
            print(f"[ok] {tag}: {len(feats)} markers placed -> {jf}")

    up = os.path.join(args.outdir, "unplaced_markers.tsv")
    with open(up, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["drug", "gene", "marker", "evidence", "reason_or_notes", "doi"])
        for m in unplaced:
            w.writerow([m["drug"], m["gene"], m["marker"], m["evidence"], m["notes"], m["doi"]])

    print(f"\ncalls:    {calls_path}")
    print(f"unplaced: {up}  <- review these by hand, they are not failures of the data")


if __name__ == "__main__":
    main()
