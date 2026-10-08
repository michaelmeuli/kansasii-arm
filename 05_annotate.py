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

Sample/gene QC: a call is only trusted when the sample is a non-contaminated member of the
M. kansasii complex (species + contamination_flag from screening_map_results.csv) AND the
gene was extracted with >= --min-ident % BLAST identity and >= --min-cov coverage of the
reference gene. Otherwise the `qc` column says why, `differs` is NOCALL and aa_call is empty.
marker_calls.pass.tsv holds only the rows with qc == ok.

Usage:
  python3 05_annotate.py            # defaults: $KANSASII_ROOT/output/arm/{refs,work}
  python3 05_annotate.py --aln-dir DIR --refs DIR --markers markers.tsv \
      --resolved targets.resolved.tsv --outdir DIR
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
from pathlib import Path
from typing import Any

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
ARM = Path(os.environ.get("ARM_DIR", ROOT / "output" / "arm"))
HERE = Path(__file__).resolve().parent

COMPLEX = {"kansasii", "persicum", "pseudokansasii", "innocens", "attenuatum", "ostraviense", "gastri"}

CODONS: dict[str, str] = {}
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


def read_fasta(path: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    name: str | None = None
    buf: list[str] = []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        out.append((name, "".join(buf)))
    return out


def translate(nt: str) -> str:
    nt = nt.upper().replace("-", "")
    return "".join(CODONS.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def mafft_pair(a_name: str, a_seq: str, b_name: str, b_seq: str) -> tuple[str, str]:
    """Align two sequences, return the aligned pair."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False, encoding="utf-8") as tf:
        tf.write(f">{a_name}\n{a_seq}\n>{b_name}\n{b_seq}\n")
        p = tf.name
    res = subprocess.run(["mafft", "--auto", "--quiet", "--preservecase", p],
                         capture_output=True, text=True, check=True)
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False, encoding="utf-8") as tf:
        tf.write(res.stdout)
        q = tf.name
    recs = dict(read_fasta(q))
    os.unlink(p)
    os.unlink(q)
    return recs[a_name], recs[b_name]


def transfer(ext_seq: str, mkn_seq: str, ext_pos: int) -> tuple[int | None, tuple[str | None, str | None]]:
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


def ungapped_to_column(aligned_seq: str) -> dict[int, int]:
    """Map 1-based ungapped position -> 0-based alignment column."""
    m: dict[int, int] = {}
    u = 0
    for col, ch in enumerate(aligned_seq):
        if ch != "-":
            u += 1
            m[u] = col
    return m


def load_meta(path: str) -> dict[str, tuple[str, str]]:
    """sample id -> (species, contamination_flag) from screening_map_results.csv."""
    out: dict[str, tuple[str, str]] = {}
    if not os.path.exists(path):
        print(f"[warn] no sample metadata at {path}: species/contamination QC skipped", file=sys.stderr)
        return out
    with open(path, encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh):
            out[row["PROBENNUMMER"]] = (row.get("species", "").strip(), row.get("contamination_flag", "").strip())
    return out


def load_hit(samples_dir: str, sample: str, tag: str) -> tuple[float, float] | None:
    """(percent identity, fraction of the reference gene covered) of the extracted gene, or None."""
    path = os.path.join(samples_dir, sample, "blast", f"{tag}.hit")
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    with open(path, encoding="utf-8") as fh:
        f = fh.readline().split()
    return float(f[3]), float(f[4]) / float(f[5])


def sample_qc(sample: str, tag: str, meta: dict[str, tuple[str, str]], samples_dir: str,
              min_ident: float, min_cov: float) -> str:
    """'ok' or a ';'-joined list of reasons this sample/gene must not be called."""
    why: list[str] = []
    species, contam = meta.get(sample, ("", ""))
    if contam:
        why.append("contaminated")
    if species and species not in COMPLEX:
        why.append(f"non_complex({species})")
    hit = load_hit(samples_dir, sample, tag)
    if hit is None:
        if os.path.isdir(os.path.join(samples_dir, sample, "blast")):   # BLAST QC exists only in asm mode
            why.append("no_hit")
    else:
        pid, cov = hit
        if pid < min_ident:
            why.append(f"low_identity({pid:.1f})")
        if cov < min_cov:
            why.append(f"partial({cov:.0%})")
    return ";".join(why) or "ok"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln-dir", default=str(ARM / "work" / "aln"))
    ap.add_argument("--refs", default=str(ARM / "refs"))
    ap.add_argument("--markers", default=str(HERE / "markers.tsv"))
    ap.add_argument("--resolved", default=str(ARM / "work" / "targets.resolved.tsv"))
    ap.add_argument("--outdir", default=str(ARM / "work" / "report"))
    ap.add_argument("--samples", default=str(ARM / "work" / "samples"))
    ap.add_argument("--sample-meta", default=str(ROOT / "output" / "screening_map_results.csv"))
    ap.add_argument("--min-ident", type=float, default=80.0,
                    help="minimum BLAST %% identity of the extracted gene to the reference gene "
                         "(a floor against wrong-gene hits; species gating does the real work, and "
                         "complex species differ a lot at e.g. gid: 88%% persicum, 80%% attenuatum)")
    ap.add_argument("--min-cov", type=float, default=0.9,
                    help="minimum fraction of the reference gene covered by the extracted gene")
    ap.add_argument("--min-cov-gene", action="append", default=[], metavar="GENE=FRAC",
                    help="per-gene override of --min-cov (repeatable); built-in: gyrA=0.6, because "
                         "the reference gyrA carries a ~1.26 kb insertion that most samples lack")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    meta = load_meta(args.sample_meta)
    cov_gene = {"gyrA": 0.6}
    for item in args.min_cov_gene:
        g, _, v = item.partition("=")
        cov_gene[g] = float(v)
    qc_cache: dict[tuple[str, str], str] = {}
    numbering: dict[str, str] = {}
    with open(args.resolved, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            numbering[r["gene"]] = r["numbering_ref"]

    markers: defaultdict[str, list[dict[str, str]]] = defaultdict(list)
    with open(args.markers, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            markers[f[1]].append(dict(zip(
                ["drug", "gene", "kind", "marker", "position", "numbering_ref",
                 "evidence_species", "evidence", "notes", "doi"], f)))

    calls_path = os.path.join(args.outdir, "marker_calls.tsv")
    unplaced: list[dict[str, str]] = []
    pass_path = os.path.join(args.outdir, "marker_calls.pass.tsv")
    with open(calls_path, "w", newline="", encoding="utf-8") as cf, \
         open(pass_path, "w", newline="", encoding="utf-8") as pf:
        cw = csv.writer(cf, delimiter="\t", lineterminator="\n")
        pw = csv.writer(pf, delimiter="\t", lineterminator="\n")
        cw.writerow(["gene", "drug", "marker", "evidence", "numbering_ref",
                     "ext_position", "mkn_position", "aln_column",
                     "ref_state", "sample", "sample_state", "differs",
                     "ref_aa", "sample_aa", "marker_wt_aa", "marker_mut_aa",
                     "ref_aa_check", "aa_call", "species", "contamination_flag", "hit_pident",
                     "hit_cov", "qc", "doi"])
        pw.writerow(["gene", "drug", "marker", "evidence", "numbering_ref",
                     "ext_position", "mkn_position", "aln_column",
                     "ref_state", "sample", "sample_state", "differs",
                     "ref_aa", "sample_aa", "marker_wt_aa", "marker_mut_aa",
                     "ref_aa_check", "aa_call", "species", "contamination_flag", "hit_pident",
                     "hit_cov", "qc", "doi"])

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

            feats: list[tuple[dict[str, str], int, int, str, int]] = []
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
                mkn_lo: int | None = None
                mkn_hi: int | None = None
                for p in (int(lo), int(hi)):
                    mkn_p, (ec, mc) = transfer(ext_seq, mkn_cmp, p)
                    if mkn_p is None:
                        continue
                    if mkn_lo is None:
                        mkn_lo, ref_state_ext = mkn_p, (ec, mc)
                    mkn_hi = mkn_p
                if mkn_lo is None or mkn_hi is None:
                    unplaced.append({**m, "notes": m["notes"] + " [no alignment anchor]"})
                    continue

                # protein position -> first nt of that codon
                def to_nt(p: int) -> int:
                    return (p - 1) * 3 + 1 if nref.endswith("_PROT") else p

                nt_lo, nt_hi = to_nt(mkn_lo), to_nt(mkn_hi) + (2 if nref.endswith("_PROT") else 0)
                c_lo, c_hi = colmap.get(nt_lo), colmap.get(nt_hi)
                if c_lo is None:
                    unplaced.append({**m, "notes": m["notes"] + " [outside alignment]"})
                    continue
                c_hi = c_hi if c_hi is not None else c_lo

                ref_state = ref_aln[c_lo:c_hi + 1]
                feats.append((m, c_lo, c_hi, ref_state, mkn_lo))

                # Coding markers are compared at the amino-acid level so synonymous
                # changes (SYN) are not reported as differences.
                use_aa = nref.endswith("_PROT") and "-" not in ref_state and len(ref_state) % 3 == 0
                ref_aa = translate(ref_state) if use_aa else ""
                mm = re.match(r"^([A-Z])\d+([A-Z*])(?![A-Za-z0-9])", m["marker"])
                wt_aa, mut_aa = (mm.group(1), mm.group(2)) if mm and use_aa else ("", "")
                # does the reference itself carry the wild-type or the mutant residue?
                # "neither" means the position transfer is probably off.
                ref_check = ""
                if wt_aa and len(ref_aa) == 1:
                    ref_check = "wt" if ref_aa == wt_aa else ("mut" if ref_aa == mut_aa else "neither")

                for sname, saln in recs[1:]:
                    s_state = saln[c_lo:c_hi + 1]
                    # missing data (depth-masked N or alignment gap) must never read as
                    # a result: decide NOCALL before comparing to the reference
                    s_up = s_state.upper()
                    s_aa = ""
                    aa_call = ""
                    sid = sname.split("__", 1)[-1]
                    if (sid, tag) not in qc_cache:
                        qc_cache[(sid, tag)] = sample_qc(
                            sid, tag, meta, args.samples, args.min_ident,
                            cov_gene.get(re.sub(r"_c\d+$", "", tag), args.min_cov))
                    qc = qc_cache[(sid, tag)]
                    species, contam = meta.get(sid, ("", ""))
                    hit = load_hit(args.samples, sid, tag)
                    if qc != "ok":
                        # wrong species / contaminated / poorly extracted gene: whatever the
                        # alignment shows at this column is not a trustworthy genotype
                        differs = "NOCALL"
                    elif "N" in s_up or "-" in s_up:
                        differs = "NOCALL"
                    elif use_aa:
                        s_aa = translate(s_up)
                        if s_aa != ref_aa:
                            differs = "YES"
                        else:
                            differs = "SYN" if s_up != ref_state.upper() else "no"
                        if ref_check == "neither":
                            # the reference has neither the marker's wt nor mut residue, so
                            # the numbering does not fit this gene: do not make a call
                            aa_call = "UNVERIFIED"
                        elif wt_aa and wt_aa != mut_aa and len(s_aa) == 1:  # synonymous markers: no aa call
                            aa_call = "MUT" if s_aa == mut_aa else ("wt" if s_aa == wt_aa else "OTHER")
                    else:
                        differs = "YES" if s_up != ref_state.upper() else "no"
                    out_row = [tag, m["drug"], m["marker"], m["evidence"], nref,
                               pos, mkn_lo, c_lo + 1, ref_state, sname, s_state,
                               differs, ref_aa, s_aa, wt_aa, mut_aa, ref_check, aa_call,
                               species, contam, f"{hit[0]:.2f}" if hit else "",
                               f"{hit[1]:.3f}" if hit else "", qc, m["doi"]]
                    cw.writerow(out_row)
                    if qc == "ok":
                        pw.writerow(out_row)

            # ---- Jalview features file -------------------------------------
            jf = os.path.join(args.outdir, f"{tag}.jalview_features.txt")
            with open(jf, "w", encoding="utf-8") as fh:
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
    with open(up, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["drug", "gene", "marker", "evidence", "reason_or_notes", "doi"])
        for m in unplaced:
            w.writerow([m["drug"], m["gene"], m["marker"], m["evidence"], m["notes"], m["doi"]])

    print(f"\ncalls:    {calls_path}  (qc == ok only: {pass_path})")
    print(f"unplaced: {up}  <- review these by hand, they are not failures of the data")


if __name__ == "__main__":
    main()
