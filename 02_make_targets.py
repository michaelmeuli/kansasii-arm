#!/usr/bin/env python3
"""
02_make_targets.py - resolve target gene coordinates from the GCF_000157895.3 GFF.

Nothing is hardcoded: every coordinate below comes out of the assembly's own
annotation, so the output stays correct if you move to a different assembly
version or a different M. kansasii genome.

Outputs (into --outdir):
  targets.resolved.tsv  gene, contig, start, end, strand, locus_tag, product, copies
  targets.bed           exact gene bounds, 0-based half-open (for samtools/bedtools)
  targets.padded.bed    gene bounds + per-gene pad (promoters, 5' uORFs)
  targets.missing.txt   anything that did not resolve, with a hint
  ref_genes/<gene>.fa   reference gene sequence, already strand-corrected

Usage:
  python3 02_make_targets.py        # defaults: $KANSASII_ROOT/output/arm/{refs,work}
  python3 02_make_targets.py --gff G.gff --fasta G.fna --targets targets.tsv --outdir DIR
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any
from urllib.parse import unquote

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
ARM = Path(os.environ.get("ARM_DIR", ROOT / "output" / "arm"))
HERE = Path(__file__).resolve().parent

COMP = str.maketrans("ACGTacgtNnRYKMSWByrkmswbVvDdHh", "TGCAtgcaNnYRMKSWVrymkswvBbHhDd")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name: str | None = None
    buf: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def parse_attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for kv in field.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = unquote(v)
    return out


def load_targets(path: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "gene": f[0],
                    "drugs": f[1],
                    "aliases": {a.strip().lower() for a in f[2].split(",") if a.strip()},
                    "product_regex": None if f[3].strip().upper() == "NONE" else re.compile(f[3], re.I),
                    "ftypes": {t.strip() for t in f[4].split(",")},
                    "pad": int(f[5]),
                    "numbering_ref": f[6].strip(),
                }
            )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", default=str(ARM / "refs" / "reference.gff"))
    ap.add_argument("--fasta", default=str(ARM / "refs" / "reference.fna"))
    ap.add_argument("--targets", default=str(HERE / "targets.tsv"))
    ap.add_argument("--outdir", default=str(ARM / "work"))
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, "ref_genes"), exist_ok=True)
    targets = load_targets(args.targets)
    genome = read_fasta(args.fasta)

    # pass 1: collect every candidate feature per target
    hits: dict[str, list[dict[str, Any]]] = defaultdict(list)
    with open(args.gff, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            contig, _, ftype, start, end, _, strand, _, attr = f[:9]
            a = parse_attrs(attr)
            names = {a.get("gene", "").lower(), a.get("Name", "").lower(),
                     a.get("gene_synonym", "").lower(), a.get("locus_tag", "").lower()}
            names.discard("")
            product = a.get("product", "")

            for t in targets:
                if ftype not in t["ftypes"]:
                    continue
                by_name = bool(names & t["aliases"])
                by_prod = bool(t["product_regex"] and product and t["product_regex"].search(product))
                if by_name or by_prod:
                    hits[t["gene"]].append(
                        {
                            "contig": contig,
                            "start": int(start),       # GFF is 1-based inclusive
                            "end": int(end),
                            "strand": strand,
                            "locus_tag": a.get("locus_tag", a.get("ID", "")),
                            "product": product,
                            "matched_by": "name" if by_name else "product",
                        }
                    )

    resolved_path = os.path.join(args.outdir, "targets.resolved.tsv")
    bed_path = os.path.join(args.outdir, "targets.bed")
    pad_path = os.path.join(args.outdir, "targets.padded.bed")
    miss_path = os.path.join(args.outdir, "targets.missing.txt")

    with open(resolved_path, "w", newline="", encoding="utf-8") as rf, \
         open(bed_path, "w", encoding="utf-8") as bf, \
         open(pad_path, "w", encoding="utf-8") as pf, \
         open(miss_path, "w", encoding="utf-8") as mf:

        w = csv.writer(rf, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "copy", "contig", "start_1based", "end_1based", "strand",
                    "length_bp", "locus_tag", "matched_by", "numbering_ref", "drugs", "product"])

        for t in targets:
            gene = t["gene"]
            fs = hits.get(gene, [])
            # de-duplicate identical intervals (gene + CDS rows overlap)
            seen: set[tuple[str, int, int, str]] = set()
            uniq: list[dict[str, Any]] = []
            for h in fs:
                key = (h["contig"], h["start"], h["end"], h["strand"])
                if key not in seen:
                    seen.add(key)
                    uniq.append(h)
            uniq.sort(key=lambda h: (h["contig"], h["start"]))

            if not uniq:
                mf.write(f"{gene}\tNOT FOUND\tcheck aliases {sorted(t['aliases'])} "
                         f"or product regex against: grep -i '{gene}' {args.gff}\n")
                print(f"[MISS] {gene}", file=sys.stderr)
                continue

            if len(uniq) > 1:
                print(f"[WARN] {gene}: {len(uniq)} copies found - rRNA operons and paralogues "
                      f"are real; each is written as {gene}_c1, {gene}_c2, ...", file=sys.stderr)

            for i, h in enumerate(uniq, 1):
                tag = gene if len(uniq) == 1 else f"{gene}_c{i}"
                seq = genome[h["contig"]][h["start"] - 1:h["end"]]
                if h["strand"] == "-":
                    seq = revcomp(seq)

                w.writerow([tag, i, h["contig"], h["start"], h["end"], h["strand"],
                            h["end"] - h["start"] + 1, h["locus_tag"], h["matched_by"],
                            t["numbering_ref"], t["drugs"], h["product"]])

                bf.write(f"{h['contig']}\t{h['start']-1}\t{h['end']}\t{tag}\t0\t{h['strand']}\n")
                ps = max(0, h["start"] - 1 - t["pad"])
                pe = min(len(genome[h["contig"]]), h["end"] + t["pad"])
                pf.write(f"{h['contig']}\t{ps}\t{pe}\t{tag}_pad{t['pad']}\t0\t{h['strand']}\n")

                with open(os.path.join(args.outdir, "ref_genes", f"{tag}.fa"), "w", encoding="utf-8") as gf:
                    gf.write(f">REF_{tag}|{h['locus_tag']}|{h['contig']}:{h['start']}-{h['end']}({h['strand']})\n")
                    for j in range(0, len(seq), 70):
                        gf.write(seq[j:j + 70] + "\n")

    print(f"\nwrote {resolved_path}")
    print(f"wrote {bed_path} and {pad_path}")
    print(f"unresolved targets listed in {miss_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
