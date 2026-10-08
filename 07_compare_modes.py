#!/usr/bin/env python3
"""07_compare_modes.py - compare assembly-based and read-mapping marker calls.

Joins marker_calls.tsv of both modes on (gene, marker, sample) and writes
  mode_comparison.tsv       one row per key with both modes' differs / aa_call / qc / state
  mode_disagreements.tsv    rows where the two modes disagree (classified, see below)
and prints a summary. A call counts as "called" when differs is YES / SYN / no; NOCALL is not
a call and never counts as agreement or disagreement. Classes:
  agree             both called, same differs
  disagree          both called, differs differ        <- investigate
  asm_only/bam_only only one mode produced a call
  both_nocall       neither called
Rows where either mode failed qc are tagged qc_either=1 (disagreements there are expected).
"""
from __future__ import annotations

import argparse
import csv
import os
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
CALLED = {"YES", "SYN", "no"}
Row = dict[str, str]


def load(path: Path) -> dict[tuple[str, str, str], Row]:
    with open(path, encoding="utf-8", newline="") as fh:
        return {(r["gene"], r["marker"], r["sample"]): r for r in csv.DictReader(fh, delimiter="\t")}


def get(r: Row | None, field: str) -> str:
    return r[field] if r else ""


def classify(a: Row | None, b: Row | None) -> str:
    ca = a is not None and a["differs"] in CALLED
    cb = b is not None and b["differs"] in CALLED
    if ca and cb:
        assert a is not None and b is not None
        return "agree" if a["differs"] == b["differs"] else "disagree"
    if ca:
        return "asm_only"
    if cb:
        return "bam_only"
    return "both_nocall"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--asm", default=str(ROOT / "output/arm/work/report/marker_calls.tsv"))
    ap.add_argument("--bam", default=str(ROOT / "output/arm_bam/work/report/marker_calls.tsv"))
    ap.add_argument("--outdir", default=str(ROOT / "output/arm_bam/work/report"))
    args = ap.parse_args()

    asm, bam = load(Path(args.asm)), load(Path(args.bam))
    cols = ["gene", "marker", "sample", "class", "qc_either", "asm_differs", "bam_differs",
            "asm_aa_call", "bam_aa_call", "asm_state", "bam_state", "asm_qc", "bam_qc", "species"]
    rows: list[list[str]] = []
    for key in sorted(asm.keys() | bam.keys()):
        a, b = asm.get(key), bam.get(key)
        qc_either = int(any(r is not None and r["qc"] != "ok" for r in (a, b)))
        rows.append([*key, classify(a, b), str(qc_either), get(a, "differs"), get(b, "differs"),
                     get(a, "aa_call"), get(b, "aa_call"), get(a, "sample_state"), get(b, "sample_state"),
                     get(a, "qc"), get(b, "qc"), get(a, "species") or get(b, "species")])

    out = Path(args.outdir)
    for name, sel in (("mode_comparison.tsv", rows),
                      ("mode_disagreements.tsv", [r for r in rows if r[3] in ("disagree", "asm_only", "bam_only")])):
        with open(out / name, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(cols)
            w.writerows(sel)

    n = Counter((r[3], r[4]) for r in rows)
    print(f"{len(rows)} gene x marker x sample keys ({len(asm)} asm, {len(bam)} bam)")
    for cls in ("agree", "disagree", "asm_only", "bam_only", "both_nocall"):
        print(f"  {cls:12s} qc_ok={n[(cls, '0')]:5d}  qc_failed_in_a_mode={n[(cls, '1')]:5d}")
    print(f"wrote {out}/mode_comparison.tsv and mode_disagreements.tsv")


if __name__ == "__main__":
    main()
