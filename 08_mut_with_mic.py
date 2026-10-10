#!/usr/bin/env python3
"""08_mut_with_mic.py - predicted-resistant (aa_call = MUT) samples that also have a MIC.

Takes marker_calls.pass.tsv, keeps rows with aa_call == MUT (qc ok), and joins them on the sample id
(PROBENNUMMER) to the parsed phenotype tables:
  MHK   output/mic/mhk/mic_parsed.csv     Sensititre MIC (mg/L)
  MGIT  output/mic/mgit/mgit_parsed.csv   MGIT S/R per concentration (mg/L)
Writes mut_with_mic.tsv, one row per (marker call x MIC row). Samples with no MIC in either
source are dropped (listed in the summary). mic_matches_drug = 1 where the MIC antibiotic is the
drug the marker predicts (e.g. streptomycin); the MGIT panel has no streptomycin.
"""
from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("KANSASII_ROOT", "/shares/sander.imm.uzh/MM/kansasii"))
Row = dict[str, str]

COLS = ["gene", "marker", "drug", "sample", "species", "aa_call", "ref_aa_check", "discriminating", "evidence",
        "mic_source", "mic_drug", "mic_value_mg_l", "mic_int_erg", "mic_matches_drug", "doi"]


def read(path: Path, delimiter: str) -> list[Row]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def sample_id(r: Row) -> str:
    return r["sample"].split("__", 1)[1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calls", default=str(ROOT / "output/arm/work/report/marker_calls.pass.tsv"))
    ap.add_argument("--mhk", default=str(ROOT / "output/mic/mhk/mic_parsed.csv"))
    ap.add_argument("--mgit", default=str(ROOT / "output/mic/mgit/mgit_parsed.csv"))
    ap.add_argument("--out", default=str(ROOT / "output/arm/work/report/mut_with_mic.tsv"))
    args = ap.parse_args()

    mut = [r for r in read(Path(args.calls), "\t") if r["aa_call"] == "MUT" and r["qc"] == "ok"]

    mic: dict[str, list[tuple[str, str, str, str]]] = defaultdict(list)  # sample -> (source, drug, value, int_erg)
    for r in read(Path(args.mhk), ","):
        mic[r["PROBENNUMMER"]].append(("MHK", r["antibiotic"], r["mhk_parsed"], r["int_erg"]))
    for r in read(Path(args.mgit), ","):
        mic[r["PROBENNUMMER"]].append(("MGIT", r["antibiotic"], r["concentration_mg_l"], r["int_erg"]))

    rows: list[list[str]] = []
    for r in sorted(mut, key=lambda x: (x["gene"], x["marker"], sample_id(x))):
        s = sample_id(r)
        for src, drug, value, int_erg in mic.get(s, []):
            if not value:
                continue
            match = int(drug.lower() == r["drug"].lower())
            rows.append([r["gene"], r["marker"], r["drug"], s, r["species"], r["aa_call"], r["ref_aa_check"],
                         r["discriminating"], r["evidence"], src, drug, value, int_erg, str(match), r["doi"]])

    out = Path(args.out)
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(COLS)
        w.writerows(rows)

    print(f"{len(mut)} MUT marker calls; {len(rows)} rows written to {out}")
    for gene, marker in sorted({(r["gene"], r["marker"]) for r in mut}):
        called = {sample_id(r) for r in mut if r["gene"] == gene and r["marker"] == marker}
        kept = {r[3] for r in rows if r[0] == gene and r[1] == marker}
        str_ = {r[3] for r in rows if r[0] == gene and r[1] == marker and r[13] == "1"}
        mhk = {r[3] for r in rows if r[0] == gene and r[1] == marker and r[9] == "MHK"}
        mgit = {r[3] for r in rows if r[0] == gene and r[1] == marker and r[9] == "MGIT"}
        print(f"  {gene} {marker}: {len(called)} MUT samples; {len(kept)} with a MIC "
              f"(MHK {len(mhk)}, MGIT {len(mgit)}), {len(str_)} with a MIC for the predicted drug")
        dropped = sorted(called - kept)
        if len(dropped) <= 12:
            print(f"    no MIC: {', '.join(dropped) or '-'}")


if __name__ == "__main__":
    main()
