# M. kansasii AMR target extraction → MAFFT → Jalview

Reference: **GCF_000157895.3 / ASM15789v2**, *Mycobacterium kansasii* ATCC 12478,
the designated RefSeq reference for the species.

| | |
|---|---|
| chromosome | `NC_022663.1` (= NZ_CP006835.1), 6,432,277 bp |
| plasmid | `NC_022654.1` (= NZ_CP006836.1; pMK12478), 144,951 bp |
| total | 6,577,228 bp |
| genes in RefSeq GFF | 5,817 |
| genome paper | Wang et al. 2015, [10.1093/gbe/evv035](https://doi.org/10.1093/gbe/evv035) |

## Why there are no hardcoded coordinates

Every start/end in this pipeline is resolved from the assembly's own GFF at
runtime by `02_make_targets.py`. Hardcoded coordinates silently rot across
assembly versions, and a one-off coordinate error in a resistance target is
invisible until it produces a wrong call. The resolver matches on `gene=`,
`Name=`, `locus_tag=` and a product regex, warns when a target has multiple
copies, and writes anything it cannot resolve to `work/targets.missing.txt`
with the grep command to investigate.

## Run

Paths: reference downloads and all outputs live outside the repo in
`$KANSASII_ROOT/output/arm/{refs,work}` (default root
`/shares/sander.imm.uzh/MM/kansasii`; override the whole dir with `ARM_DIR`).

**All Mkan329 assemblies** (inputs: `runs/mkan329/assembly/results/*/1_unicycler/*.fasta`):

```bash
conda env create -f environment.yml     # once: env_arm (samtools bcftools bedtools blast mafft)
sbatch submit_arm.sbatch                # 01 -> 02 -> 03 per new sample -> 04 -> 05
```

**Step by step** (deps: samtools bcftools bedtools blast mafft python3; macse optional):

```bash
./01_fetch_reference.sh                       # assembly + external numbering refs
python3 02_make_targets.py                    # resolve coordinates from the GFF

# per sample, either mode:
./03_sample_seqs.sh bam  PT001 mapped/PT001.bam
./03_sample_seqs.sh asm  PT002 assemblies/PT002.fasta

./04_align.sh 8
python3 05_annotate.py
```

Type check: `mypy` from the repo root (strict; see `pyproject.toml`).

Then in Jalview: **File → Input Alignment** `$ARM_DIR/work/aln/<gene>.codon.aln.fasta`,
then **File → Load Features** `$ARM_DIR/work/report/<gene>.jalview_features.txt`.
Features are coloured by evidence level (red = validated in *M. kansasii*,
orange = validated in other NTM, amber = other species, green = DO_NOT_CALL
species polymorphisms you should expect to see and ignore).

For scripted analysis use `$ARM_DIR/work/report/marker_calls.tsv`, one row per
marker × sample with `differs` = YES / SYN / no / NOCALL:
coding markers are compared at the amino-acid level (YES = residue differs from the
reference, SYN = synonymous nucleotide change only), other markers at the nucleotide level;
NOCALL = depth-masked N or alignment gap. Reference strain ATCC 12478 is *not* always the
wild type, so for substitution markers named like `K43R` read `aa_call` (MUT / wt / OTHER,
relative to the marker's own wild-type and mutant residue) and `ref_aa_check` (does the
reference carry wt, mut, or neither - `neither` means the position transfer is suspect).

## Coordinate transfer — the thing to understand

Almost nothing in `markers.tsv` is in *M. kansasii* coordinates:

- rrl / rrs positions (2058, 2059, 1408, 903–906, 514) are **E. coli numbering**
- rpoB, rpsL, gid, rplC, rplD, gyrB codons are **M. tuberculosis numbering**
- gyrA D95/A91 are **M. avium numbering** (= Mtb D94/A90)

`05_annotate.py` converts each one by pairwise-aligning the M. kansasii gene to
the external reference and walking the alignment. The `mkn_position` column in
`marker_calls.tsv` is the converted position — check a couple by hand the first
time you run this. The consensus numbering paper for mycobacterial rpoB is
[10.1016/j.cmi.2016.09.006](https://doi.org/10.1016/j.cmi.2016.09.006).

## Caveats that will cost you if ignored

1. **Codon-aware alignment matters.** A plain nucleotide MAFFT can open a 1 or
   2 bp gap and shift every downstream codon. `04_align.sh` builds
   `*.codon.aln.fasta` for CDS targets; `05_annotate.py` prefers it.
2. **rRNA copy number.** Slow-growing mycobacteria typically carry a single rrn
   operon, which is why acquired rrl/rrs point mutations are clinically
   meaningful here. If `02_make_targets.py` warns about multiple copies, stop
   and check — mixed operons change how you interpret a heterozygous call.
3. **Heteroresistance.** `03_sample_seqs.sh` in bam mode writes
   `mixed_sites.tsv` for any position between 10% and 90% alt fraction. The
   consensus collapses these; emerging resistance often appears as a minority
   allele first.
4. **Depth masking.** Positions below `MIN_DEPTH` (default 10) become N, which
   `05_annotate.py` reports as `NOCALL` rather than wild-type. Do not let a
   low-coverage N read as susceptible.
5. **Fluoroquinolones will mostly come back negative.** 17/85 ciprofloxacin-
   resistant *M. kansasii* isolates had no gyrA/gyrB mutation
   ([10.1128/AAC.01788-17](https://doi.org/10.1128/AAC.01788-17)). Efflux is the
   likely driver and is invisible to this pipeline.
6. **Doxycycline and SXT are not genotypable today.** No validated markers
   exist, despite 96% and 49% pooled phenotypic resistance.
7. This is a research tool. Nothing here is validated for clinical reporting;
   phenotypic DST remains the reference method.

## Streptomycin (newly added)

Not in the source meta-analysis but in CLSI's secondary panel for *M. kansasii*.
Targets added: `rpsL`, `gid` (gidB/rsmG), plus the streptomycin positions in
`rrs`. Note the rrs cluster II positions 1401/1402 sit inside the amikacin
helix-44 window — a variant there should be interpreted for both drugs. The
only *M. kansasii*-specific observation is rrs A128G in a single isolate with
streptomycin MIC >64 mg/L, which was never functionally confirmed.

## Windows checkout

- `.gitattributes` forces LF line endings, so a Windows checkout (even with
  `core.autocrlf=true`) keeps scripts runnable. Recommended: `git config core.autocrlf false`
  and `git config core.longpaths true`.
- Scripts default to the cluster root `/shares/sander.imm.uzh/MM/kansasii`; set the
  `KANSASII_ROOT` environment variable to point elsewhere (e.g. a mapped drive).
- The shell steps and sbatch only run on the cluster (or WSL).
