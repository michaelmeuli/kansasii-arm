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

Then in Jalview: **File → Input Alignment** `$ARM_DIR/work/aln/<gene>.codon.aln.fasta`
(coding genes) or `$ARM_DIR/work/aln/<gene>.nt.aln.fasta` (non-coding genes such as
`rrl` and `rrs`, which have no codon alignment),
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

**Sample/gene QC.** A call is only trusted when the sample is a non-contaminated member of the
*M. kansasii* complex (`species` and `contamination_flag` from `output/screening_map_results.csv`,
override with `--sample-meta`) and the gene was extracted with >= `--min-ident` (default 80) % BLAST
identity (a floor against wrong-gene hits; the species check does the real work) and >= `--min-cov` (default 0.9, `gyrA` 0.6 because most samples lack the
reference's ~1.26 kb insertion; override with `--min-cov-gene GENE=FRAC`) coverage of the reference gene. Otherwise the `qc` column says why
(`contaminated`, `non_complex(<species>)`, `low_identity(<pct>)`, `partial(<cov>)`, `no_hit`), `differs` is
`NOCALL` and `aa_call` is empty. `marker_calls.pass.tsv` contains only the rows with `qc == ok`; use it for
analysis. Without the metadata file only the identity/coverage checks run.

### `marker_calls.tsv` columns

| Column | Meaning |
|---|---|
| `gene`, `drug` | Gene and the drug the marker is associated with |
| `marker` | Literature mutation (wild-type aa, position, mutant aa, e.g. `E92D`) or nucleotide marker |
| `evidence` | Support level (`VALIDATED`, `VALIDATED_OTHER_NTM`, `VALIDATED_OTHER_SPECIES`, `SUGGESTIVE`, `DO_NOT_CALL`, ...); sets the Jalview colour |
| `numbering_ref` | External reference the marker's numbering comes from (e.g. `MTB_GIDB_PROT`) |
| `ext_position` | Position in that external numbering |
| `mkn_position` | Same residue in *M. kansasii* (ATCC 12478) coordinates |
| `aln_column` | Column in the multi-sample MAFFT alignment |
| `ref_state` | Reference codon/base at that column |
| `sample` | `<gene>__<sample id>` |
| `sample_state` | Sample codon/base at that column |
| `differs` | `YES` / `SYN` / `no` / `NOCALL` (see above) |
| `ref_aa`, `sample_aa` | Translated reference / sample residue (coding markers) |
| `marker_wt_aa`, `marker_mut_aa` | Wild-type and mutant residue (or base, for single-nucleotide markers such as rRNA `A2058G`) the marker expects |
| `ref_aa_check` | Does the reference carry `wt`, `mut` or `neither` (`neither` = suspect transfer); also set for single-nucleotide markers. `mut` means the reference strain itself carries the mutant allele (e.g. rrs `A514C`), so `MUT` is then the species baseline, not resistance |
| `aa_call` | `wt` / `MUT` / `OTHER` relative to the marker's own alleles (amino acid, or base for single-nucleotide markers); `OTHER` = changed, but not to the marker's mutant (`differs = YES` alone does not mean the marker mutation); empty if QC failed or for region scans |
| `species`, `contamination_flag` | Sample metadata from `screening_map_results.csv` |
| `hit_pident`, `hit_cov` | BLAST identity (%) and coverage (fraction) of the extracted gene |
| `qc` | `ok` or the reason the row is masked (see QC above) |
| `mixed` | Read-mapping mode only: 10-90 % alt-allele sites inside the marker window as `genomic_pos:ref>alt@alt_fraction` (plus-strand genomic coordinates), `;`-separated. The consensus collapses these to the majority base, so a non-empty `mixed` means heteroresistance, a mixed culture or a second gene copy; empty in assembly mode |
| `doi` | Source paper for the marker |

Markers with `numbering_ref = NONE` cannot be placed automatically and are listed in
`unplaced_markers.tsv` (same directory) for manual curation.

## Read-mapping mode (genotyping from reads instead of assemblies)

The default run extracts genes from Unicycler assemblies. A mixed culture or a non-complex
species can give a misassembled or non-orthologous gene copy, which shows up as false
"mutations" (this is how apparent rpoB S531/H526 changes in 035/127/129 arose). Read mapping is
independent of the assembly, so the same 54 markers can also be genotyped from the trimmed reads:

```bash
./run_arm_bam.sh prepare                    # once: shared refs/ symlink, targets, bwa index
sbatch submit_arm_bam.sbatch                # array, one task per sample (06_bam_sample.sh)
sbatch --wrap '...; ./run_arm_bam.sh finish'  # after the array: 04_align + 05_annotate (needs env_arm)
./07_compare_modes.py                       # asm vs bam
```

- `06_bam_sample.sh ID`: `bwa mem` of `runs/mkan329/assembly/results/<ID>/0_trimming/<ID>_r{1,2}.fastq.gz`
  to `refs/reference.fna`, then `03_sample_seqs.sh bam` (depth < `MIN_DEPTH` masked to N, minority
  alleles in `mixed_sites.tsv`). The full BAM is discarded; `bam/<ID>.targets.bam` keeps only the
  target regions.
- Outputs go to `$KANSASII_ROOT/output/arm_bam` (`ARM_DIR`), so `output/arm` is untouched; `refs/`
  is shared via symlink. `marker_calls.tsv` has the same columns, except that `hit_pident` and
  `hit_cov` are empty (no BLAST step), so only the species/contamination QC applies.
- The consensus is built per target region. A whole-genome consensus shifts coordinates after
  every upstream indel, so slicing it with reference coordinates gives the wrong codons.
- `07_compare_modes.py` joins both `marker_calls.tsv` on gene x marker x sample and writes
  `mode_comparison.tsv` and `mode_disagreements.tsv` (classes `agree`, `disagree`, `asm_only`,
  `bam_only`, `both_nocall`; NOCALL is never a call). On the 131 Mkan329 samples: 3,550 shared
  calls, 0 disagreements. Reads fill gaps where the assembly extraction was partial (rrl, rrs,
  gyrA); the assembly gives rrl/rrs calls where read mapping is NOCALL (mapping to the rRNA genes
  is not deep or unique enough).
- Limits: reads are mapped to the ATCC 12478 reference, so divergent complex species map worse and
  non-complex species mostly do not map; the consensus carries SNPs and small indels only.

## Predicted resistance vs MIC

`08_mut_with_mic.py` keeps the `aa_call = MUT` (qc ok) rows of `marker_calls.pass.tsv` and joins
them by sample id to the parsed MHK (`output/mic/mhk/mic_parsed.csv`) and MGIT
(`output/mic/mgit/mgit_parsed.csv`) tables. Output `mut_with_mic.tsv` (next to `marker_calls.tsv`)
has one row per marker call x MIC row; `mic_matches_drug = 1` marks the MIC of the drug the
marker predicts. Samples without any MIC are dropped (listed in the printed summary). The MGIT
panel has no streptomycin, so for the streptomycin markers only MHK rows can match. Both
streptomycin markers are reference-polarity alleles (the ATCC 12478 reference carries the mutant
residue/base), so `MUT` means "same as the reference", which is why rrs A514C is `MUT` in every sample.

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

### `numbering_ref` values

`numbering_ref` names the external sequence a marker's `position` is counted in.
`01_fetch_reference.sh` downloads each one from NCBI into `$ARM_DIR/refs` as `<name>.fa`.
`MTB_GIDB_PROT`, for example, is *M. tuberculosis* H37Rv GidB (RsmG, Rv3919c,
`NP_218436.1`), so gid E92D means residue 92 of that protein.

| `numbering_ref` | Source | Used for |
|---|---|---|
| `MTB_RPOB_PROT`, `MTB_GYRB_PROT`, `MTB_RPSL_PROT`, `MTB_GIDB_PROT`, `MTB_RPLC_PROT`, `MTB_RPLD_PROT` | Mtb H37Rv proteins (RefSeq `NP_…`) | rpoB, gyrB, rpsL, gid, rplC, rplD |
| `MAV_GYRA_PROT` | *M. avium* GyrA (`WP_011723278.1`) | gyrA D95/A91 fluoroquinolone numbering |
| `ECOLI_16S_NT`, `ECOLI_23S_NT` | *E. coli* K-12 rRNA (`J01695.2` 1268..2809, `V00331.1`) | rrs, rrl |
| `NONE` | no transfer | region scans, negative findings, study-specific coordinates |

### Evidence levels

The `evidence` column in `markers.tsv` sets the Jalview colour (`COLOURS` in
`05_annotate.py`). `evidence_species` says where the evidence was observed and `doi`
gives the paper.

| Level | Meaning |
|---|---|
| `VALIDATED` | Confirmed in *M. kansasii* (e.g. rpoB codons 513, 516, 526, 531). Red. |
| `VALIDATED_OTHER_NTM` | Confirmed in other NTM, e.g. MAC or *M. abscessus* (rrl A2058G, rrs A1408G, gyrA D95G). Orange. |
| `VALIDATED_OTHER_SPECIES` | Confirmed outside NTM, mostly Mtb (gid loss of function, E92D, R20P). Amber. |
| `VALIDATED_N1` | Seen in *M. kansasii*, but in a single isolate (rrl A2266C). Orange. |
| `SUGGESTIVE`, `SUGGESTIVE_N1` | Associated but not functionally confirmed (gyrB in *M. simiae*, rrs A128G in one isolate). Grey. |
| `REGION` | A scan window, not a single marker: report any non-synonymous change (rpoB RRDR, gyrA QRDR, rrs helix 44, rrl domain V). |
| `DO_NOT_CALL` | Known species polymorphism, not resistance (gid A205A, MAC/MABC rrl and gyr variants). Green. |
| `NEGATIVE_FINDING` | Published evidence that a gene does not explain the phenotype (no gyrA mutations in 17 CIP-R *M. kansasii*). |
| `GAP` | No known mechanism for the drug. |
| `HYPERMUTATOR` | nucS loss of function: a risk flag, not a resistance call. |
| `SCREENED_NEGATIVE` | Screened with no mutations found (whiB7). |

The gid streptomycin markers are `VALIDATED_OTHER_SPECIES`: they are transferred from Mtb,
not *M. kansasii* evidence.

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
