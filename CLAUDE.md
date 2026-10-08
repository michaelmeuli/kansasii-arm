# CLAUDE.md

See README.md. Easy to get wrong:

- **Environment:** `env_arm` (`environment.yml`); system python has no samtools/mafft/blast.
  Type check with `~/data/conda/envs/env_typecheck/bin/mypy` (strict, py3.11).
- **Paths:** `env.sh` (shell) and the `ARM`/`ROOT` constants (python) derive
  `$KANSASII_ROOT/output/arm` (override `ARM_DIR`). refs/ and work/ are NOT in the repo.
  `targets.tsv` / `markers.tsv` are resolved relative to the script, not the cwd.
- **Real runs go through `sbatch submit_arm.sh`-style jobs** (`submit_arm.sbatch`), not the
  login node. Input assemblies: `runs/mkan329/assembly/results/*/1_unicycler/*.fasta`.
- **No hardcoded coordinates:** everything comes from the GFF in `02_make_targets.py`.
- **Marker numbering is external** (E. coli rRNA, Mtb/M. avium protein); `05_annotate.py`
  transfers it by alignment. Spot-check `mkn_position` after any reference change.
- **`NOCALL` is not wild-type** (depth-masked N or gap). Research use only, not clinical.
- `markers.tsv` / `targets.tsv` are tab-separated; keep them that way.
