#!/usr/bin/env bash
# 01_fetch_reference.sh - download GCF_000157895.3 (M. kansasii ATCC 12478, ASM15789v2)
# and the external sequences used for coordinate transfer.
#
# Assembly facts (verified): ASM15789v2, designated RefSeq reference for M. kansasii.
#   chromosome NZ_CP006835.1  6,432,277 bp
#   plasmid    NZ_CP006836.1    144,951 bp  (pMK12478)
#   total 6,577,228 bp, 5,817 genes in the RefSeq GFF
set -euo pipefail

ACC=GCF_000157895.3
ASM=ASM15789v2
BASE=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/157/895/${ACC}_${ASM}
REFDIR=${1:-refs}

mkdir -p "$REFDIR"
cd "$REFDIR"

for ext in genomic.fna.gz genomic.gff.gz protein.faa.gz cds_from_genomic.fna.gz; do
  f="${ACC}_${ASM}_${ext}"
  [[ -f "${f%.gz}" ]] && { echo "have ${f%.gz}"; continue; }
  echo "fetching $f"
  curl -fsSL -O "${BASE}/${f}"
  gunzip -f "$f"
done

ln -sf "${ACC}_${ASM}_genomic.fna" reference.fna
ln -sf "${ACC}_${ASM}_genomic.gff" reference.gff
samtools faidx reference.fna

echo
echo "Contigs in reference:"
cut -f1,2 reference.fna.fai

# --- external numbering references -------------------------------------------
# Published markers are expressed in E. coli rRNA numbering or in M. tuberculosis /
# M. avium protein numbering. 05_annotate.py transfers them onto M. kansasii by
# alignment. Fetch them here; adjust accessions if you prefer different references.
fetch_efetch () {  # $1=accession $2=db $3=outfile
  [[ -f "$3" ]] && { echo "have $3"; return; }
  echo "fetching $1 -> $3"
  curl -fsSL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$2&id=$1&rettype=fasta&retmode=text" -o "$3"
}

# E. coli K-12 MG1655 rRNA (rrnB operon) - the source of "E. coli numbering"
fetch_efetch "J01695.2"      nuccore ECOLI_16S_NT.fa   # contains 16S; trim to rrsB if needed
fetch_efetch "V00331.1"      nuccore ECOLI_23S_NT.fa   # 23S rRNA rrlB

# M. tuberculosis H37Rv proteins (Mtb codon numbering)
fetch_efetch "NP_215181.1"   protein MTB_RPOB_PROT.fa  # RpoB  Rv0667
fetch_efetch "NP_215566.1"   protein MTB_GYRB_PROT.fa  # GyrB  Rv0005
fetch_efetch "NP_216136.1"   protein MTB_RPSL_PROT.fa  # RpsL  Rv0682
fetch_efetch "NP_218389.1"   protein MTB_GIDB_PROT.fa  # GidB  Rv3919c
fetch_efetch "NP_215193.1"   protein MTB_RPLC_PROT.fa  # RplC  Rv0701
fetch_efetch "NP_215194.1"   protein MTB_RPLD_PROT.fa  # RplD  Rv0702

# M. avium GyrA - the numbering used by the D95/A91 fluoroquinolone literature
fetch_efetch "WP_003876997.1" protein MAV_GYRA_PROT.fa

echo
echo "NOTE: verify each accession above is the protein you expect before trusting a"
echo "      transferred codon number. Print headers with: head -1 refs/*_PROT.fa"
