#!/usr/bin/env bash
# 03_sample_seqs.sh - pull per-gene sequences out of your NGS data.
#
# Two modes:
#   bam  : reads already mapped to refs/reference.fna -> masked haploid consensus
#   asm  : de novo contigs -> best blastn hit per gene, strand-corrected
#
# Usage:
#   ./03_sample_seqs.sh bam  SAMPLE_ID  aligned.bam
#   ./03_sample_seqs.sh asm  SAMPLE_ID  contigs.fasta
#
# Requires: samtools, bcftools (bam mode); blast+ (asm mode); seqkit optional.
set -euo pipefail

MODE=${1:?mode: bam|asm}
SAMPLE=${2:?sample id}
INPUT=${3:?bam or contigs fasta}

REF=refs/reference.fna
BED=work/targets.bed
OUT=work/samples/$SAMPLE
MIN_DEPTH=${MIN_DEPTH:-10}          # positions below this are masked to N
MIN_IDENT=${MIN_IDENT:-85}          # asm mode: minimum blast identity
mkdir -p "$OUT"

if [[ ! -s "$BED" ]]; then echo "run 02_make_targets.py first" >&2; exit 1; fi

case "$MODE" in
bam)
  samtools index -c "$INPUT" 2>/dev/null || samtools index "$INPUT"

  # Mask anything with insufficient depth so you never align a low-coverage
  # artefact and call it a resistance mutation.
  samtools depth -a -b "$BED" "$INPUT" \
    | awk -v d="$MIN_DEPTH" '$3 < d {print $1"\t"$2-1"\t"$2}' \
    | bedtools merge -i - > "$OUT/lowcov.bed" || : > "$OUT/lowcov.bed"

  bcftools mpileup -f "$REF" -R "$BED" -a AD,DP -d 2000 -Q 20 -q 20 -Ou "$INPUT" \
    | bcftools call -m --ploidy 1 -Oz -o "$OUT/calls.vcf.gz"
  bcftools index "$OUT/calls.vcf.gz"

  # Heteroresistance matters clinically - keep the mixed sites visible.
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%AD]\n' "$OUT/calls.vcf.gz" \
    | awk -F'\t' '{split($5,a,","); tot=a[1]+a[2]; if (tot>0 && a[2]/tot>0.10 && a[2]/tot<0.90) print $0"\tMIXED"}' \
    > "$OUT/mixed_sites.tsv" || : > "$OUT/mixed_sites.tsv"

  bcftools consensus -f "$REF" -m "$OUT/lowcov.bed" -M N \
    -H A "$OUT/calls.vcf.gz" > "$OUT/consensus.fna"
  samtools faidx "$OUT/consensus.fna"

  # Slice the gene regions out of the consensus using the same BED.
  bedtools getfasta -fi "$OUT/consensus.fna" -bed "$BED" -s -name \
    > "$OUT/genes.raw.fa"
  ;;

asm)
  mkdir -p "$OUT/blast"
  makeblastdb -in "$INPUT" -dbtype nucl -out "$OUT/blast/db" >/dev/null
  : > "$OUT/genes.raw.fa"
  for g in work/ref_genes/*.fa; do
    tag=$(basename "$g" .fa)
    blastn -query "$g" -db "$OUT/blast/db" -max_target_seqs 5 \
      -outfmt '6 sseqid sstart send pident length qlen bitscore' \
      | sort -k7,7gr | head -1 > "$OUT/blast/$tag.hit"
    if [[ ! -s "$OUT/blast/$tag.hit" ]]; then
      echo "[MISS] $tag: no blast hit in $SAMPLE" >&2; continue
    fi
    read -r sid ss se pid alen qlen bits < "$OUT/blast/$tag.hit"
    awk -v p="$pid" -v a="$alen" -v q="$qlen" -v m="$MIN_IDENT" -v t="$tag" \
      'BEGIN{ if (p<m) print "[WARN] "t": identity "p"% below "m"%" > "/dev/stderr";
              if (a < 0.9*q) print "[WARN] "t": hit covers only "a"/"q" bp" > "/dev/stderr" }'
    if (( ss <= se )); then
      samtools faidx "$INPUT" "$sid:$ss-$se"
    else
      samtools faidx -i "$INPUT" "$sid:$se-$ss"   # -i = reverse complement
    fi | sed "1s|.*|>${tag}|" >> "$OUT/genes.raw.fa"
  done
  ;;
*)
  echo "mode must be bam or asm" >&2; exit 1 ;;
esac

# Normalise headers to  <gene>__<sample>  so the alignment step can group them.
python3 - "$OUT/genes.raw.fa" "$SAMPLE" "$OUT" <<'PY'
import os, re, sys
src, sample, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(os.path.join(out, "genes"), exist_ok=True)
name, buf = None, []
def flush():
    if not name: return
    tag = re.split(r'[:\s(]', name)[0].lstrip('>')
    with open(os.path.join(out, "genes", f"{tag}.fa"), "w") as fh:
        fh.write(f">{tag}__{sample}\n")
        s = "".join(buf)
        for i in range(0, len(s), 70): fh.write(s[i:i+70] + "\n")
for line in open(src):
    if line.startswith(">"):
        flush(); name, buf = line.strip(), []
    else: buf.append(line.strip())
flush()
print(f"per-gene sequences for {sample} in {out}/genes/")
PY
