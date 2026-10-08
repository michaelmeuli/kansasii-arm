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
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

MODE=${1:?mode: bam|asm}
SAMPLE=${2:?sample id}
INPUT=${3:?bam or contigs fasta}

REF=$ARM/refs/reference.fna
BED=$ARM/work/targets.bed
OUT=$ARM/work/samples/$SAMPLE
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

  # Build the consensus per target region: a whole-genome consensus shifts coordinates after
  # every upstream indel, so slicing it with the reference BED would cut the wrong bases.
  : > "$OUT/genes.raw.fa"
  while IFS=$'\t' read -r chr s e name _ strand; do
    reg="$chr:$((s+1))-$e"
    samtools faidx "$REF" "$reg" \
      | bcftools consensus -m "$OUT/lowcov.bed" -M N -H A "$OUT/calls.vcf.gz" > "$OUT/frag.fa" 2>/dev/null
    samtools faidx "$OUT/frag.fa"
    frag=$(cut -f1 "$OUT/frag.fa.fai")
    if [[ "$strand" == "-" ]]; then samtools faidx -i "$OUT/frag.fa" "$frag"; else samtools faidx "$OUT/frag.fa" "$frag"; fi \
      | sed "1s|.*|>${name}::${reg}|" >> "$OUT/genes.raw.fa"
  done < "$BED"
  rm -f "$OUT/frag.fa" "$OUT/frag.fa.fai"
  ;;

asm)
  mkdir -p "$OUT/blast"
  makeblastdb -in "$INPUT" -dbtype nucl -out "$OUT/blast/db" >/dev/null
  : > "$OUT/genes.raw.fa"
  for g in "$ARM"/work/ref_genes/*.fa; do
    tag=$(basename "$g" .fa)
    blastn -query "$g" -db "$OUT/blast/db" -max_target_seqs 5 \
      -outfmt '6 sseqid sstart send pident length qlen bitscore' \
      | sort -k7,7gr \
      | awk -F'\t' -v OFS='\t' '
          # Merge split HSPs (e.g. an intein present in the reference but not in the
          # sample splits one gene into two hits): start from the best hit, then add
          # same-contig, same-strand hits while the merged span stays <= 1.5 x query length.
          NR==1 { sid=$1; plus=($2<=$3); lo=($2<$3?$2:$3); hi=($2<$3?$3:$2)
                  pid=$4; alen=$5; q=$6; bits=$7; next }
          $1==sid && (($2<=$3)==plus) {
            l=($2<$3?$2:$3); h=($2<$3?$3:$2)
            nlo=(l<lo?l:lo); nhi=(h>hi?h:hi)
            if (nhi-nlo+1 <= 1.5*q) { lo=nlo; hi=nhi; alen+=$5; bits+=$7 }
          }
          END { if (NR) { if (alen>q) alen=q
                  print sid, (plus?lo:hi), (plus?hi:lo), pid, alen, q, bits } }' \
      > "$OUT/blast/$tag.hit"
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
    with open(os.path.join(out, "genes", f"{tag}.fa"), "w", encoding="utf-8") as fh:
        fh.write(f">{tag}__{sample}\n")
        s = "".join(buf)
        for i in range(0, len(s), 70): fh.write(s[i:i+70] + "\n")
for line in open(src, encoding="utf-8"):
    if line.startswith(">"):
        flush(); name, buf = line.strip(), []
    else: buf.append(line.strip())
flush()
print(f"per-gene sequences for {sample} in {out}/genes/")
PY
