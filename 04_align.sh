#!/usr/bin/env bash
# 04_align.sh - build one MAFFT alignment per gene: reference + every sample.
#
# Usage: ./04_align.sh [THREADS]
#
# For coding genes this also produces a codon-aware alignment (translate -> align
# protein -> back-translate), which is what you want before calling codon numbers:
# a plain nucleotide alignment can open gaps out of frame and shift every
# downstream codon by one or two.
set -euo pipefail

THREADS=${1:-4}
ALN=work/aln
mkdir -p "$ALN"

# Which genes are protein-coding, taken from targets.tsv feature type
mapfile -t CODING < <(awk -F'\t' '!/^#/ && $5 ~ /CDS/ {print $1}' targets.tsv)
is_coding () { for g in "${CODING[@]}"; do [[ "${1%%_c[0-9]*}" == "$g" ]] && return 0; done; return 1; }

for refgene in work/ref_genes/*.fa; do
  tag=$(basename "$refgene" .fa)
  in="$ALN/$tag.input.fa"
  cat "$refgene" > "$in"
  found=0
  for s in work/samples/*/genes/"$tag".fa; do
    [[ -e "$s" ]] || continue
    cat "$s" >> "$in"; found=$((found+1))
  done
  if (( found == 0 )); then
    echo "[skip] $tag: reference only, no sample sequences"; continue
  fi

  # --adjustdirectionaccurately is intentionally NOT used: 03 already
  # strand-corrects, and direction flipping would silently rename sequences.
  mafft --auto --thread "$THREADS" --preservecase "$in" > "$ALN/$tag.nt.aln.fasta" 2> "$ALN/$tag.mafft.log"

  if is_coding "$tag"; then
    if command -v macse >/dev/null 2>&1; then
      macse -prog alignSequences -seq "$in" \
        -out_NT "$ALN/$tag.codon.aln.fasta" -out_AA "$ALN/$tag.aa.aln.fasta" \
        >> "$ALN/$tag.mafft.log" 2>&1 || echo "[warn] macse failed for $tag, using nt alignment"
    else
      # translatorx-style fallback with MAFFT on the protein level
      python3 tools/backtranslate.py --nt "$in" --threads "$THREADS" \
        --out-nt "$ALN/$tag.codon.aln.fasta" --out-aa "$ALN/$tag.aa.aln.fasta"
    fi
  fi
  echo "[ok] $tag  ($((found)) samples + reference)"
done

echo
echo "Alignments in $ALN/"
echo "  *.nt.aln.fasta     plain nucleotide - open these in Jalview"
echo "  *.codon.aln.fasta  codon-aware nucleotide (coding genes) - use for codon calls"
echo "  *.aa.aln.fasta     protein alignment"
