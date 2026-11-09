#!/usr/bin/env bash
# Build HapMap3 SNP table (GRCh38 positions + FinnGen p) for gene-level enrichment.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ENDO_VTE_DATA:-$ROOT/data}"

extract () {
  f="$1"; out="$2"
  zcat "$f" | awk -v FS='\t' -v OFS='\t' '
    NR==FNR { hm[$1]=1; next }
    FNR==1 { next }
    {
      split($5, a, ","); snp=a[1]
      if (!(snp in hm) || (snp in seen)) next
      ch=$1
      if (ch=="X") ch=23; else if (ch=="Y" || ch=="MT") next
      p=$7+0; if (p<=0 || p>=1) next
      seen[snp]=1
      print snp, ch, $2, p
    }' hm3_snps.txt - | gzip -c > "$out"
  echo "$out $(zcat "$out" | wc -l)"
}

extract finngen_R13_N14_ENDOMETRIOSIS.gz snp_endo.tsv.gz
extract finngen_R13_I9_VTE.gz snp_vte.tsv.gz
