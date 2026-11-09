#!/usr/bin/env bash
# Convert FinnGen R13 sumstats to LDSC format on HapMap3 SNPs.
# A1 = ALT (effect allele), A2 = REF; Z = beta/sebeta; N = cases+controls.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ENDO_VTE_DATA:-$ROOT/data}"

make_ldsc () {
  f="$1"; n="$2"; out="$3"
  zcat "$f" | awk -v FS='\t' -v OFS='\t' -v N="$n" '
    NR==FNR { hm[$1]=1; next }
    FNR==1 { next }
    {
      split($5, a, ","); snp=a[1]
      if (!(snp in hm) || (snp in seen)) next
      a1=$4; a2=$3
      if (a1 !~ /^[ACGT]$/ || a2 !~ /^[ACGT]$/) next
      se=$10+0; if (se<=0) next
      p=$7+0; if (p<=0 || p>=1) next
      af=$11+0; if (af<0.01 || af>0.99) next
      z=$9/se
      if (z != z) next
      seen[snp]=1
      print snp, a1, a2, N, sprintf("%.6g", z), sprintf("%.6g", p)
    }' hm3_snps.txt - | gzip -c > "$out"
  echo "$out $(gzip -dc "$out" | wc -l) SNPs"
}

make_ldsc finngen_R13_N14_ENDOMETRIOSIS.gz 143156 endo_ldsc.sumstats.gz
make_ldsc finngen_R13_I9_VTE.gz 500186 vte_ldsc.sumstats.gz
