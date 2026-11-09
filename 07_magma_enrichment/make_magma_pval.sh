#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ENDO_VTE_DATA:-$ROOT/data}"
make_p () {
  f="$1"; out="$2"
  { echo -e "SNP\tP"; zcat "$f" | awk -v FS='\t' -v OFS='\t' '
    NR==FNR { b[$1]=1; next }
    FNR==1 { next }
    {
      split($5,a,","); s=a[1]
      if (!(s in b) || (s in seen)) next
      p=$7+0; if (p<=0 || p>=1) next
      seen[s]=1; print s, p
    }' bim_snps.txt -; } | gzip -c > "$out"
  echo "$out $(zcat "$out" | wc -l)"
}
make_p finngen_R13_N14_ENDOMETRIOSIS.gz endo_magma.pval.gz
make_p finngen_R13_I9_VTE.gz vte_magma.pval.gz
