#!/usr/bin/env bash
# Real-data download for Fig 2 (V6.1). All files public.
set -u
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D="${ENDO_VTE_DATA:-$ROOT/data}"
mkdir -p "$D"
cd "$D"

dl () {
  url="$1"; out="$2"
  echo "[$(date +%T)] START $out"
  curl -L --retry 5 --retry-delay 5 -C - --max-time 3600 -o "$out" "$url" \
    && echo "[$(date +%T)] DONE $out $(stat -c%s "$out")" \
    || echo "[$(date +%T)] FAIL $out"
}

dl "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_N14_ENDOMETRIOSIS.gz" finngen_R13_N14_ENDOMETRIOSIS.gz
dl "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_I9_VTE.gz" finngen_R13_I9_VTE.gz
dl "https://zenodo.org/records/10515792/files/1000G_Phase3_ldscores.tgz?download=1" 1000G_Phase3_ldscores.tgz
dl "https://zenodo.org/records/10515792/files/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1" 1000G_Phase3_weights_hm3_no_MHC.tgz
dl "https://zenodo.org/records/10515792/files/hm3_no_MHC.list.txt?download=1" hm3_no_MHC.list.txt
dl "https://zenodo.org/records/10515792/files/1000G_Phase3_plinkfiles.tgz?download=1" 1000G_Phase3_plinkfiles.tgz
echo "[$(date +%T)] ALL DOWNLOADS FINISHED"
