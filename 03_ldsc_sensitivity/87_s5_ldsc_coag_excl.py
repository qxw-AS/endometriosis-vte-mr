#!/usr/bin/env python3
"""S5: LDSC rg after excluding F5/F2/ABO +/- 500 kb or +/- 1 Mb.

Uses the existing munged FinnGen R13 LDSC sumstats and 1000G EUR LD scores.
SNP positions come from the 1000G PLINK .bim files used for LD scoring.
"""
import os
import re
import subprocess

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "ldsc")
BIM = os.path.join(ROOT, "data", "1000g_eur", "1000G_EUR_Phase3_plink")
os.makedirs(OUT, exist_ok=True)

GENES = {
    "F5": ("1", 169511950, 169586481),
    "F2": ("11", 46719181, 46719771),
    "ABO": ("9", 133250400, 133275201),
}
WINDOWS = {"500kb": 500_000, "1Mb": 1_000_000}

# Build one SNP -> (chr, pos) map from the per-chromosome bim files.
bims = []
for fn in sorted(os.listdir(BIM)):
    if fn.endswith(".bim"):
        bims.append(pd.read_csv(os.path.join(BIM, fn), sep="\t", header=None,
                                names=["chr", "SNP", "cm", "pos", "A1", "A2"],
                                dtype={"chr": str}))
pos = pd.concat(bims, ignore_index=True)[["SNP", "chr", "pos"]]
pos = pos.drop_duplicates("SNP")
print(f"bim SNPs: {len(pos):,}")

rows = []
for label, pad in WINDOWS.items():
    intervals = [(c, max(1, s - pad), e + pad) for c, s, e in GENES.values()]
    keep = pd.Series(True, index=pos.index)
    for c, s, e in intervals:
        keep &= ~((pos["chr"] == c) & (pos["pos"].between(s, e)))
    drop_snps = set(pos.loc[~keep, "SNP"])
    print(f"{label}: excluding {len(drop_snps):,} SNPs in "
          f"{(list(zip(GENES, intervals)))}")

    files = []
    for trait in ("endo", "vte"):
        d = pd.read_csv(os.path.join(ROOT, "data", f"{trait}_ldsc.sumstats.gz"),
                        sep="\t")
        before = len(d)
        d = d[~d["SNP"].isin(drop_snps)]
        path = os.path.join(OUT, f"{trait}_ldsc_excl_{label}.sumstats.gz")
        d.to_csv(path, sep="\t", index=False, compression="gzip")
        files.append(path)
        print(f"  {trait}: {before:,} -> {len(d):,} SNPs")

    out_prefix = os.path.join(OUT, f"rg_excl_{label}")
    cmd = [
        os.path.join(ROOT, "data", "ldsc"), "rg",
        "--rg", ",".join(files),
        "--ref-ld-chr", os.path.join(ROOT, "data", "ref_ld", "LDscore."),
        "--w-ld-chr", os.path.join(ROOT, "data", "ref_w",
                                   "1000G_Phase3_weights_hm3_no_MHC",
                                   "weights.hm3_noMHC."),
        "--out", out_prefix,
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    log = res.stdout + "\n" + res.stderr
    with open(out_prefix + ".log", "w") as fh:
        fh.write(log)
    print(f"\n===== {label} =====")
    print(log[-1800:])
    rows.append(dict(window=label, log=log))

# Parse the two runs into a compact table.
def grab(pat, s, default=None):
    m = re.search(pat, s)
    return m.group(1) if m else default

summary = []
for r in rows:
    s = r["log"]
    summary.append(dict(
        window=r["window"],
        rg=grab(r"Genetic Correlation: ([\-\d.]+)", s),
        rg_se=grab(r"Genetic Correlation: [\-\d.]+ \(([\d.]+)\)", s),
        p=grab(r"P: ([\d.eE\-]+)", s),
        # the last "Intercept:" is the genetic-covariance intercept (gcov_int);
        # the first two are the per-trait h2 intercepts
        gcov_int=(re.findall(r"Intercept: ([\-\d.]+) \(", s) or [None])[-1],
    ))
summary = pd.DataFrame(summary)
summary.to_csv(os.path.join(OUT, "s5_rg_exclusion_summary.tsv"), sep="\t", index=False)
print("\n=== S5 summary ===")
print(summary.to_string(index=False))
