#!/usr/bin/env python3
"""M4 step 1: extract genome-wide significant SNPs for MR instrument clumping."""
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "out", "mr")
os.makedirs(OUT, exist_ok=True)

FILES = {
    "endo": os.path.join(DATA, "finngen_R13_N14_ENDOMETRIOSIS.gz"),
    "vte": os.path.join(DATA, "finngen_R13_I9_VTE.gz"),
}

for trait, path in FILES.items():
    keep = []
    n = 0
    for chunk in pd.read_csv(path, sep="\t", chunksize=2_000_000,
                             dtype={"#chrom": str, "rsids": str},
                             low_memory=False):
        n += len(chunk)
        sub = chunk[chunk["pval"] < 5e-8].copy()
        if not sub.empty:
            keep.append(sub)
    df = pd.concat(keep, ignore_index=True)
    df["SNP"] = df["rsids"].astype(str).str.split(",").str[0]
    df = df.rename(columns={"#chrom": "chr", "pval": "p", "sebeta": "se",
                            "af_alt": "af"})
    df = df[["SNP", "chr", "pos", "ref", "alt", "beta", "se", "p", "af", "rsids"]]
    df.to_csv(os.path.join(OUT, f"{trait}_sig.tsv"), sep="\t", index=False,
              float_format="%.10g")
    df[["SNP", "p"]].to_csv(os.path.join(OUT, f"{trait}_clump.p"), sep="\t",
                            index=False, header=False, float_format="%.10g")
    print(trait, "scanned", n, "significant", len(df),
          "unique SNP", df["SNP"].nunique())
