#!/usr/bin/env python3
"""M3 stage 1b: stream both FinnGen R13 files and keep SNPs in candidate regions."""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "out", "m3")

FILES = {
    "endo": os.path.join(DATA, "finngen_R13_N14_ENDOMETRIOSIS.gz"),
    "vte": os.path.join(DATA, "finngen_R13_I9_VTE.gz"),
}
regions = pd.read_csv(os.path.join(OUT, "candidate_regions.tsv"), sep="\t",
                      dtype={"chr": str})
# interval lookup per chromosome
iv = {c: g[["start", "end"]].values for c, g in regions.groupby("chr")}

cols = ["#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"]
for trait, path in FILES.items():
    keep = []
    n = 0
    for chunk in pd.read_csv(path, sep="\t", chunksize=2_000_000,
                             usecols=cols, dtype={"#chrom": str, "rsids": str},
                             low_memory=False):
        n += len(chunk)
        mask = np.zeros(len(chunk), dtype=bool)
        for c, arr in iv.items():
            sub = (chunk["#chrom"].values == c)
            if not sub.any():
                continue
            pos = chunk["pos"].values
            inside = np.zeros(len(chunk), dtype=bool)
            for s, e in arr:
                inside |= (pos >= s) & (pos <= e)
            mask |= sub & inside
        if mask.any():
            keep.append(chunk[mask])
    df = pd.concat(keep, ignore_index=True)
    df = df[(df["pval"] > 0) & (df["sebeta"] > 0) &
            (df["af_alt"] > 0) & (df["af_alt"] < 1)].copy()
    df["rsids"] = df["rsids"].astype(str).str.split(",").str[0]
    df.to_csv(os.path.join(OUT, f"region_snps_{trait}.tsv.gz"), sep="\t",
              index=False, float_format="%.8g", compression="gzip")
    print(f"{trait}: scanned {n:,} variants, kept {len(df):,} in "
          f"{len(regions)} regions")
