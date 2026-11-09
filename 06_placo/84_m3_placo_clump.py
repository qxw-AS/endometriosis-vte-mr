#!/usr/bin/env python3
"""M3 stage 4b: clump PLACO-significant SNPs into loci (PLINK 1.9, 1000G EUR).

PLACO (Ray 2020) clumps genome-wide significant variants with r2 < 0.1 within
1 Mb. PLINK resolves hg19 positions from the 1000G reference by rsID; PLACO
output (hg38) is mapped back through the hit table for region assignment.
"""
import os
import subprocess

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m3")
PLINK = os.path.join(ROOT, "data", "plink")
BFILE = os.path.join(ROOT, "data", "g1000", "1000G.EUR.QC.merged")

hits = pd.read_csv(os.path.join(OUT, "placo_hits_p1e5.tsv.gz"), sep="\t")
sig = hits[(hits["placo_p"] < 5e-8) & hits["rsid"].notna()].copy()
sig = sig.sort_values("placo_p").drop_duplicates("rsid")
sig[["rsid", "placo_p"]].rename(columns={"rsid": "SNP", "placo_p": "P"}).to_csv(
    os.path.join(OUT, "placo_clump.p"), sep="\t", index=False, header=True,
    float_format="%.6g")
print(f"PLACO genome-wide significant SNPs for clumping: {len(sig):,}")

subprocess.run([PLINK, "--bfile", BFILE, "--clump",
                os.path.join(OUT, "placo_clump.p"),
                "--clump-p1", "5e-8", "--clump-p2", "1",
                "--clump-r2", "0.1", "--clump-kb", "1000",
                "--out", os.path.join(OUT, "placo_clump"), "--allow-no-sex"],
               check=True, capture_output=True, text=True)

cl = pd.read_csv(os.path.join(OUT, "placo_clump.clumped"), sep=r"\s+")
cl = cl[["CHR", "SNP", "BP", "P", "TOTAL"]].copy()
cl["SNP"] = cl["SNP"].astype(str)

pos = hits.dropna(subset=["rsid"]).drop_duplicates("rsid").set_index("rsid")
cl["chr38"] = cl["SNP"].map(pos["chr"])
cl["pos38"] = cl["SNP"].map(pos["pos"])
cl["lead_placo_p"] = cl["SNP"].map(pos["placo_p"])
cl["beta_endo"] = cl["SNP"].map(pos["beta_e"])
cl["beta_vte"] = cl["SNP"].map(pos["beta_v"])

regions = pd.read_csv(os.path.join(OUT, "candidate_regions.tsv"), sep="\t",
                      dtype={"chr": str})


def assign(r):
    if pd.isna(r["chr38"]):
        return None, None
    c = str(int(r["chr38"]))
    g = regions[(regions["chr"] == c) & (regions["start"] <= r["pos38"]) &
                (regions["end"] >= r["pos38"])]
    if g.empty:
        return None, None
    return g.iloc[0]["region"], bool(g.iloc[0]["mhc"])


cl[["region", "mhc"]] = cl.apply(lambda r: pd.Series(assign(r)), axis=1)
cl = cl.sort_values("P")
cl.to_csv(os.path.join(OUT, "placo_loci.tsv"), sep="\t", index=False,
          float_format="%.6g")
print(f"{len(cl)} PLACO loci (r2<0.1, 1 Mb); "
      f"{cl['region'].notna().sum()} overlap a candidate region "
      f"({int(cl['mhc'].fillna(False).sum())} MHC)")
print(cl[["CHR", "SNP", "BP", "P", "TOTAL", "chr38", "pos38", "region",
          "beta_endo", "beta_vte"]].to_string(index=False,
                                              float_format=lambda v: f"{v:.3g}"))
