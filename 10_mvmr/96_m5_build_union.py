#!/usr/bin/env python3
"""M5 step 1: build the union of endometriosis, BMI and smoking instruments.

Instrument sources
- endometriosis: 36 harmonised FinnGen N14 lead SNPs (M4)
- BMI: FinnGen R13 BMI_IRN pheweb peaks (p<5e-8, peak=True)
- smoking: FinnGen R13 EVER_SMOKER pheweb peaks (p<5e-8, peak=True)

The union is LD-clumped against 1000G EUR (r2<0.001, 10 Mb) so that the
MVMR-IVW independence assumption is not violated by cross-phenotype LD.
"""
import gzip
import json
import os
import subprocess

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m5")
os.makedirs(OUT, exist_ok=True)

PEAK_FILES = {
    "BMI": "/tmp/man_BMI_IRN.json",
    "SMOKE": "/tmp/man_EVER_SMOKER.json",
    "PE": "/tmp/man_I9_PULMEMB.json",
    "DVT": "/tmp/man_DVT.json",
}


def load_peaks(name):
    j = json.loads(gzip.decompress(open(PEAK_FILES[name], "rb").read()))
    d = pd.DataFrame(j["unbinned_variants"])
    d = d[(d["pval"] < 5e-8) & (d["peak"] == True)].copy()
    d["SNP"] = d["rsids"].astype(str).str.split(",").str[0]
    d["chr"] = d["chrom"].astype(int)
    d["pos"] = d["pos"].astype(int)
    d["source"] = name
    return d[["SNP", "chr", "pos", "ref", "alt", "beta", "sebeta", "pval",
              "af_alt", "source"]]


endo = pd.read_csv(os.path.join(ROOT, "out", "mr", "harmonised_forward.tsv"),
                   sep="\t")
endo = endo[endo["use"]].copy()
endo["P"] = 2 * stats.norm.sf(np.abs(endo["b_exp"] / endo["se_exp"]))
endo["source"] = "ENDO"
endo = endo.rename(columns={"chr": "chr", "pos": "pos", "SNP": "SNP",
                            "ea": "alt", "nea": "ref", "P": "pval"})
endo["beta"] = endo["b_exp"]
endo["sebeta"] = endo["se_exp"]

parts = [load_peaks("BMI"), load_peaks("SMOKE")]
parts.append(endo[["SNP", "chr", "pos", "ref", "alt", "beta", "sebeta",
                   "pval", "source"]])
u = pd.concat(parts, ignore_index=True)
u["chr"] = u["chr"].astype(int)
u["pos"] = u["pos"].astype(int)
u = u.sort_values("pval").drop_duplicates("SNP", keep="first")
u["P"] = u["pval"].clip(lower=1e-300)
u[["SNP", "P"]].to_csv(os.path.join(OUT, "union_assoc.txt"), sep="\t",
                       index=False)
u.to_csv(os.path.join(OUT, "union_all.tsv"), sep="\t", index=False)

cmd = [os.path.join(ROOT, "data", "plink"),
       "--bfile", os.path.join(ROOT, "data", "g1000", "1000G.EUR.QC.merged"),
       "--clump", os.path.join(OUT, "union_assoc.txt"),
       "--clump-snp-field", "SNP", "--clump-field", "P",
       "--clump-p1", "5e-8", "--clump-r2", "0.001", "--clump-kb", "10000",
       "--out", os.path.join(OUT, "union_clump")]
print("running:", " ".join(cmd))
r = subprocess.run(cmd, capture_output=True, text=True)
print(r.stdout[-3000:])
print(r.stderr[-1000:])

cl = pd.read_csv(os.path.join(OUT, "union_clump.clumped"), sep=r"\s+")
print("clumped leads:", len(cl))
print(cl[["CHR", "SNP", "BP", "P"]].head(10).to_string(index=False))
print("source counts in clumped:")
print(cl.merge(u[["SNP", "source"]], on="SNP", how="left")["source"]
      .value_counts(dropna=False).to_string())
