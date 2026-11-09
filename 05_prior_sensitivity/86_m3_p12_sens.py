#!/usr/bin/env python3
"""M3 p12 prior sensitivity at fixed sdY = 0.15.

Re-runs the single-causal-variant coloc.abf model over all 58 candidate
regions with p12 = 1e-5 (primary), 5e-6 and 1e-6.  p1 = p2 = 1e-4.
Outputs all five posteriors plus region-level pass counts.
"""
import os

import numpy as np
import pandas as pd
from scipy.special import logsumexp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m3")

P1 = P2 = 1e-4
P12_GRID = [1e-5, 5e-6, 1e-6]
SDY = 0.15
MIN_MAF = 0.01

regions = pd.read_csv(os.path.join(OUT, "candidate_regions.tsv"), sep="\t",
                      dtype={"chr": str})


def labf(beta, se, sdy=SDY):
    v = se ** 2
    w = sdy ** 2
    z = beta / se
    return 0.5 * (np.log(v) - np.log(v + w) + z ** 2 * w / (v + w))


def posteriors(l1, l2, p12):
    lh = np.array([
        0.0,
        np.log(P1) + logsumexp(l1),
        np.log(P2) + logsumexp(l2),
        np.log(P1) + np.log(P2) + logsumexp(l1) + logsumexp(l2),
        np.log(p12) + logsumexp(l1 + l2),
    ])
    return np.exp(lh - logsumexp(lh))


def load(trait):
    d = pd.read_csv(os.path.join(OUT, f"region_snps_{trait}.tsv.gz"), sep="\t",
                    dtype={"#chrom": str, "rsids": str},
                    usecols=["#chrom", "pos", "ref", "alt", "rsids", "pval",
                             "beta", "sebeta", "af_alt"])
    d = d.rename(columns={"#chrom": "chr", "pval": "p", "sebeta": "se"})
    d["maf"] = np.minimum(d["af_alt"], 1 - d["af_alt"])
    return d.drop_duplicates(["chr", "pos", "ref", "alt"])


endo = load("endo")
vte = load("vte")
m = endo.merge(vte, on=["chr", "pos", "ref", "alt"], suffixes=("_e", "_v"))
print(f"shared SNP keys: {len(m):,}")

rows = []
for _, r in regions.iterrows():
    g = m[(m["chr"] == r["chr"]) & (m["pos"] >= r["start"]) &
          (m["pos"] <= r["end"])]
    g = g[(g["maf_e"] >= MIN_MAF) & (g["maf_v"] >= MIN_MAF)]
    out = dict(region=r["region"], chr=r["chr"], n_snp=len(g), mhc=bool(r["mhc"]))
    if len(g) < 2:
        for p12 in P12_GRID:
            out.update({f"p12_{p12:g}_PP_H{i}": np.nan for i in range(5)})
        rows.append(out)
        continue
    l1 = labf(g["beta_e"].values, g["se_e"].values)
    l2 = labf(g["beta_v"].values, g["se_v"].values)
    for p12 in P12_GRID:
        pp = posteriors(l1, l2, p12)
        out.update({f"p12_{p12:g}_PP_H{i}": pp[i] for i in range(5)})
    rows.append(out)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "coloc_p12_sensitivity.tsv"), sep="\t",
           index=False, float_format="%.6g")

# Region-to-gene mapping for the three loci highlighted in the paper.
coag = pd.read_csv(os.path.join(OUT, "coag_loci_coloc.tsv"), sep="\t")
print("\ncoag columns:", coag.columns.tolist())
print(coag[["gene", "region", "PP.H3", "PP.H4"]].to_string(index=False))

print("\n=== p12 sensitivity: F5 / PROCR / ABO ===")
for gene, region in coag.set_index("gene")["region"].items():
    hit = res[res.region == region]
    if hit.empty:
        print(f"\n{gene} ({region}): not in candidate regions")
        continue
    row = hit.iloc[0]
    print(f"\n{gene} ({region})")
    for p12 in P12_GRID:
        print(f"  p12={p12:g}: H3={row[f'p12_{p12:g}_PP_H3']:.3f} "
              f"H4={row[f'p12_{p12:g}_PP_H4']:.3f}")

print("\n=== counts passing thresholds ===")
for p12 in P12_GRID:
    h4 = res[f"p12_{p12:g}_PP_H4"]
    h3 = res[f"p12_{p12:g}_PP_H3"]
    print(f"p12={p12:g}: H4>=0.8 {int((h4 >= 0.8).sum())}/58; "
          f"H4 0.5-0.8 {int(((h4 >= 0.5) & (h4 < 0.8)).sum())}; "
          f"H3>=0.8 {int((h3 >= 0.8).sum())}/58")
