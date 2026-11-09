#!/usr/bin/env python3
"""M3 stage 3: coloc.abf (Wakefield ABF) in pure Python.

Implements the standard single-causal-variant coloc.abf model:
  lABF_i = 0.5 * (log(V_i) - log(V_i + W) + z_i^2 * W / (V_i + W))
  H0: 1
  H1: p1 * sum_i ABF1_i
  H2: p2 * sum_i ABF2_i
  H3: p1 * p2 * sum_i ABF1_i * sum_j ABF2_j
  H4: p12 * sum_i ABF1_i * ABF2_i
with W = sdY^2, priors p1 = p2 = 1e-4, p12 = 1e-5.

Primary sdY = 0.15 per protocol; sensitivity 0.10 / 0.20 / 0.30.
LAVA (M2) was not run, so regions are the protocol fallback lead +/- 500 kb.
"""
import os

import numpy as np
import pandas as pd
from scipy.special import logsumexp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m3")

P1, P2, P12 = 1e-4, 1e-4, 1e-5
LOG_P1, LOG_P2, LOG_P12 = np.log(P1), np.log(P2), np.log(P12)
SDY_MAIN = 0.15
SDY_GRID = [0.10, 0.15, 0.20, 0.30]
MIN_MAF = 0.01

regions = pd.read_csv(os.path.join(OUT, "candidate_regions.tsv"), sep="\t",
                      dtype={"chr": str})


def labf(beta, se, sdy):
    v = se ** 2
    w = sdy ** 2
    z = beta / se
    return 0.5 * (np.log(v) - np.log(v + w) + z ** 2 * w / (v + w))


def posteriors(l1, l2, p1=P1, p2=P2, p12=P12):
    lh = np.array([
        0.0,
        np.log(p1) + logsumexp(l1),
        np.log(p2) + logsumexp(l2),
        np.log(p1) + np.log(p2) + logsumexp(l1) + logsumexp(l2),
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
print(f"region SNPs: endo {len(endo):,}, vte {len(vte):,}")

key = ["chr", "pos", "ref", "alt"]
m = endo.merge(vte, on=key, suffixes=("_e", "_v"))
print(f"shared SNP keys: {len(m):,}")

res, sens = [], []
for _, r in regions.iterrows():
    g = m[(m["chr"] == r["chr"]) & (m["pos"] >= r["start"]) & (m["pos"] <= r["end"])]
    g = g[(g["maf_e"] >= MIN_MAF) & (g["maf_v"] >= MIN_MAF)].copy()
    out = dict(region=r["region"], chr=r["chr"], start=r["start"], end=r["end"],
               n_leads_endo=r["n_endo"], n_leads_vte=r["n_vte"],
               mhc=bool(r["mhc"]), n_snp=len(g))
    if len(g) < 2:
        out.update({f"PP.H{i}": np.nan for i in range(5)})
        out["note"] = "too few shared SNPs"
        res.append(out)
        continue
    # allele key already identical, so betas are on the same ALT allele
    pp = {}
    for sdy in SDY_GRID:
        l1 = labf(g["beta_e"].values, g["se_e"].values, sdy)
        l2 = labf(g["beta_v"].values, g["se_v"].values, sdy)
        pp[sdy] = posteriors(l1, l2)
        sens.append(dict(region=r["region"], sdy=sdy,
                         PP_H0=pp[sdy][0], PP_H1=pp[sdy][1], PP_H2=pp[sdy][2],
                         PP_H3=pp[sdy][3], PP_H4=pp[sdy][4]))
    l1 = labf(g["beta_e"].values, g["se_e"].values, SDY_MAIN)
    l2 = labf(g["beta_v"].values, g["se_v"].values, SDY_MAIN)
    p = pp[SDY_MAIN]
    i1 = int(np.argmax(l1))
    i2 = int(np.argmax(l2))
    i4 = int(np.argmax(l1 + l2))
    top = g.iloc[i4]
    out.update({f"PP.H{i}": p[i] for i in range(5)})
    def label(row):
        s = row["rsids_e"]
        return s if isinstance(s, str) else f"{row['chr']}:{int(row['pos'])}"
    out.update(dict(
        top_endo=label(g.iloc[i1]), top_vte=label(g.iloc[i2]),
        top_shared=label(top), top_shared_pos=int(top["pos"]),
        beta_endo=float(top["beta_e"]), beta_vte=float(top["beta_v"]),
        p_endo=float(top["p_e"]), p_vte=float(top["p_v"]),
        direction_concordant=bool(np.sign(top["beta_e"]) == np.sign(top["beta_v"])),
        same_top=bool(g.iloc[i1]["rsids_e"] == g.iloc[i2]["rsids_v"]),
        note="",
    ))
    res.append(out)

res = pd.DataFrame(res)
res.to_csv(os.path.join(OUT, "coloc_results.tsv"), sep="\t", index=False,
           float_format="%.6g")
sens = pd.DataFrame(sens)
sens.to_csv(os.path.join(OUT, "coloc_sensitivity.tsv"), sep="\t", index=False,
            float_format="%.6g")

pd.set_option("display.width", 250)
print("\n=== coloc.abf (sdY=0.15) ===")
print(res[["region", "chr", "n_snp", "mhc", "PP.H3", "PP.H4", "top_shared",
           "beta_endo", "beta_vte", "direction_concordant"]]
      .sort_values("PP.H4", ascending=False).head(25)
      .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print(f"\nPP.H4 >= 0.8: {(res['PP.H4'] >= 0.8).sum()} of {len(res)} regions "
      f"(MHC flagged {(res['mhc'] & (res['PP.H4'] >= 0.8)).sum()})")
print(f"PP.H4 0.5-0.8: {((res['PP.H4'] >= 0.5) & (res['PP.H4'] < 0.8)).sum()}")
print("\n=== sdY sensitivity (PP.H4) ===")
piv = sens.pivot(index="region", columns="sdy", values="PP_H4")
print(piv.round(3).to_string())
