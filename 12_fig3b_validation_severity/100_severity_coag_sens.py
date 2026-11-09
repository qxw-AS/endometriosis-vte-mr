#!/usr/bin/env python3
"""Fill remaining numbers for Table III and the severity analysis.

1. MR-PRESSO global p for the forward coagulation-excluded set (34 SNPs).
2. Severity estimates after removing the 9-gene coagulation windows:
   PE-specific and DVT-specific -> endometriosis.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RNG = np.random.default_rng(20261108)
sent = pd.read_csv(os.path.join(ROOT, "out", "fig3b_sentinel_or.tsv"), sep="\t")


def ivw(bx, by, sy):
    w = 1 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se = np.sqrt(1 / np.sum(w * bx ** 2))
    q = float(np.sum(w * (by - b * bx) ** 2))
    dfree = len(bx) - 1
    se_re = se * np.sqrt(max(1.0, q / dfree))
    return dict(beta=b, se=se_re, p=2 * stats.norm.sf(abs(b / se_re)),
                q=q, df=dfree, i2=max(0.0, (q - dfree) / q) if q > 0 else 0.0)


def presso(bx, by, sy, n_sim=1000):
    w = 1 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    rss = float(np.sum(w * (by - b * bx) ** 2))
    sim = np.empty(n_sim)
    for i in range(n_sim):
        y = b * bx + RNG.normal(0, sy)
        bb = np.sum(w * bx * y) / np.sum(w * bx ** 2)
        sim[i] = np.sum(w * (y - bb * bx) ** 2)
    return (np.sum(sim >= rss) + 1) / (n_sim + 1)


def coag_mask(chrom, pos):
    chrom = pd.to_numeric(chrom, errors="coerce").astype(int)
    pos = pd.to_numeric(pos, errors="coerce").astype(int)
    m = pd.Series(False, index=chrom.index)
    for _, s in sent.iterrows():
        m |= (chrom == int(s["chrom"])) & \
            pos.between(int(s["pos"]) - 1_000_000, int(s["pos"]) + 1_000_000)
    return m


# 1. forward PRESSO
d = pd.read_csv(os.path.join(ROOT, "out", "mr", "harmonised_forward.tsv"),
                sep="\t")
d = d[d["use"]].reset_index(drop=True)
m = coag_mask(d["chr"], d["pos"])
r = d[~m]
pr = presso(r["b_exp"].values, r["b_out"].values, r["se_out"].values)
iv = ivw(r["b_exp"].values, r["b_out"].values, r["se_out"].values)
print(f"forward coag-excluded n={len(r)} OR={np.exp(iv['beta']):.3f} "
      f"P={iv['p']:.3g} I2={iv['i2']:.2f} PRESSO_global_p={pr:.4f}")

# 2. severity sensitivity
s = pd.read_csv(os.path.join(ROOT, "out", "sev", "severity_instruments.tsv"),
                sep="\t")
s["chr"] = pd.to_numeric(s["chrom"], errors="coerce").astype(int)
s["coag"] = coag_mask(s["chr"], s["pos"])
PE, DVT, ENDO = "I9_PULMEMB", "I9_PHLETHROMBDVTLOW", "N14_ENDOMETRIOSIS"
pe_only = (s[f"p_{PE}"] < 5e-8) & (s[f"p_{DVT}"] >= 5e-8)
dvt_only = (s[f"p_{DVT}"] < 5e-8) & (s[f"p_{PE}"] >= 5e-8)
for lab, m in (("PE_specific", pe_only), ("DVT_specific", dvt_only)):
    for tag, mm in (("all", m), ("minus_coag", m & ~s["coag"])):
        g = s[mm]
        if len(g) < 3:
            print(f"{lab} {tag}: n={len(g)} too few")
            continue
        bx = g[f"beta_{PE if lab.startswith('PE') else DVT}"].values
        by = g[f"beta_{ENDO}"].values
        sy = g[f"se_{ENDO}"].values
        iv = ivw(bx, by, sy)
        print(f"{lab} {tag}: n={len(g)} OR={np.exp(iv['beta']):.3f} "
              f"({np.exp(iv['beta']-1.96*iv['se']):.3f}-"
              f"{np.exp(iv['beta']+1.96*iv['se']):.3f}) P={iv['p']:.3g} "
              f"I2={iv['i2']:.2f}")

# collinearity between subtype effects
print(f"\ncorr(beta_PE, beta_DVT) over 73 leads: "
      f"{np.corrcoef(s[f'beta_{PE}'], s[f'beta_{DVT}'])[0,1]:.3f}")
