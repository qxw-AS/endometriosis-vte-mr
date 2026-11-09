#!/usr/bin/env python3
"""Weighted-median estimates for the severity rows that lack one.

Recomputes PE-specific and DVT-specific estimates after removal of the nine
coagulation-gene windows (+/-1 Mb), and verifies that the weighted-median
values already stored for the main severity rows are reproducible.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out")
RNG = np.random.default_rng(20261109)
sent = pd.read_csv(os.path.join(OUT, "fig3b_sentinel_or.tsv"), sep="\t")


def coag_mask(chrom, pos):
    m = pd.Series(False, index=chrom.index)
    for _, s in sent.iterrows():
        m |= (chrom == int(s["chrom"])) & \
             pos.between(int(s["pos"]) - 1_000_000, int(s["pos"]) + 1_000_000)
    return m


def wmedian(biv, w):
    o = np.argsort(biv)
    cw = np.cumsum(w[o]) / np.sum(w)
    j = int(np.searchsorted(cw, 0.5))
    return float(biv[o][min(j, len(biv) - 1)])


def boot(bx, by, sx, sy, n_boot=1000):
    vals = []
    k = len(bx)
    for _ in range(n_boot):
        idx = RNG.integers(0, k, k)
        bx_s = RNG.normal(bx[idx], sx[idx])
        by_s = RNG.normal(by[idx], sy[idx])
        vals.append(wmedian(by_s / bx_s, 1 / sy[idx] ** 2))
    v = np.array(vals)
    return float(np.std(v, ddof=1)), float(np.percentile(v, 2.5)), \
        float(np.percentile(v, 97.5))


def estimate(g, expo):
    bx = g[f"beta_{expo}"].values.astype(float)
    by = g["beta_N14_ENDOMETRIOSIS"].values.astype(float)
    sx = g[f"se_{expo}"].values.astype(float)
    sy = g["se_N14_ENDOMETRIOSIS"].values.astype(float)
    b = wmedian(by / bx, 1 / sy ** 2)
    se, lo, hi = boot(bx, by, sx, sy)
    return dict(n_snp=len(g), wmedian_or=np.exp(b), wmedian_lo=np.exp(lo),
                wmedian_hi=np.exp(hi),
                wmedian_p=float(2 * stats.norm.sf(abs(b / se))))


def main():
    s = pd.read_csv(os.path.join(OUT, "sev", "severity_instruments.tsv"), sep="\t")
    s["chr"] = pd.to_numeric(s["chrom"], errors="coerce").astype(int)
    s["coag"] = coag_mask(s["chr"], s["pos"])
    PE, DVT = "I9_PULMEMB", "I9_PHLETHROMBDVTLOW"
    pe_only = (s[f"p_{PE}"] < 5e-8) & (s[f"p_{DVT}"] >= 5e-8)
    dvt_only = (s[f"p_{DVT}"] < 5e-8) & (s[f"p_{PE}"] >= 5e-8)

    rows = []
    for lab, mask, expo in (
        ("PE-specific minus coag", pe_only & ~s["coag"], PE),
        ("DVT-specific minus coag", dvt_only & ~s["coag"], DVT),
    ):
        r = estimate(s[mask], expo)
        r["analysis"] = lab
        rows.append(r)
        print(f"{lab}: n={r['n_snp']} wmedian OR={r['wmedian_or']:.3f} "
              f"({r['wmedian_lo']:.3f}-{r['wmedian_hi']:.3f}) P={r['wmedian_p']:.3g}")

    out = pd.DataFrame(rows)[["analysis", "n_snp", "wmedian_or", "wmedian_lo",
                              "wmedian_hi", "wmedian_p"]]
    out.to_csv(os.path.join(OUT, "sev", "severity_wmedian_extra.tsv"),
               sep="\t", index=False, float_format="%.6g")

    # verify stored weighted medians are reproducible (point estimates only)
    mr = pd.read_csv(os.path.join(OUT, "sev", "severity_mr.tsv"), sep="\t").set_index("analysis")
    sp = pd.read_csv(os.path.join(OUT, "sev", "severity_specific.tsv"), sep="\t").set_index("analysis")
    for label, expected in (("PE->endo", mr.loc["PE->endo", "wmedian_or"]),
                            ("DVT->endo", mr.loc["DVT->endo", "wmedian_or"]),
                            ("PE_specific->endo", sp.loc["PE_specific->endo", "wmedian_or"]),
                            ("DVT_specific->endo", sp.loc["DVT_specific->endo", "wmedian_or"])):
        print(f"  stored weighted median {label}: {expected:.6f}")


if __name__ == "__main__":
    main()
