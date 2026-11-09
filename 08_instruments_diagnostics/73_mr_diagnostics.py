#!/usr/bin/env python3
"""Per-SNP diagnostics for the M4 estimates: Wald ratios, Q contributions,
MR-PRESSO outlier p, and fixed-effect IVW for comparison."""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "mr")


def analyse(file, tag):
    d = pd.read_csv(os.path.join(OUT, file), sep="\t")
    d = d[d["use"]].reset_index(drop=True)
    bx, by, sx, sy = (d["b_exp"].values, d["b_out"].values,
                      d["se_exp"].values, d["se_out"].values)
    w = 1 / sy ** 2
    b_ivw = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se_fe = np.sqrt(1 / np.sum(w * bx ** 2))
    q = w * (by - b_ivw * bx) ** 2
    df = len(bx) - 1
    se_re = se_fe * np.sqrt(max(1, q.sum() / df))
    biv = by / bx
    se_iv = sy / np.abs(bx)
    out_p = []
    for j in range(len(bx)):
        m = np.arange(len(bx)) != j
        bj = np.sum(w[m] * bx[m] * by[m]) / np.sum(w[m] * bx[m] ** 2)
        rss_j = np.sum(w[m] * (by[m] - bj * bx[m]) ** 2)
        out_p.append(stats.chi2.sf(max(q.sum() - rss_j, 0), 1))
    d["wald_or"] = np.exp(biv)
    d["wald_lo"] = np.exp(biv - 1.96 * se_iv)
    d["wald_hi"] = np.exp(biv + 1.96 * se_iv)
    d["q_contrib"] = q
    d["q_pct"] = 100 * q / q.sum()
    d["presso_outlier_p"] = out_p
    d["locus"] = np.where((d["chr"].astype(str) == "1") &
                          (d["pos"].between(169549811 - 1e6, 169549811 + 1e6)), "F5",
                 np.where((d["chr"].astype(str) == "9") &
                          (d["pos"].between(136130000, 136160000)), "ABO",
                 np.where((d["chr"].astype(str) == "6") &
                          (d["pos"].between(25e6, 35e6)), "MHC", "")))
    d.sort_values("q_pct", ascending=False).to_csv(
        os.path.join(OUT, f"diagnostics_{tag}.tsv"), sep="\t", index=False,
        float_format="%.6g")
    print(f"\n=== {tag} === n={len(d)}  IVW fixed OR {np.exp(b_ivw):.3f} "
          f"({np.exp(b_ivw - 1.96 * se_fe):.3f}-{np.exp(b_ivw + 1.96 * se_fe):.3f}), "
          f"P={2 * stats.norm.sf(abs(b_ivw / se_fe)):.3g}; random-effects SE x"
          f"{se_re / se_fe:.2f}")
    print("top Q contributors:")
    print(d.sort_values("q_contrib", ascending=False)
          [["SNP", "chr", "pos", "locus", "wald_or", "wald_lo", "wald_hi",
            "q_pct", "presso_outlier_p"]].head(8).to_string(index=False,
                                                           float_format=lambda v: f"{v:.3g}"))
    print(f"Q total {q.sum():.0f}, df {df}, I2 {max(0, (q.sum() - df) / q.sum()):.3f}; "
          f"SNPs with outlier p<0.05: {(np.array(out_p) < 0.05).sum()}")


analyse("harmonised_forward.tsv", "forward")
analyse("harmonised_reverse.tsv", "reverse")
