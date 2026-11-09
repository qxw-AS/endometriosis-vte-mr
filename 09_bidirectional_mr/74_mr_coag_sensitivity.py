#!/usr/bin/env python3
"""M4 sensitivity: remove the protocol-defined 9-gene coagulation set (M6.1)
from the instrument set, then re-estimate.

Exclusion windows are built from the GRCh38 sentinel SNP of each gene
(out/fig3b_sentinel_or.tsv) +/-1 Mb, which also captures proxies such as
rs1209731 (r2 = 0.77 with rs6025 in 1000G EUR) that clumping kept because
rs6025 itself is not genome-wide significant for endometriosis.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "mr")
sent = pd.read_csv(os.path.join(ROOT, "out", "fig3b_sentinel_or.tsv"), sep="\t")
print(sent[["gene", "SNP", "chrom", "pos"]].to_string(index=False))


def ivw(bx, by, sy):
    w = 1 / sy ** 2
    b = (w * bx * by).sum() / (w * bx ** 2).sum()
    se = np.sqrt(1 / (w * bx ** 2).sum())
    q = (w * (by - b * bx) ** 2).sum()
    df = len(bx) - 1
    se_re = se * np.sqrt(max(1, q / df))
    return dict(or_=np.exp(b), lo=np.exp(b - 1.96 * se_re), hi=np.exp(b + 1.96 * se_re),
                p=2 * stats.norm.sf(abs(b / se_re)), q=q, df=df,
                i2=max(0, (q - df) / q) if q > 0 else 0)


def egger(bx, by, sy):
    X = np.column_stack([np.ones_like(bx), bx])
    w = 1 / sy ** 2
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    coef = cov @ (X.T @ (w * by))
    resid = by - X @ coef
    s2 = (w * resid ** 2).sum() / (len(bx) - 2)
    se = np.sqrt(np.diag(cov) * max(1, s2))
    return dict(or_=np.exp(coef[1]), lo=np.exp(coef[1] - 1.96 * se[1]),
                hi=np.exp(coef[1] + 1.96 * se[1]),
                p=2 * stats.t.sf(abs(coef[1] / se[1]), len(bx) - 2),
                int_p=2 * stats.t.sf(abs(coef[0] / se[0]), len(bx) - 2))


def wmedian(biv, w):
    o = np.argsort(biv)
    cw = np.cumsum(w[o]) / w[o].sum()
    j = np.searchsorted(cw, 0.5)
    return biv[o][min(j, len(biv) - 1)]


def wmode(biv, w, bw=0.9):
    dens = [np.sum(w * stats.norm.pdf((b - biv) / bw)) for b in biv]
    return biv[int(np.argmax(dens))]


rng = np.random.default_rng(7)
res = []
for tag, file in (("forward endo->VTE", "harmonised_forward.tsv"),
                  ("reverse VTE->endo", "harmonised_reverse.tsv")):
    d = pd.read_csv(os.path.join(OUT, file), sep="\t")
    d = d[d["use"]].copy()
    d["chr"] = d["chr"].astype(str).str.replace("chr", "")
    d["pos"] = d["pos"].astype(float)
    hit = pd.Series(False, index=d.index)
    for _, r in sent.iterrows():
        hit |= (d["chr"] == str(r["chrom"])) & \
               (d["pos"].between(r["pos"] - 1_000_000, r["pos"] + 1_000_000))
    rem = d[~hit]
    print(f"\n{tag}: {len(d)} -> {len(rem)} after removing {int(hit.sum())} "
          f"coagulation-locus SNPs ({', '.join(d.loc[hit, 'SNP'])})")
    bx, by, sy = rem["b_exp"].values, rem["b_out"].values, rem["se_out"].values
    iv, eg = ivw(bx, by, sy), egger(bx, by, sy)
    biv, w = by / bx, 1 / sy ** 2
    wm, wmo = wmedian(biv, w), wmode(biv, w)
    bs = []
    for _ in range(1000):
        idx = rng.integers(0, len(bx), len(bx))
        b = rng.normal(by[idx], rem["se_out"].values[idx]) / \
            rng.normal(bx[idx], rem["se_exp"].values[idx])
        bs.append((wmedian(b, 1 / rem["se_out"].values[idx] ** 2),
                   wmode(b, 1 / rem["se_out"].values[idx] ** 2)))
    bs = np.array(bs)
    res.append(dict(analysis=tag, n_snp=len(rem),
                    ivw_or=iv["or_"], ivw_lo=iv["lo"], ivw_hi=iv["hi"], ivw_p=iv["p"],
                    egger_or=eg["or_"], egger_p=eg["p"], egger_int_p=eg["int_p"],
                    wmedian_or=np.exp(wm),
                    wmedian_p=2 * stats.norm.sf(abs(wm / bs[:, 0].std(ddof=1))),
                    wmode_or=np.exp(wmo),
                    wmode_p=2 * stats.norm.sf(abs(wmo / bs[:, 1].std(ddof=1))),
                    Q=iv["q"], Q_df=iv["df"], I2=iv["i2"]))
out = pd.DataFrame(res)
out.to_csv(os.path.join(OUT, "mr_coag_excluded.tsv"), sep="\t", index=False,
           float_format="%.6g")
print("\n" + out.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
