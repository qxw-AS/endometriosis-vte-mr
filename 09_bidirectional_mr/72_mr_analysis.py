#!/usr/bin/env python3
"""M4 step 3: bidirectional two-sample MR, implemented from summary statistics.

Estimators: IVW (multiplicative random effects), MR-Egger, weighted median,
weighted mode. Diagnostics: Cochran Q / I^2, I^2_GX, mean/min F, Steiger,
MR-PRESSO global and outlier tests, leave-one-out, and locus-drop sensitivity.
No R / TwoSampleMR available in this environment; formulas follow Bowden 2016,
Hartwig 2017, Verbanck 2018 and the TwoSampleMR implementations.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "mr")
RNG = np.random.default_rng(20261108)
N = {"endo": 20534 + 122622, "vte": 27928 + 472258}


def ivw(bx, by, sy, se_x=None):
    w = 1.0 / sy ** 2
    denom = np.sum(w * bx ** 2)
    beta = np.sum(w * bx * by) / denom
    se = np.sqrt(1.0 / denom)
    q = np.sum(w * (by - beta * bx) ** 2)
    df = len(bx) - 1
    se_rand = se * np.sqrt(max(1.0, q / df)) if df > 0 else se
    p = 2 * stats.norm.sf(abs(beta / se_rand))
    return dict(beta=beta, se=se_rand, p=p, q=q, q_df=df,
                i2=max(0.0, (q - df) / q) if q > 0 else 0.0)


def egger(bx, by, sy):
    X = np.column_stack([np.ones_like(bx), bx])
    w = 1.0 / sy ** 2
    xtwx = X.T @ (X * w[:, None])
    cov = np.linalg.inv(xtwx)
    coef = cov @ (X.T @ (w * by))
    resid = by - X @ coef
    df = len(bx) - 2
    sigma2 = np.sum(w * resid ** 2) / df
    scale = max(1.0, np.sqrt(sigma2))
    se = np.sqrt(np.diag(cov)) * scale
    p_slope = 2 * stats.t.sf(abs(coef[1] / se[1]), df)
    p_int = 2 * stats.t.sf(abs(coef[0] / se[0]), df)
    return dict(beta=coef[1], se=se[1], p=p_slope,
                intercept=coef[0], intercept_se=se[0], intercept_p=p_int)


def weighted_median(biv, w):
    order = np.argsort(biv)
    b, ww = biv[order], w[order]
    cw = np.cumsum(ww) / np.sum(ww)
    j = int(np.searchsorted(cw, 0.5))
    return float(b[min(j, len(b) - 1)])


def weighted_mode(biv, w, bw=0.9):
    dens = np.array([np.sum(w * stats.norm.pdf((bi - biv) / bw)) for bi in biv])
    return float(biv[int(np.argmax(dens))])


def boot_est(bx, by, sx, sy, fn, n_boot=1000):
    vals = []
    k = len(bx)
    for _ in range(n_boot):
        idx = RNG.integers(0, k, k)
        bx_s = RNG.normal(bx[idx], sx[idx])
        by_s = RNG.normal(by[idx], sy[idx])
        biv = by_s / bx_s
        w = 1.0 / sy[idx] ** 2
        try:
            vals.append(fn(biv, w))
        except Exception:
            continue
    vals = np.array(vals)
    return float(np.std(vals, ddof=1)), float(np.percentile(vals, 2.5)), \
        float(np.percentile(vals, 97.5))


def presso(bx, by, sy, n_sim=1000):
    w = 1.0 / sy ** 2
    beta = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    resid = by - beta * bx
    rss_obs = np.sum(resid ** 2 / sy ** 2)
    sim = np.empty(n_sim)
    for i in range(n_sim):
        y = beta * bx + RNG.normal(0, sy)
        bb = np.sum(w * bx * y) / np.sum(w * bx ** 2)
        sim[i] = np.sum((y - bb * bx) ** 2 / sy ** 2)
    p_global = (np.sum(sim >= rss_obs) + 1) / (n_sim + 1)
    out_p = []
    for j in range(len(bx)):
        mask = np.arange(len(bx)) != j
        bj = np.sum(w[mask] * bx[mask] * by[mask]) / np.sum(w[mask] * bx[mask] ** 2)
        rss_j = np.sum(w[mask] * (by[mask] - bj * bx[mask]) ** 2)
        stat = rss_obs - rss_j
        out_p.append(stats.chi2.sf(max(stat, 0), 1))
    return p_global, np.array(out_p)


def steiger(bx, sx, by, sy, n_x, n_y):
    r2x = bx ** 2 / (bx ** 2 + sx ** 2 * n_x)
    r2y = by ** 2 / (by ** 2 + sy ** 2 * n_y)
    j = int(np.argmax(r2x))
    zx = np.arctanh(np.sqrt(min(r2x[j], 0.999999)))
    zy = np.arctanh(np.sqrt(min(r2y[j], 0.999999)))
    se = np.sqrt(1 / (n_x - 3) + 1 / (n_y - 3))
    z = (zx - zy) / se
    return 2 * stats.norm.sf(abs(z)), r2x.sum(), r2y.sum(), j


def run(name, h):
    d = h[h["use"]].reset_index(drop=True)
    bx, by = d["b_exp"].values, d["b_out"].values
    sx, sy = d["se_exp"].values, d["se_out"].values
    res = {"analysis": name, "n_snp": len(d),
           "mean_F": float(np.mean((bx / sx) ** 2)),
           "min_F": float(np.min((bx / sx) ** 2))}
    iv = ivw(bx, by, sy)
    res.update(ivw_beta=iv["beta"], ivw_se=iv["se"], ivw_p=iv["p"],
               ivw_or=np.exp(iv["beta"]),
               ivw_lo=np.exp(iv["beta"] - 1.96 * iv["se"]),
               ivw_hi=np.exp(iv["beta"] + 1.96 * iv["se"]),
               Q=iv["q"], Q_df=iv["q_df"], I2=iv["i2"])
    qgx = np.sum((bx / sx) ** 2)
    res["I2_GX"] = max(0.0, (qgx - (len(bx) - 1)) / qgx)

    eg = egger(bx, by, sy)
    res.update(egger_beta=eg["beta"], egger_se=eg["se"], egger_p=eg["p"],
               egger_or=np.exp(eg["beta"]),
               egger_lo=np.exp(eg["beta"] - 1.96 * eg["se"]),
               egger_hi=np.exp(eg["beta"] + 1.96 * eg["se"]),
               egger_intercept=eg["intercept"], egger_intercept_p=eg["intercept_p"])

    biv = by / bx
    w = 1.0 / sy ** 2
    wm = weighted_median(biv, w)
    wm_se, wm_lo, wm_hi = boot_est(bx, by, sx, sy, weighted_median)
    res.update(wmedian_beta=wm, wmedian_se=wm_se,
               wmedian_or=np.exp(wm), wmedian_lo=np.exp(wm_lo), wmedian_hi=np.exp(wm_hi),
               wmedian_p=2 * stats.norm.sf(abs(wm / wm_se)))
    wmo = weighted_mode(biv, w)
    wmo_se, wmo_lo, wmo_hi = boot_est(bx, by, sx, sy, weighted_mode)
    res.update(wmode_beta=wmo, wmode_se=wmo_se,
               wmode_or=np.exp(wmo), wmode_lo=np.exp(wmo_lo), wmode_hi=np.exp(wmo_hi),
               wmode_p=2 * stats.norm.sf(abs(wmo / wmo_se)))

    pg, outp = presso(bx, by, sy)
    j = int(np.argmin(outp))
    res.update(presso_global_p=pg, presso_min_outlier_p=float(outp[j]),
               presso_worst_snp=d.loc[j, "SNP"], n_presso_outlier=int(np.sum(outp < 0.05)))

    st_p, r2x, r2y, js = steiger(bx, sx, by, sy, N[d["exposure"][0]], N[d["outcome"][0]])
    res.update(steiger_p=st_p, r2_exposure=r2x, r2_outcome=r2y,
               steiger_snp=d.loc[js, "SNP"], steiger_pass=st_p < 0.001)

    # leave-one-out (IVW)
    loo = []
    for j in range(len(bx)):
        mask = np.arange(len(bx)) != j
        r = ivw(bx[mask], by[mask], sy[mask])
        loo.append((d.loc[j, "SNP"], np.exp(r["beta"]), r["p"]))
    loo = pd.DataFrame(loo, columns=["SNP", "OR", "p"])
    loo["shift_from_full"] = loo["OR"] - res["ivw_or"]
    loo = loo.reindex(loo["shift_from_full"].abs().sort_values(ascending=False).index)
    res["loo_top"] = "; ".join(f"{r.SNP} {r.OR:.2f} (shift {r.shift_from_full:+.2f})"
                               for r in loo.head(3).itertuples())
    return res, d, loo


def main():
    out = []
    for name, file in (("forward endo->VTE", "harmonised_forward.tsv"),
                       ("reverse VTE->endo", "harmonised_reverse.tsv")):
        h = pd.read_csv(os.path.join(OUT, file), sep="\t")
        r, d, loo = run(name, h)
        r["direction"] = "forward" if "forward" in name else "reverse"
        out.append(r)
        loo.to_csv(os.path.join(OUT, f"loo_{r['direction']}.tsv"), sep="\t", index=False)

        # --- sensitivity: drop the F5 (rs6025 +/-1 Mb on chr1) locus
        f5 = d[(d["chr"].astype(str) == "1") &
               (d["pos"].astype(float).between(169549811 - 1_000_000, 169549811 + 1_000_000))]
        if len(f5):
            r2, _, _ = run(name + " | drop F5 locus", d[~d["SNP"].isin(f5["SNP"])])
            r2["direction"] = r["direction"] + "_dropF5"
            out.append(r2)
        # --- sensitivity: drop MHC (chr6:25-35 Mb)
        mhc = d[(d["chr"].astype(str) == "6") & (d["pos"].astype(float).between(25e6, 35e6))]
        if len(mhc):
            r3, _, _ = run(name + " | drop MHC", d[~d["SNP"].isin(mhc["SNP"])])
            r3["direction"] = r["direction"] + "_dropMHC"
            out.append(r3)

    res = pd.DataFrame(out)
    res.to_csv(os.path.join(OUT, "mr_results.tsv"), sep="\t", index=False,
               float_format="%.6g")
    pd.set_option("display.width", 200)
    cols = ["analysis", "n_snp", "mean_F", "ivw_or", "ivw_lo", "ivw_hi", "ivw_p",
            "egger_or", "egger_p", "egger_intercept_p", "wmedian_or", "wmedian_p",
            "wmode_or", "wmode_p", "Q", "I2", "I2_GX", "presso_global_p",
            "steiger_p", "steiger_pass"]
    print(res[cols].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("\nleave-one-out top movers:")
    for r in out:
        print(f"  {r['analysis']}: {r.get('loo_top', '')}")


if __name__ == "__main__":
    main()
