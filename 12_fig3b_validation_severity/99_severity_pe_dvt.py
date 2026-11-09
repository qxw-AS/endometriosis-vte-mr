#!/usr/bin/env python3
"""Step 3: reverse signal stratified by VTE severity (PE vs DVT).

Primary analysis
- Build LD-independent lead sets for pulmonary embolism (I9_PULMEMB) and
  deep-vein thrombosis (I9_PHLETHROMBDVTLOW) from FinnGen R13 pheweb
  manhattan peaks (peak=True, p<5e-8).
- Estimate the causal effect of genetic liability to each subtype on
  endometriosis (N14_ENDOMETRIOSIS) with IVW / MR-Egger / weighted median.
- Compare the two subtype effects formally with an MVMR of PE and DVT on
  the union of subtype instruments (difference b_PE - b_DVT, covariance
  taken from the joint fit).

Secondary analysis
- Split the 44 VTE instruments by whether they are also PE or DVT leads
  (p<5e-8) and repeat VTE->endo in each subset.

Caveat: all traits are FinnGen R13, so the design is one-sample.
"""
import gzip
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import log, sqrt

import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "sev")
CACHE = os.path.join(OUT, "cache")
os.makedirs(CACHE, exist_ok=True)

PE, DVT, ENDO, VTE = ("I9_PULMEMB", "I9_PHLETHROMBDVTLOW",
                      "N14_ENDOMETRIOSIS", "I9_VTE")
PHENOS = [PE, DVT, ENDO, VTE]
BASE = "https://r13.finngen.fi/api/variant"
RNG = np.random.default_rng(20261108)


def fetch_peaks(pheno):
    url = f"https://r13.finngen.fi/api/manhattan/pheno/{pheno}"
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=180)
            if r.status_code == 200:
                raw = r.content
                if raw[:2] == b"\x1f\x8b":
                    raw = gzip.decompress(raw)
                j = json.loads(raw)
                d = pd.DataFrame(j["unbinned_variants"])
                d = d[(d["peak"] == True) & (d["pval"] < 5e-8)].copy()
                d["pheno"] = pheno
                return d
        except Exception as exc:
            print(f"peaks {pheno} attempt {attempt}: {exc}", flush=True)
            time.sleep(2)
    raise RuntimeError(f"could not fetch peaks for {pheno}")


def get_one(vid, pheno):
    fn = os.path.join(CACHE, f"{vid}__{pheno}.json")
    if os.path.exists(fn) and os.path.getsize(fn) > 0:
        try:
            return vid, pheno, json.load(open(fn))
        except Exception:
            pass
    url = f"{BASE}/{vid}/{pheno}"
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                j = r.json()
                json.dump(j, open(fn, "w"))
                return vid, pheno, j
            if r.status_code == 404:
                json.dump({"results": None}, open(fn, "w"))
                return vid, pheno, {"results": None}
        except Exception:
            time.sleep(0.5 * (attempt + 1))
    return vid, pheno, {"results": None}


def parse(j):
    r = (j or {}).get("results")
    if not r:
        return None
    beta, p, mlogp = r.get("beta"), r.get("pval"), r.get("mlogp")
    if beta is None or p is None:
        return None
    beta, p = float(beta), float(p)
    if p > 0:
        z = stats.norm.isf(p / 2)
    elif mlogp:
        z = sqrt(2 * float(mlogp) * log(10))
    else:
        return None
    if not np.isfinite(z) or z <= 0 or beta == 0:
        return None
    return dict(beta=beta, se=abs(beta) / z, p=p, maf=r.get("maf"),
                n_case=r.get("n_case"), n_control=r.get("n_control"))


def ivw(bx, by, sy):
    w = 1 / sy ** 2
    denom = np.sum(w * bx ** 2)
    b = np.sum(w * bx * by) / denom
    se = np.sqrt(1 / denom)
    q = float(np.sum(w * (by - b * bx) ** 2))
    dfree = len(bx) - 1
    se_re = se * np.sqrt(max(1.0, q / dfree)) if dfree > 0 else se
    return dict(beta=b, se=se_re, p=2 * stats.norm.sf(abs(b / se_re)),
                q=q, q_df=dfree,
                i2=max(0.0, (q - dfree) / q) if q > 0 else 0.0)


def egger(bx, by, sy):
    X = np.column_stack([np.ones_like(bx), bx])
    w = 1 / sy ** 2
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    coef = cov @ (X.T @ (w * by))
    resid = by - X @ coef
    dfree = len(bx) - 2
    s2 = float(np.sum(w * resid ** 2)) / dfree
    se = np.sqrt(np.diag(cov)) * max(1.0, sqrt(s2))
    return dict(beta=coef[1], se=se[1],
                p=2 * stats.t.sf(abs(coef[1] / se[1]), dfree),
                intercept=coef[0],
                intercept_p=2 * stats.t.sf(abs(coef[0] / se[0]), dfree))


def weighted_median(biv, w):
    o = np.argsort(biv)
    cw = np.cumsum(w[o]) / np.sum(w)
    j = int(np.searchsorted(cw, 0.5))
    return float(biv[o][min(j, len(biv) - 1)])


def boot_ci(bx, by, sx, sy, fn, n_boot=1000):
    vals = []
    k = len(bx)
    for _ in range(n_boot):
        idx = RNG.integers(0, k, k)
        bx_s = RNG.normal(bx[idx], sx[idx])
        by_s = RNG.normal(by[idx], sy[idx])
        try:
            vals.append(fn(by_s / bx_s, 1 / sy[idx] ** 2))
        except Exception:
            pass
    v = np.array(vals)
    return float(np.std(v, ddof=1)), float(np.percentile(v, 2.5)), \
        float(np.percentile(v, 97.5))


def presso(bx, by, sy, n_sim=1000):
    w = 1 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    rss = float(np.sum(w * (by - b * bx) ** 2))
    sim = np.empty(n_sim)
    for i in range(n_sim):
        y = b * bx + RNG.normal(0, sy)
        bb = np.sum(w * bx * y) / np.sum(w * bx ** 2)
        sim[i] = np.sum(w * (y - bb * bx) ** 2)
    p = (np.sum(sim >= rss) + 1) / (n_sim + 1)
    outp = []
    for j in range(len(bx)):
        m = np.arange(len(bx)) != j
        bj = np.sum(w[m] * bx[m] * by[m]) / np.sum(w[m] * bx[m] ** 2)
        rss_j = float(np.sum(w[m] * (by[m] - bj * bx[m]) ** 2))
        outp.append(stats.chi2.sf(max(rss - rss_j, 0), 1))
    return p, np.array(outp)


def run_mr(d, exp_col, out_col, label):
    d = d.reset_index(drop=True)
    bx = d[f"beta_{exp_col}"].values.astype(float)
    by = d[f"beta_{out_col}"].values.astype(float)
    sx = d[f"se_{exp_col}"].values.astype(float)
    sy = d[f"se_{out_col}"].values.astype(float)
    r = dict(analysis=label, n_snp=len(d),
             mean_F=float(np.mean((bx / sx) ** 2)),
             min_F=float(np.min((bx / sx) ** 2)))
    iv = ivw(bx, by, sy)
    r.update(ivw_beta=iv["beta"], ivw_se=iv["se"], ivw_p=iv["p"],
             ivw_or=np.exp(iv["beta"]),
             ivw_lo=np.exp(iv["beta"] - 1.96 * iv["se"]),
             ivw_hi=np.exp(iv["beta"] + 1.96 * iv["se"]),
             Q=iv["q"], Q_df=iv["q_df"], I2=iv["i2"])
    if len(bx) >= 3:
        eg = egger(bx, by, sy)
        r.update(egger_beta=eg["beta"], egger_or=np.exp(eg["beta"]),
                 egger_lo=np.exp(eg["beta"] - 1.96 * eg["se"]),
                 egger_hi=np.exp(eg["beta"] + 1.96 * eg["se"]),
                 egger_p=eg["p"], egger_intercept=eg["intercept"],
                 egger_intercept_p=eg["intercept_p"])
        w = 1 / sy ** 2
        wm = weighted_median(by / bx, w)
        se_wm, lo_wm, hi_wm = boot_ci(bx, by, sx, sy, weighted_median)
        r.update(wmedian_or=np.exp(wm), wmedian_lo=np.exp(lo_wm),
                 wmedian_hi=np.exp(hi_wm),
                 wmedian_p=2 * stats.norm.sf(abs(wm / se_wm)))
        pg, _ = presso(bx, by, sy)
        r["presso_global_p"] = pg
    biv = by / bx
    r["wald_below_1"] = int((np.exp(biv) < 1).sum())
    r["wald_dir_p"] = stats.binomtest(int((np.exp(biv) < 1).sum()),
                                      len(biv), 0.5).pvalue
    return r


def mvmr(df, exposures, out_col):
    X = df[[f"beta_{e}" for e in exposures]].values.astype(float)
    y = df[f"beta_{out_col}"].values.astype(float)
    sy = df[f"se_{out_col}"].values.astype(float)
    w = 1 / sy ** 2
    xtwx = X.T @ (X * w[:, None])
    beta = np.linalg.solve(xtwx, X.T @ (w * y))
    resid = y - X @ beta
    dfree = len(y) - len(exposures)
    s2 = max(1.0, float(np.sum(w * resid ** 2)) / dfree)
    cov = np.linalg.inv(xtwx) * s2
    se = np.sqrt(np.diag(cov))
    out = {}
    for k, e in enumerate(exposures):
        out[e] = dict(beta=beta[k], se=se[k], p=2 * stats.norm.sf(abs(beta[k] / se[k])),
                      or_=np.exp(beta[k]),
                      lo=np.exp(beta[k] - 1.96 * se[k]),
                      hi=np.exp(beta[k] + 1.96 * se[k]))
    # conditional F
    for k, e in enumerate(exposures):
        others = [j for j in range(len(exposures)) if j != k]
        Xo = np.column_stack([np.ones(len(y)), X[:, others]])
        coef = np.linalg.solve(Xo.T @ (Xo * w[:, None]), Xo.T @ (w * X[:, k]))
        rr = X[:, k] - Xo @ coef
        out[e]["cond_F"] = float(np.mean(rr ** 2 /
                                         df[f"se_{e}"].values.astype(float) ** 2))
    if len(exposures) == 2:
        diff = beta[0] - beta[1]
        se_diff = sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1])
        out["diff"] = dict(beta=diff, se=se_diff,
                           p=2 * stats.norm.sf(abs(diff / se_diff)),
                           or_=np.exp(diff),
                           lo=np.exp(diff - 1.96 * se_diff),
                           hi=np.exp(diff + 1.96 * se_diff))
    return out


def main():
    # ---------- subtype lead sets ----------
    pe = fetch_peaks(PE)
    dvt = fetch_peaks(DVT)
    print(f"PE peaks (p<5e-8): {len(pe)}; DVT peaks (p<5e-8): {len(dvt)}",
          flush=True)
    pe["subtype"] = "PE"
    dvt["subtype"] = "DVT"
    lead = pd.concat([pe, dvt], ignore_index=True)
    lead = lead.drop_duplicates(["chrom", "pos", "ref", "alt"])
    lead["vid"] = (lead["chrom"].astype(str) + "-" + lead["pos"].astype(str) +
                   "-" + lead["ref"].astype(str) + "-" + lead["alt"].astype(str))

    def alt_of(v):
        p = v.split("-")
        return "-".join(p[:2] + [p[3], p[2]])

    st2ph = {"PE": PE, "DVT": DVT}
    tasks = []
    for _, r in lead.iterrows():
        tasks.append((r["vid"], st2ph[r["subtype"]]))
        tasks.append((alt_of(r["vid"]), st2ph[r["subtype"]]))
    pres = {}
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs = [ex.submit(get_one, v, ph) for v, ph in tasks]
        for f in as_completed(futs):
            vid, ph, j = f.result()
            pres[(vid, ph)] = j
    chosen = {}
    for _, r in lead.iterrows():
        v, ph = r["vid"], st2ph[r["subtype"]]
        if (pres.get((v, ph)) or {}).get("results"):
            chosen[v] = (v, False)
        elif (pres.get((alt_of(v), ph)) or {}).get("results"):
            chosen[v] = (alt_of(v), True)
    print(f"resolved leads {len(chosen)}/{len(lead)}", flush=True)

    rest = [(c, ph) for c, _ in chosen.values() for ph in PHENOS]
    res = dict(pres)
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs = [ex.submit(get_one, v, ph) for v, ph in rest]
        for f in as_completed(futs):
            vid, ph, j = f.result()
            res[(vid, ph)] = j

    rows = []
    for _, r in lead.iterrows():
        if r["vid"] not in chosen:
            continue
        cv, flip = chosen[r["vid"]]
        row = dict(SNP=(r.get("rsids") or f"{r['chrom']}:{r['pos']}"),
                   chrom=r["chrom"], pos=r["pos"], subtype=r["subtype"],
                   manhattan_beta=float(r["beta"]),
                   manhattan_p=float(r["pval"]))
        ok = True
        for ph in PHENOS:
            pr = parse(res.get((cv, ph)))
            if pr is None:
                ok = False
                break
            row[f"beta_{ph}"] = pr["beta"] * (-1 if flip else 1)
            row[f"se_{ph}"] = pr["se"]
            row[f"p_{ph}"] = pr["p"]
        if ok:
            rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(OUT, "severity_instruments.tsv"), sep="\t",
             index=False, float_format="%.6g")
    print(f"\ninstruments with complete effects: {len(d)} "
          f"(PE {int((d.subtype=='PE').sum())}, "
          f"DVT {int((d.subtype=='DVT').sum())})")
    st2ph = {"PE": PE, "DVT": DVT}
    print("manhattan vs fetched beta (correlation, slope):")
    for st, ph in st2ph.items():
        g = d[d["subtype"] == st]
        if len(g) > 2:
            cor = np.corrcoef(g["manhattan_beta"], g[f"beta_{ph}"])[0, 1]
            slope = np.polyfit(g["manhattan_beta"], g[f"beta_{ph}"], 1)[0]
            print(f"  {st}: n={len(g)} r={cor:.4f} slope={slope:.4f}")

    # ---------- primary: subtype -> endo ----------
    res_rows = []
    for st, ph in (("PE", PE), ("DVT", DVT)):
        g = d[(d["subtype"] == st) & (d[f"p_{ph}"] < 5e-8)].copy()
        print(f"\n{st} instruments for {st}->endo: {len(g)}", flush=True)
        if len(g) >= 3:
            res_rows.append(run_mr(g, ph, ENDO, f"{st}->endo"))
    if res_rows:
        rr = pd.DataFrame(res_rows)
        rr.to_csv(os.path.join(OUT, "severity_mr.tsv"), sep="\t", index=False,
                  float_format="%.6g")
        print(rr[["analysis", "n_snp", "ivw_or", "ivw_lo", "ivw_hi", "ivw_p",
                  "I2", "egger_or", "egger_p", "egger_intercept_p",
                  "wmedian_or", "wmedian_lo", "wmedian_hi", "wald_below_1"]]
              .to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    # ---------- formal comparison 1: subtype-specific (non-overlapping) ----------
    pe_only = d[(d[f"p_{PE}"] < 5e-8) & (d[f"p_{DVT}"] >= 5e-8)].copy()
    dvt_only = d[(d[f"p_{DVT}"] < 5e-8) & (d[f"p_{PE}"] >= 5e-8)].copy()
    shared = d[(d[f"p_{PE}"] < 5e-8) & (d[f"p_{DVT}"] < 5e-8)].copy()
    print(f"\nsubtype-specific leads: PE-only {len(pe_only)}, "
          f"DVT-only {len(dvt_only)}, shared {len(shared)}", flush=True)
    spec = []
    if len(pe_only) >= 3:
        spec.append(run_mr(pe_only, PE, ENDO, "PE_specific->endo"))
    if len(dvt_only) >= 3:
        spec.append(run_mr(dvt_only, DVT, ENDO, "DVT_specific->endo"))
    if len(spec) == 2:
        a, b = spec
        diff = a["ivw_beta"] - b["ivw_beta"]
        se_diff = sqrt(a["ivw_se"] ** 2 + b["ivw_se"] ** 2)
        pdiff = 2 * stats.norm.sf(abs(diff / se_diff))
        spec.append(dict(analysis="PE_specific minus DVT_specific",
                         n_snp=len(pe_only) + len(dvt_only),
                         ivw_beta=diff, ivw_se=se_diff, ivw_p=pdiff,
                         ivw_or=np.exp(diff),
                         ivw_lo=np.exp(diff - 1.96 * se_diff),
                         ivw_hi=np.exp(diff + 1.96 * se_diff)))
        sr = pd.DataFrame(spec)
        sr.to_csv(os.path.join(OUT, "severity_specific.tsv"), sep="\t",
                  index=False, float_format="%.6g")
        print(sr[["analysis", "n_snp", "ivw_or", "ivw_lo", "ivw_hi",
                  "ivw_p", "I2"]].to_string(
                      index=False, float_format=lambda v: f"{v:.3g}"))

    # ---------- formal comparison 2: MVMR PE + DVT -> endo ----------
    union = d[(d[f"p_{PE}"] < 5e-8) | (d[f"p_{DVT}"] < 5e-8)].copy()
    union = union.drop_duplicates(["chrom", "pos"])
    print(f"\nunion instruments for MVMR: {len(union)}", flush=True)
    if len(union) > 5:
        mm = mvmr(union, [PE, DVT], ENDO)
        mrows = []
        for k, v in mm.items():
            mrows.append(dict(term=k, beta=v["beta"], se=v["se"], p=v["p"],
                              OR=v["or_"], lo=v["lo"], hi=v["hi"],
                              cond_F=v.get("cond_F", np.nan)))
        mdf = pd.DataFrame(mrows)
        mdf.to_csv(os.path.join(OUT, "severity_mvmr.tsv"), sep="\t",
                   index=False, float_format="%.6g")
        print(mdf.to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    # ---------- secondary: split the 44 VTE instruments ----------
    mx = pd.read_csv(os.path.join(ROOT, "out", "m5", "instrument_matrix.tsv"),
                     sep="\t")
    mx = mx[mx["group"] == "VTE"].copy()
    mx["pe_lead"] = mx["p_I9_PULMEMB"] < 5e-8
    mx["dvt_lead"] = mx["p_I9_PHLETHROMBDVTLOW"] < 5e-8
    print("\nVTE instruments by subtype membership:")
    print(pd.crosstab(mx["pe_lead"], mx["dvt_lead"]).to_string())
    srows = []
    for lab, q in (("VTE_pe_lead", mx["pe_lead"] & ~mx["dvt_lead"]),
                   ("VTE_dvt_lead", mx["dvt_lead"] & ~mx["pe_lead"]),
                   ("VTE_both", mx["pe_lead"] & mx["dvt_lead"]),
                   ("VTE_neither", ~mx["pe_lead"] & ~mx["dvt_lead"])):
        g = mx[q]
        if len(g) >= 3:
            srows.append(run_mr(g, VTE, ENDO, lab))
    if srows:
        sr = pd.DataFrame(srows)
        sr.to_csv(os.path.join(OUT, "severity_vte_split.tsv"), sep="\t",
                  index=False, float_format="%.6g")
        print(sr[["analysis", "n_snp", "ivw_or", "ivw_lo", "ivw_hi", "ivw_p",
                  "I2"]].to_string(index=False,
                                   float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
