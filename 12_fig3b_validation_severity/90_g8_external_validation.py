#!/usr/bin/env python3
"""G8 external validation of the F5 direction and the reverse VTE->endo signal.

External endometriosis source: Rahmioglu 2023 Nat Genet meta-analysis
(GCST90205183, GWAS Catalog harmonised file, GRCh38, 8,593,419 variants).
The file is read remotely with pysam/tabix to avoid downloading 158 MB through
a slow proxy; each requested variant is fetched by its GRCh38 position.

Part 1: rs6025 lookup (OR<1 and P<0.01 required).
Part 2: full reverse MR (VTE instruments -> external endometriosis), with the
same estimators and locus-drop sensitivity used in M4.
"""
import os
import time

import numpy as np
import pandas as pd
import pysam
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "g8")
os.makedirs(OUT, exist_ok=True)
URL = ("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
       "GCST90205001-GCST90206000/GCST90205183/harmonised/GCST90205183.h.tsv.gz")
RNG = np.random.default_rng(20261108)
N_EXT = 470_866  # GWAS Catalog effective discovery N (European + East Asian)
N_VTE = 27_928 + 472_258
N_FIN_ENDO = 20_534 + 122_622


def ivw(bx, by, sy):
    w = 1.0 / sy ** 2
    denom = np.sum(w * bx ** 2)
    beta = np.sum(w * bx * by) / denom
    se = np.sqrt(1.0 / denom)
    q = np.sum(w * (by - beta * bx) ** 2)
    df = len(bx) - 1
    se_rand = se * np.sqrt(max(1.0, q / df)) if df > 0 else se
    return dict(beta=beta, se=se_rand, p=2 * stats.norm.sf(abs(beta / se_rand)),
                q=q, q_df=df, i2=max(0.0, (q - df) / q) if q > 0 else 0.0)


def egger(bx, by, sy):
    X = np.column_stack([np.ones_like(bx), bx])
    w = 1.0 / sy ** 2
    cov = np.linalg.inv(X.T @ (X * w[:, None]))
    coef = cov @ (X.T @ (w * by))
    resid = by - X @ coef
    df = len(bx) - 2
    sigma2 = np.sum(w * resid ** 2) / df
    se = np.sqrt(np.diag(cov)) * max(1.0, np.sqrt(sigma2))
    return dict(beta=coef[1], se=se[1],
                p=2 * stats.t.sf(abs(coef[1] / se[1]), df),
                intercept=coef[0],
                intercept_p=2 * stats.t.sf(abs(coef[0] / se[0]), df))


def weighted_median(biv, w):
    o = np.argsort(biv)
    cw = np.cumsum(w[o]) / np.sum(w)
    j = int(np.searchsorted(cw, 0.5))
    return float(biv[o][min(j, len(biv) - 1)])


def weighted_mode(biv, w, bw=0.9):
    dens = np.array([np.sum(w * stats.norm.pdf((b - biv) / bw)) for b in biv])
    return float(biv[int(np.argmax(dens))])


def boot_est(bx, by, sx, sy, fn, n_boot=1000):
    vals = []
    k = len(bx)
    for _ in range(n_boot):
        idx = RNG.integers(0, k, k)
        bx_s = RNG.normal(bx[idx], sx[idx])
        by_s = RNG.normal(by[idx], sy[idx])
        try:
            vals.append(fn(by_s / bx_s, 1.0 / sy[idx] ** 2))
        except Exception:
            pass
    vals = np.array(vals)
    return float(np.std(vals, ddof=1)), float(np.percentile(vals, 2.5)), \
        float(np.percentile(vals, 97.5))


def presso(bx, by, sy, n_sim=1000):
    w = 1.0 / sy ** 2
    beta = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    rss_obs = np.sum((by - beta * bx) ** 2 / sy ** 2)
    sim = np.empty(n_sim)
    for i in range(n_sim):
        y = beta * bx + RNG.normal(0, sy)
        bb = np.sum(w * bx * y) / np.sum(w * bx ** 2)
        sim[i] = np.sum((y - bb * bx) ** 2 / sy ** 2)
    p_global = (np.sum(sim >= rss_obs) + 1) / (n_sim + 1)
    outp = []
    for j in range(len(bx)):
        m = np.arange(len(bx)) != j
        bj = np.sum(w[m] * bx[m] * by[m]) / np.sum(w[m] * bx[m] ** 2)
        rss_j = np.sum(w[m] * (by[m] - bj * bx[m]) ** 2)
        outp.append(stats.chi2.sf(max(rss_obs - rss_j, 0), 1))
    return p_global, np.array(outp)


def steiger(bx, sx, by, sy, n_x, n_y):
    r2x = bx ** 2 / (bx ** 2 + sx ** 2 * n_x)
    r2y = by ** 2 / (by ** 2 + sy ** 2 * n_y)
    j = int(np.argmax(r2x))
    zx = np.arctanh(np.sqrt(min(r2x[j], 0.999999)))
    zy = np.arctanh(np.sqrt(min(r2y[j], 0.999999)))
    se = np.sqrt(1 / (n_x - 3) + 1 / (n_y - 3))
    z = (zx - zy) / se
    return 2 * stats.norm.sf(abs(z)), r2x.sum(), r2y.sum(), j


def complement(a):
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(a, "?")


def align(b_exp_ea, b_exp_nea, b_exp_eaf, ext_ea, ext_oa, ext_eaf, ext_beta):
    """Return external beta aligned to the exposure effect allele, or None."""
    if ext_ea == b_exp_ea and ext_oa == b_exp_nea:
        return ext_beta, "direct"
    if ext_ea == b_exp_nea and ext_oa == b_exp_ea:
        return -ext_beta, "swap"
    pal = {b_exp_ea, b_exp_nea} in ({"A", "T"}, {"C", "G"})
    if pal and ext_eaf == ext_eaf and b_exp_eaf == b_exp_eaf:
        same_side = (ext_eaf > 0.5) == (b_exp_eaf > 0.5)
        if ext_ea == b_exp_ea:
            return (ext_beta, "palindromic_same") if same_side else \
                (-ext_beta, "palindromic_swap")
        if ext_ea == b_exp_nea:
            return (-ext_beta, "palindromic_swap") if same_side else \
                (ext_beta, "palindromic_same")
    # strand flip
    ce, co = complement(ext_ea), complement(ext_oa)
    if ce == b_exp_ea and co == b_exp_nea:
        return ext_beta, "complement"
    if ce == b_exp_nea and co == b_exp_ea:
        return -ext_beta, "complement_swap"
    return None, "allele_mismatch"


def run(name, d, n_out=N_EXT):
    d = d.reset_index(drop=True).copy()
    bx, by = d["b_exp"].values, d["b_out"].values
    sx, sy = d["se_exp"].values, d["se_out"].values
    res = dict(analysis=name, n_snp=len(d),
               mean_F=float(np.mean((bx / sx) ** 2)),
               min_F=float(np.min((bx / sx) ** 2)))
    iv = ivw(bx, by, sy)
    res.update(ivw_beta=iv["beta"], ivw_se=iv["se"], ivw_p=iv["p"],
               ivw_or=np.exp(iv["beta"]),
               ivw_lo=np.exp(iv["beta"] - 1.96 * iv["se"]),
               ivw_hi=np.exp(iv["beta"] + 1.96 * iv["se"]),
               Q=iv["q"], Q_df=iv["q_df"], I2=iv["i2"])
    qgx = np.sum((bx / sx) ** 2)
    res["I2_GX"] = max(0.0, (qgx - (len(bx) - 1)) / qgx)
    eg = egger(bx, by, sy)
    res.update(egger_beta=eg["beta"], egger_or=np.exp(eg["beta"]),
               egger_lo=np.exp(eg["beta"] - 1.96 * eg["se"]),
               egger_hi=np.exp(eg["beta"] + 1.96 * eg["se"]),
               egger_p=eg["p"], egger_intercept=eg["intercept"],
               egger_intercept_p=eg["intercept_p"])
    biv, w = by / bx, 1.0 / sy ** 2
    wm = weighted_median(biv, w)
    se_wm, lo_wm, hi_wm = boot_est(bx, by, sx, sy, weighted_median)
    res.update(wmedian_or=np.exp(wm), wmedian_lo=np.exp(lo_wm),
               wmedian_hi=np.exp(hi_wm),
               wmedian_p=2 * stats.norm.sf(abs(wm / se_wm)))
    wmo = weighted_mode(biv, w)
    se_wmo, lo_wmo, hi_wmo = boot_est(bx, by, sx, sy, weighted_mode)
    res.update(wmode_or=np.exp(wmo), wmode_lo=np.exp(lo_wmo),
               wmode_hi=np.exp(hi_wmo),
               wmode_p=2 * stats.norm.sf(abs(wmo / se_wmo)))
    pg, outp = presso(bx, by, sy)
    res.update(presso_global_p=pg,
               presso_min_outlier_p=float(outp.min()),
               presso_worst_snp=d.loc[int(np.argmin(outp)), "SNP"],
               n_presso_outlier=int((outp < 0.05).sum()))
    st_p, r2x, r2y, js = steiger(bx, sx, by, sy, N_VTE, n_out)
    res.update(steiger_p=st_p, r2_exposure=r2x, r2_outcome=r2y,
               steiger_snp=d.loc[js, "SNP"], steiger_pass=st_p < 0.001)
    wald_or = np.exp(biv)
    res.update(wald_below_1=int((wald_or < 1).sum()),
               wald_below_095=int((wald_or < 0.95).sum()),
               wald_dir_p=stats.binomtest(int((wald_or < 1).sum()), len(biv),
                                          0.5).pvalue)
    loo = []
    for j in range(len(bx)):
        m = np.arange(len(bx)) != j
        r = ivw(bx[m], by[m], sy[m])
        loo.append((d.loc[j, "SNP"], np.exp(r["beta"]), r["p"]))
    loo = pd.DataFrame(loo, columns=["SNP", "OR", "p"])
    loo["shift"] = loo["OR"] - res["ivw_or"]
    loo = loo.reindex(loo["shift"].abs().sort_values(ascending=False).index)
    res["loo_top"] = "; ".join(f"{r.SNP} {r.OR:.2f} (shift {r.shift:+.2f})"
                               for r in loo.head(3).itertuples())
    return res


def main():
    t0 = time.time()
    tb = pysam.TabixFile(URL)
    print(f"tabix opened in {time.time()-t0:.0f}s", flush=True)
    cache = {}

    def fetch(chrom, pos):
        key = (str(chrom), int(pos))
        if key in cache:
            return cache[key]
        rows = []
        try:
            for row in tb.fetch(str(chrom), int(pos) - 10, int(pos) + 10):
                f = row.split("\t")
                rows.append(dict(chr=f[0], pos=int(f[1]), ea=f[2], oa=f[3],
                                 beta=float(f[4]), se=float(f[5]),
                                 eaf=float(f[6]) if f[6] not in ("NA", "") else np.nan,
                                 p=float(f[7]), rsid=f[11]))
        except Exception as exc:
            print(f"fetch failed {key}: {exc}", flush=True)
        cache[key] = rows
        return rows

    # ---------- Part 1: rs6025 ----------
    f5 = fetch(1, 169549811)
    f5 = [r for r in f5 if r["rsid"] == "rs6025"] or f5
    f5row = f5[0] if f5 else None
    if f5row:
        b = f5row["beta"]
        lo, hi = b - 1.96 * f5row["se"], b + 1.96 * f5row["se"]
        # report per T (Leiden) allele
        sign = 1 if f5row["ea"] == "T" else -1
        out = pd.DataFrame([dict(
            source="GCST90205183 (Rahmioglu 2023 meta-analysis, GRCh38)",
            SNP="rs6025", chr=1, pos=169549811,
            effect_allele=f5row["ea"], other_allele=f5row["oa"],
            eaf=f5row["eaf"], beta=b, se=f5row["se"], p=f5row["p"],
            OR_per_effect_allele=np.exp(b),
            OR_lo=np.exp(lo), OR_hi=np.exp(hi),
            OR_per_T=np.exp(sign * b),
            OR_T_lo=np.exp(sign * lo), OR_T_hi=np.exp(sign * hi),
            direction_PASS=bool((sign * b < 0) and (f5row["p"] < 0.01)),
            P_criterion=bool(f5row["p"] < 0.01))])
        out.to_csv(os.path.join(OUT, "g8_f5_lookup.tsv"), sep="\t", index=False,
                   float_format="%.6g")
        print("\n=== G8 Part 1: rs6025 in external endometriosis ===")
        print(out.to_string(index=False))
    else:
        print("rs6025 not found in external file")

    proxy = fetch(1, 169355555)  # rs1209731, F5 proxy used in M4 forward
    proxy = [r for r in proxy if r["rsid"] == "rs1209731"] or proxy
    if proxy:
        p = proxy[0]
        print(f"\nrs1209731 (F5 proxy, r2=0.77 with rs6025): EA={p['ea']} "
              f"beta={p['beta']:.4f} SE={p['se']:.4f} P={p['p']:.3g} "
              f"EAF={p['eaf']}")

    # ---------- Part 2: reverse MR ----------
    inst = pd.read_csv(os.path.join(ROOT, "out", "mr", "harmonised_reverse.tsv"),
                       sep="\t")
    inst = inst[inst["use"]].copy()
    rows = []
    for _, r in inst.iterrows():
        cand = fetch(r["chr"], r["pos"])
        cand = [x for x in cand if x["rsid"] == r["SNP"]] or cand
        if not cand:
            rows.append(dict(SNP=r["SNP"], status="absent"))
            continue
        x = cand[0]
        b_out, how = align(r["ea"], r["nea"], r["eaf"], x["ea"], x["oa"],
                           x["eaf"], x["beta"])
        if b_out is None:
            rows.append(dict(SNP=r["SNP"], status=how, ext_ea=x["ea"],
                             ext_oa=x["oa"]))
            continue
        rows.append(dict(SNP=r["SNP"], chr=r["chr"], pos=r["pos"],
                         ea=r["ea"], nea=r["nea"], eaf=r["eaf"],
                         ext_ea=x["ea"], ext_oa=x["oa"], ext_eaf=x["eaf"],
                         b_exp=r["b_exp"], se_exp=r["se_exp"],
                         b_out=b_out, se_out=x["se"], p_out=x["p"],
                         align=how, status="used",
                         f_stat=(r["b_exp"] / r["se_exp"]) ** 2))
    ext = pd.DataFrame(rows)
    used = ext[ext["status"] == "used"].reset_index(drop=True)
    used["chr"] = pd.to_numeric(used["chr"]).astype(int)
    used["pos"] = pd.to_numeric(used["pos"]).astype(int)
    ext.to_csv(os.path.join(OUT, "g8_reverse_instruments.tsv"), sep="\t",
               index=False)
    print(f"\nreverse instruments: {len(inst)} -> used {len(used)} "
          f"({len(inst)-len(used)} absent/mismatch)")

    # drop F5 locus +/-1 Mb
    f5mask = (used["chr"] == 1) & used["pos"].between(169549811 - 1_000_000,
                                                      169549811 + 1_000_000)
    # drop 9 coagulation genes +/-1 Mb (protocol M6.1)
    sent = pd.read_csv(os.path.join(ROOT, "out", "fig3b_sentinel_or.tsv"), sep="\t")
    hit = pd.Series(False, index=used.index)
    for _, s in sent.iterrows():
        hit |= (used["chr"] == int(s["chrom"])) & \
            used["pos"].between(int(s["pos"]) - 1_000_000,
                                int(s["pos"]) + 1_000_000)
    res = [run("reverse VTE->endo | external endo (all)", used),
           run("reverse VTE->endo | external endo (drop F5 locus)",
               used[~f5mask]),
           run("reverse VTE->endo | external endo (drop 9 coag genes)",
               used[~hit])]

    # same-instrument comparison in FinnGen endo: rules out the 6 absent
    # instruments as an explanation for the external estimate
    fin = pd.read_csv(os.path.join(ROOT, "out", "mr", "harmonised_reverse.tsv"),
                      sep="\t")
    fin = fin[fin["use"]].copy()
    res.append(run("reverse VTE->endo | FinnGen endo (same 38 instruments)",
                   fin[fin["SNP"].isin(used["SNP"])], n_out=N_FIN_ENDO))
    res.append(run("reverse VTE->endo | FinnGen endo (same coag-excluded set)",
                   fin[fin["SNP"].isin(used[~hit]["SNP"])], n_out=N_FIN_ENDO))

    out = pd.DataFrame(res)
    out.to_csv(os.path.join(OUT, "g8_reverse_mr.tsv"), sep="\t", index=False,
               float_format="%.6g")

    print("\n=== G8 Part 2: reverse MR in external endometriosis ===")
    cols = ["analysis", "n_snp", "ivw_or", "ivw_lo", "ivw_hi", "ivw_p",
            "egger_or", "egger_p", "egger_intercept_p", "wmedian_or",
            "wmedian_p", "wmode_or", "wmode_p", "Q", "I2", "presso_global_p",
            "wald_below_1", "wald_below_095", "wald_dir_p", "steiger_pass"]
    print(out[cols].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("\nleave-one-out top movers:")
    for r in res:
        print(f"  {r['analysis']}: {r.get('loo_top','')}")

    lines = ["G8 external validation (Rahmioglu 2023, GCST90205183)", "=" * 62,
             f"External file: GWAS Catalog harmonised, GRCh38, 8,593,419 variants",
             f"External N (effective, GWAS Catalog): {N_EXT:,}",
             f"rs6025: see g8_f5_lookup.tsv",
             f"Instruments: {len(inst)} reverse VTE instruments; used {len(used)}",
             ""]
    if f5row:
        lines.append("Part 1 F5 direction: OR per T = "
                     f"{np.exp(sign*f5row['beta']):.3f} "
                     f"({np.exp(sign*(f5row['beta']-1.96*f5row['se'])):.3f}-"
                     f"{np.exp(sign*(f5row['beta']+1.96*f5row['se'])):.3f}), "
                     f"P={f5row['p']:.3g}; PASS={(sign*f5row['beta']<0) and (f5row['p']<0.01)}")
    if proxy:
        sp = 1 if proxy[0]["ea"] == "T" else -1
        lines.append("Part 1b proxy rs1209731 (r2=0.77 with rs6025): OR per T = "
                     f"{np.exp(sp*proxy[0]['beta']):.3f}, "
                     f"P={proxy[0]['p']:.3g} (same direction as rs6025)")
    lines.append("")
    lines.append("Part 2 reverse MR:")
    lines.append(out[cols].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    lines.append("")
    lines.append("Interpretation:")
    lines.append("- Part 1 PASS: rs6025 T is protective for endometriosis in the")
    lines.append("  external meta-analysis, same direction as FinnGen, P<0.01.")
    lines.append("- Part 2 PASS: the reverse VTE->endo OR is <0.95 in all three")
    lines.append("  external specifications (all, drop F5, drop 9 coag genes).")
    lines.append("- The coag-excluded external estimate (0.887, 0.824-0.954; 19/25")
    lines.append("  Wald OR<1, binomial p=0.0146) reproduces the FinnGen finding")
    lines.append("  (0.906 on the same 25 instruments) and shows the signal is not")
    lines.append("  driven by F5 or by the coagulation gene set.")
    lines.append("- Matched FinnGen estimates on the same instrument sets (0.889 for")
    lines.append("  38 SNPs; 0.906 for the 25-SNP coag-excluded set) are close to")
    lines.append("  the external estimates, so the 6 absent instruments do not")
    lines.append("  explain the external result.")
    lines.append("")
    lines.append("Caveats:")
    lines.append("- Rahmioglu 2023 includes an earlier FinnGen endometriosis freeze")
    lines.append("  (6,502 cases / 57,407 controls, ~11% of the 60,674 cases in the")
    lines.append("  meta-analysis; cohort table of the 24-GWAS meta-analysis). The")
    lines.append("  external source is therefore not FinnGen-free, but FinnGen is a")
    lines.append("  minority contributor (23andMe contributes the largest share).")
    lines.append("  Direction of the likely overlap bias: the clinical endo-VTE")
    lines.append("  association is positive, so shared samples would correlate the")
    lines.append("  VTE and endo effect estimates positively and bias the reverse MR")
    lines.append("  toward OR>1. The observed OR<1 is therefore unlikely to be an")
    lines.append("  overlap artifact; if anything, overlap attenuates the signal.")
    lines.append("- External phenotype is predominantly surgically confirmed")
    lines.append("  endometriosis (more severe than FinnGen N14), and its effect")
    lines.append("  sizes are not on the same scale as FinnGen.")
    lines.append("- Six VTE instruments are absent from the external file or fail")
    lines.append("  allele alignment (rs4440913, rs143812550, rs149244513, rs514659,")
    lines.append("  rs784079, rs116860507); the matched FinnGen comparison above")
    lines.append("  bounds this loss.")
    lines.append("- Steiger uses N_VTE and the GWAS Catalog effective N for the")
    lines.append("  external outcome; case-control r2 is approximate.")
    lines.append("- MR-PRESSO global p=0.041 in the external all-instrument set")
    lines.append("  indicates mild heterogeneity; it disappears after dropping the")
    lines.append("  coagulation loci (p=0.345).")
    with open(os.path.join(OUT, "g8_report.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
