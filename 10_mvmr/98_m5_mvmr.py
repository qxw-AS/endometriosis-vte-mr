#!/usr/bin/env python3
"""M5: multivariable MR of endometriosis, BMI and smoking on VTE.

Design
- Instruments: union of lead variants for the exposures in each model
  (FinnGen R13 pheweb peaks for BMI_IRN and EVER_SMOKER; M4-clumped
  endometriosis leads). All variants were fetched in one allele frame, so
  beta/se are aligned across phenotypes.
- Primary models exclude the protocol-defined 9 coagulation genes +/-1 Mb
  (same mask as M4), because the F5 proxy otherwise dominates the
  endo->VTE estimate.
- Estimator: MVMR-IVW with uncorrelated instruments, plus MVMR-Egger
  intercept, Cochran Q, I2, and an approximate conditional F-statistic
  (residual exposure variance / exposure SE^2).
- One-sample caveat: all traits are FinnGen R13, so overlap is maximal.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m5")
os.makedirs(OUT, exist_ok=True)

d = pd.read_csv(os.path.join(OUT, "instrument_matrix.tsv"), sep="\t")
group = {"ENDO": "MVMR_ENDO", "BMI": "MVMR_BMI", "SMOKE": "MVMR_SMOKE"}
pheno = {"ENDO": "N14_ENDOMETRIOSIS", "BMI": "BMI_IRN",
         "SMOKE": "EVER_SMOKER", "VTE": "I9_VTE"}

sent = pd.read_csv(os.path.join(ROOT, "out", "fig3b_sentinel_or.tsv"), sep="\t")
coag = pd.Series(False, index=d.index)
chr_s = d["vid"].str.split("-").str[0].astype(int)
pos_s = d["vid"].str.split("-").str[1].astype(int)
for _, s in sent.iterrows():
    coag |= (chr_s == int(s["chrom"])) & \
        pos_s.between(int(s["pos"]) - 1_000_000, int(s["pos"]) + 1_000_000)

print("total instruments", len(d), "coagulation-window", int(coag.sum()))


def fit(df, exposures, label):
    bs = [f"beta_{pheno[e]}" for e in exposures]
    ss = [f"se_{pheno[e]}" for e in exposures]
    X = df[bs].values.astype(float)
    y = df[f"beta_{pheno['VTE']}"].values.astype(float)
    sy = df[f"se_{pheno['VTE']}"].values.astype(float)
    w = 1 / sy ** 2
    XtWX = X.T @ (X * w[:, None])
    beta = np.linalg.solve(XtWX, X.T @ (w * y))
    resid = y - X @ beta
    rss = float(np.sum(w * resid ** 2))
    dfree = len(y) - len(exposures)
    sigma2 = max(1.0, rss / dfree)
    cov = np.linalg.inv(XtWX) * sigma2
    se = np.sqrt(np.diag(cov))
    q = rss
    i2 = max(0.0, (q - dfree) / q) if q > 0 else 0.0

    # MVMR-Egger: weighted regression with intercept
    Xe = np.column_stack([np.ones(len(y)), X])
    cov_e = np.linalg.inv(Xe.T @ (Xe * w[:, None]))
    coef_e = cov_e @ (Xe.T @ (w * y))
    re = y - Xe @ coef_e
    s2e = max(1.0, float(np.sum(w * re ** 2)) / (len(y) - len(exposures) - 1))
    se_e = np.sqrt(np.diag(cov_e) * s2e)
    egger_int_p = 2 * stats.t.sf(abs(coef_e[0] / se_e[0]),
                                 max(1, len(y) - len(exposures) - 1))

    # approximate conditional F for each exposure
    cond_f = {}
    for k, e in enumerate(exposures):
        others = [j for j in range(len(exposures)) if j != k]
        Xo = np.column_stack([np.ones(len(y)), X[:, others]])
        coef = np.linalg.solve(Xo.T @ (Xo * w[:, None]),
                               Xo.T @ (w * X[:, k]))
        r = X[:, k] - Xo @ coef
        cond_f[e] = float(np.mean(r ** 2 / df[ss[k]].values.astype(float) ** 2))

    out = []
    for k, e in enumerate(exposures):
        out.append(dict(analysis=label, exposure=e, n_snp=len(y),
                        beta=beta[k], se=se[k], OR=np.exp(beta[k]),
                        lo=np.exp(beta[k] - 1.96 * se[k]),
                        hi=np.exp(beta[k] + 1.96 * se[k]),
                        p=2 * stats.norm.sf(abs(beta[k] / se[k])),
                        Q=q, Q_df=dfree, I2=i2,
                        cond_F=cond_f[e],
                        egger_intercept=coef_e[0], egger_int_p=egger_int_p))
    return out


def ivw(df, e, label):
    bx = df[f"beta_{pheno[e]}"].values.astype(float)
    by = df[f"beta_{pheno['VTE']}"].values.astype(float)
    sy = df[f"se_{pheno['VTE']}"].values.astype(float)
    w = 1 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se = np.sqrt(1 / np.sum(w * bx ** 2))
    q = float(np.sum(w * (by - b * bx) ** 2))
    dfree = len(bx) - 1
    return dict(analysis=label, exposure=e, n_snp=len(bx), beta=b, se=se,
                OR=np.exp(b), lo=np.exp(b - 1.96 * se),
                hi=np.exp(b + 1.96 * se),
                p=2 * stats.norm.sf(abs(b / se)), Q=q, Q_df=dfree,
                I2=max(0.0, (q - dfree) / q) if q > 0 else 0.0,
                cond_F=np.mean(bx ** 2 / df[f"se_{pheno[e]}"].values ** 2))


rows = []
for mask_name, mask in (("coag_excluded", ~coag), ("all", pd.Series(True, index=d.index))):
    sub = d[mask].copy()
    # univariable references on exposure-specific instruments
    for e in ("ENDO", "BMI", "SMOKE"):
        g = sub[sub["group"] == group[e]]
        if len(g) >= 3:
            rows.append(ivw(g, e, f"univariable_{e}_coag_{mask_name}"))
    # MVMR models
    for exps, tag in ((("ENDO", "BMI"), "MVMR_endo_BMI"),
                      (("ENDO", "SMOKE"), "MVMR_endo_smoke"),
                      (("ENDO", "BMI", "SMOKE"), "MVMR_endo_BMI_smoke")):
        use = sub[sub["group"].isin([group[e] for e in exps])]
        rows.extend(fit(use, exps, f"{tag}_{mask_name}"))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "mvmr_results.tsv"), sep="\t", index=False,
           float_format="%.6g")
print(res[["analysis", "exposure", "n_snp", "OR", "lo", "hi", "p",
           "Q", "Q_df", "I2", "cond_F", "egger_int_p"]]
      .to_string(index=False, float_format=lambda v: f"{v:.3g}"))

# compact table for the manuscript: one row per model, OR per exposure
comp = []
for (an, exp), g in res.groupby(["analysis", "exposure"]):
    pass
piv = res[res["analysis"].str.startswith("MVMR")].copy()
piv["cell"] = piv.apply(lambda r: f"{r['OR']:.2f} ({r['lo']:.2f}-{r['hi']:.2f}; P={r['p']:.3g})", axis=1)
piv = piv.pivot_table(index="analysis", columns="exposure", values="cell",
                      aggfunc="first").reset_index()
piv.to_csv(os.path.join(OUT, "mvmr_table.tsv"), sep="\t", index=False)
print("\n" + piv.to_string(index=False))
