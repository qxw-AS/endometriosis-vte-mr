#!/usr/bin/env python3
"""P1-7: reverse MVMR of VTE on endometriosis, adjusting for BMI and smoking.

Instruments are the 44 harmonised VTE tools. Effects on VTE, BMI, ever
smoking and endometriosis are taken from the aligned FinnGen R13 matrix
built for M5. Primary model adjusts for both covariates; sensitivity
excludes the 9 coagulation-gene windows.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m5")

d = pd.read_csv(os.path.join(OUT, "instrument_matrix.tsv"), sep="\t")
pheno = {"VTE": "I9_VTE", "BMI": "BMI_IRN", "SMOKE": "EVER_SMOKER",
         "ENDO": "N14_ENDOMETRIOSIS"}

sent = pd.read_csv(os.path.join(ROOT, "out", "fig3b_sentinel_or.tsv"), sep="\t")
coag = pd.Series(False, index=d.index)
chr_s = d["vid"].str.split("-").str[0].astype(int)
pos_s = d["vid"].str.split("-").str[1].astype(int)
for _, s in sent.iterrows():
    coag |= (chr_s == int(s["chrom"])) & \
        pos_s.between(int(s["pos"]) - 1_000_000, int(s["pos"]) + 1_000_000)

base = d[d["group"] == "VTE"].copy()
print("VTE instruments", len(base), "coag-window", int(coag[base.index].sum()))


def fit(df, exposures, label, outcome="ENDO"):
    bs = [f"beta_{pheno[e]}" for e in exposures]
    ss = [f"se_{pheno[e]}" for e in exposures]
    X = df[bs].values.astype(float)
    y = df[f"beta_{pheno[outcome]}"].values.astype(float)
    sy = df[f"se_{pheno[outcome]}"].values.astype(float)
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

    Xe = np.column_stack([np.ones(len(y)), X])
    cov_e = np.linalg.inv(Xe.T @ (Xe * w[:, None]))
    coef_e = cov_e @ (Xe.T @ (w * y))
    re = y - Xe @ coef_e
    s2e = max(1.0, float(np.sum(w * re ** 2)) / (len(y) - len(exposures) - 1))
    se_e = np.sqrt(np.diag(cov_e) * s2e)
    egger_int_p = 2 * stats.t.sf(abs(coef_e[0] / se_e[0]),
                                 max(1, len(y) - len(exposures) - 1))

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
                        Q=q, Q_df=dfree, I2=i2, cond_F=cond_f[e],
                        egger_intercept=coef_e[0], egger_int_p=egger_int_p))
    return out


def ivw(df, label):
    bx = df["beta_I9_VTE"].values.astype(float)
    by = df["beta_N14_ENDOMETRIOSIS"].values.astype(float)
    sy = df["se_N14_ENDOMETRIOSIS"].values.astype(float)
    w = 1 / sy ** 2
    b = np.sum(w * bx * by) / np.sum(w * bx ** 2)
    se = np.sqrt(1 / np.sum(w * bx ** 2))
    q = float(np.sum(w * (by - b * bx) ** 2))
    dfree = len(bx) - 1
    return dict(analysis=label, exposure="VTE", n_snp=len(bx), beta=b, se=se,
                OR=np.exp(b), lo=np.exp(b - 1.96 * se),
                hi=np.exp(b + 1.96 * se),
                p=2 * stats.norm.sf(abs(b / se)), Q=q, Q_df=dfree,
                I2=max(0.0, (q - dfree) / q) if q > 0 else 0.0,
                cond_F=np.mean(bx ** 2 / sy ** 2))


rows = []
for mask_name, sub in (("all", base), ("coag_excluded", base[~coag[base.index]])):
    rows.append(ivw(sub, f"univariable_VTE_to_ENDO_{mask_name}"))
    rows.extend(fit(sub, ("VTE", "BMI"), f"MVMR_VTE_BMI_{mask_name}"))
    rows.extend(fit(sub, ("VTE", "SMOKE"), f"MVMR_VTE_SMOKE_{mask_name}"))
    rows.extend(fit(sub, ("VTE", "BMI", "SMOKE"),
                    f"MVMR_VTE_BMI_SMOKE_{mask_name}"))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "reverse_mvmr_results.tsv"), sep="\t",
           index=False, float_format="%.6g")

pd.set_option("display.width", 260)
print("\n=== reverse MVMR ===")
print(res[["analysis", "exposure", "n_snp", "OR", "lo", "hi", "p",
           "I2", "cond_F", "egger_int_p"]]
      .to_string(index=False, float_format=lambda x: f"{x:.3g}"))
