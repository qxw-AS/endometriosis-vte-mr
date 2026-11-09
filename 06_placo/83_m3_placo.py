#!/usr/bin/env python3
"""M3 stage 4: PLACO (Ray & Chatterjee 2020) composite-null pleiotropy test.

Implementation follows the published approximate p-value (S1 Appendix Eq 3):

  p = F(t / sqrt(VarZ1)) + F(t / sqrt(VarZ2)) - F(t),   t = z1*z2

where F is the two-sided tail of the normal-product distribution,
F(u) = (2/pi) * integral_{|u|}^inf K0(w) dw, and VarZ1/VarZ2 are marginal
variances of the single-trait Z-scores under the composite null (estimated
from variants not significant for both traits at p < 1e-4). Correlation
between Z-scores from shared FinnGen controls is estimated from null variants
and removed by R^-1/2 before testing (as in PLACO).

Both FinnGen files are position-sorted but contain slightly different variant
sets, so the two files are merge-joined with a streaming buffer instead of
assuming row alignment.
"""
import os

import numpy as np
import pandas as pd
from scipy import special

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "out", "m3")

FILES = {"endo": os.path.join(DATA, "finngen_R13_N14_ENDOMETRIOSIS.gz"),
         "vte": os.path.join(DATA, "finngen_R13_I9_VTE.gz")}
COLS = ["#chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"]
MIN_MAF = 0.01
CHUNK = 500_000
KEY = ["chr", "pos", "ref", "alt"]

_CHR_RANK = {str(i): i for i in range(1, 23)}
_CHR_RANK.update({"X": 23, "Y": 24, "MT": 25, "M": 25})


def rank(chrom):
    return _CHR_RANK.get(str(chrom), 26)


def prep(chunk, detail):
    d = chunk.rename(columns={"#chrom": "chr", "pval": "p", "sebeta": "se",
                              "af_alt": "af"})
    d = d[(d["se"] > 0) & (d["af"] > 0) & (d["af"] < 1)].copy()
    maf = np.minimum(d["af"], 1 - d["af"])
    d = d[maf >= MIN_MAF].copy()
    keep = ["chr", "pos", "ref", "alt", "p", "beta", "se"]
    if detail:
        d["rsid"] = d["rsids"].astype(str).str.split(",").str[0]
        keep += ["rsid"]
    return d[keep]


def iter_merged(detail=False):
    it1 = pd.read_csv(FILES["endo"], sep="\t", chunksize=CHUNK, usecols=COLS,
                      dtype={"#chrom": str, "rsids": str}, low_memory=False)
    it2 = pd.read_csv(FILES["vte"], sep="\t", chunksize=CHUNK, usecols=COLS,
                      dtype={"#chrom": str, "rsids": str}, low_memory=False)
    buf, buf_max = None, None
    for c1 in it1:
        a = prep(c1, detail)
        if a.empty:
            continue
        a_min = (rank(a["chr"].iloc[0]), int(a["pos"].iloc[0]))
        a_max = (rank(a["chr"].iloc[-1]), int(a["pos"].iloc[-1]))
        while buf is None or buf_max < a_max:
            c2 = next(it2, None)
            if c2 is None:
                break
            b = prep(c2, detail)
            if b.empty:
                continue
            buf = b if buf is None else pd.concat([buf, b], ignore_index=True)
            buf_max = (rank(b["chr"].iloc[-1]), int(b["pos"].iloc[-1]))
        if buf is None:
            break
        m = a.merge(buf, on=KEY, suffixes=("_e", "_v"), how="inner")
        if detail and "rsid_e" in m.columns:
            m = m.rename(columns={"rsid_e": "rsid"})
        yield m
        # discard buffer rows before the next possible endo chunk start
        brank = buf["chr"].map(rank).astype(np.int64)
        keep = (brank * 10**12 + buf["pos"].astype(np.int64)) >= \
               (a_min[0] * 10**12 + a_min[1])
        buf = buf[keep]
        if buf.empty:
            buf, buf_max = None, None


# ---------- normal-product two-sided tail on a grid ----------
_u = np.linspace(0, 80, 80_001)
_k0 = special.k0(_u)
_k0[0] = special.k0(1e-6)                          # K0(0) = inf; regularise
_cum = np.concatenate([[0.0], np.cumsum(
    (_k0[1:] + _k0[:-1]) / 2 * np.diff(_u))])      # integral_0^{u_i} K0
_total = _cum[-1] + special.k0(80.0)               # integral_0^inf K0
_tail = np.maximum(_total - _cum, 0.0)             # integral_{u_i}^inf K0
_logF = np.log(2 / np.pi) + np.log(np.maximum(_tail, 1e-300))


def F(u):
    a = np.abs(np.asarray(u, dtype=np.float64))
    out = np.empty_like(a)
    # interpolation is accurate while the cumulative tail is well above the
    # ~1e-13 cancellation floor of the cumulative integral; beyond that use
    # the asymptotic integral_{u}^inf K0 ~ K0(u)
    inr = a <= 15
    out[inr] = np.exp(np.interp(a[inr], _u, _logF))
    big = a[~inr]
    out[~inr] = np.exp(np.log(2 / np.pi) + np.log(special.k0e(big)) - big)
    return np.clip(out, 0, 1)


_rng = np.random.default_rng(20261108)
_z1, _z2 = _rng.standard_normal(400_000), _rng.standard_normal(400_000)
_p = F(_z1 * _z2)
print(f"[self-test] global null: P(p<0.05)={np.mean(_p < 0.05):.4f} "
      f"(expected 0.050), P(p<1e-4)={np.mean(_p < 1e-4):.2e} (expected 1e-4)")

# ---------- pass 1 ----------
acc = {k: [] for k in ("chr", "pos", "z1", "z2", "p1", "p2")}
for m in iter_merged(detail=False):
    if m.empty:
        continue
    acc["chr"].append(m["chr"].values.astype("int16"))
    acc["pos"].append(m["pos"].values.astype("int32"))
    acc["z1"].append((m["beta_e"] / m["se_e"]).values.astype("float32"))
    acc["z2"].append((m["beta_v"] / m["se_v"]).values.astype("float32"))
    acc["p1"].append(m["p_e"].values.astype("float32"))
    acc["p2"].append(m["p_v"].values.astype("float32"))

chrom = np.concatenate(acc["chr"])
pos = np.concatenate(acc["pos"])
z1 = np.concatenate(acc["z1"]).astype(np.float64)
z2 = np.concatenate(acc["z2"]).astype(np.float64)
p1 = np.concatenate(acc["p1"]).astype(np.float64)
p2 = np.concatenate(acc["p2"]).astype(np.float64)
del acc
print(f"common-variant intersection: {len(z1):,}")

null = (p1 > 1e-4) & (p2 > 1e-4)
rho = np.corrcoef(z1[null], z2[null])[0, 1]
print(f"Z correlation from null variants (n={int(null.sum()):,}): rho={rho:.4f}")

d1, d2 = 1 / np.sqrt(1 + rho), 1 / np.sqrt(1 - rho)
aa, bb = (d1 + d2) / 2, (d1 - d2) / 2
u1, u2 = aa * z1 + bb * z2, bb * z1 + aa * z2
print(f"post-decorrelation corr = {np.corrcoef(u1[null], u2[null])[0, 1]:.4f}")

comp_null = ~((p1 < 1e-4) & (p2 < 1e-4))
varz1 = np.var(u1[comp_null], ddof=1)
varz2 = np.var(u2[comp_null], ddof=1)
print(f"composite-null n={int(comp_null.sum()):,}; "
      f"VarZ1={varz1:.4f}, VarZ2={varz2:.4f}")

t = u1 * u2
placo = np.clip(F(t / np.sqrt(varz1)) + F(t / np.sqrt(varz2)) - F(t), 0, 1)
print(f"PLACO p: <0.05 {np.mean(placo < 0.05):.4f}; "
      f"<5e-8 {int((placo < 5e-8).sum()):,}; <1e-5 {int((placo < 1e-5).sum()):,}")

with open(os.path.join(OUT, "placo_params.txt"), "w") as fh:
    fh.write(f"n_variants\t{len(z1)}\nrho\t{rho:.6f}\n"
             f"var_z1\t{varz1:.6f}\nvar_z2\t{varz2:.6f}\n"
             f"n_null\t{int(null.sum())}\nn_composite_null\t{int(comp_null.sum())}\n"
             f"n_p_lt_5e8\t{int((placo < 5e-8).sum())}\n"
             f"n_p_lt_1e5\t{int((placo < 1e-5).sum())}\n")

np.savez_compressed(os.path.join(OUT, "placo_all_scores.npz"),
                    chrom=chrom, pos=pos, t=t.astype("float32"),
                    placo=placo.astype("float32"))

# ---------- pass 2: details for p < 1e-5 ----------
sel = placo < 1e-5
hit_keys = set(zip(chrom[sel].tolist(), pos[sel].tolist()))
rows = []
for m in iter_merged(detail=True):
    if m.empty:
        continue
    k = list(zip(m["chr"].astype(int).tolist(), m["pos"].tolist()))
    keep = np.fromiter((kk in hit_keys for kk in k), dtype=bool, count=len(k))
    if keep.any():
        rows.append(m[keep])
hits = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
if not hits.empty:
    hits["z1"] = hits["beta_e"] / hits["se_e"]
    hits["z2"] = hits["beta_v"] / hits["se_v"]
    h1 = aa * hits["z1"].values + bb * hits["z2"].values
    h2 = bb * hits["z1"].values + aa * hits["z2"].values
    tt = h1 * h2
    hits["placo_p"] = np.clip(F(tt / np.sqrt(varz1)) + F(tt / np.sqrt(varz2))
                              - F(tt), 0, 1)
    hits["T"] = tt
    hits = hits.rename(columns={"p_e": "p_endo", "p_v": "p_vte"})
    hits = hits[["chr", "pos", "rsid", "ref", "alt", "beta_e", "se_e", "p_endo",
                 "beta_v", "se_v", "p_vte", "z1", "z2", "T", "placo_p"]]
    hits = hits.sort_values("placo_p")
    hits.to_csv(os.path.join(OUT, "placo_hits_p1e5.tsv.gz"), sep="\t",
                index=False, float_format="%.6g", compression="gzip")
    print(f"recovered {len(hits):,} hit rows; p<5e-8: "
          f"{int((hits['placo_p'] < 5e-8).sum()):,}")
    print(hits.head(20).to_string(index=False, float_format=lambda v: f"{v:.3g}"))
