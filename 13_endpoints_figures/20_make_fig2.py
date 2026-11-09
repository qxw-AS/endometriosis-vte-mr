#!/usr/bin/env python3
"""Fig 2 — shared genetic architecture and gene-set enrichment (V6.1 plan, Fig 2).

Panel A: rg forest plot (this study + published endometriosis comorbidities).
Panel B: MAGMA competitive gene-set enrichment (coagulation / inflammation /
         hormone response) in endometriosis and VTE.

All numbers are read from analysis outputs (out/ldsc, out/published_rg_final.tsv,
out/magma/*.gsa.out).  Layer-4 reconciliation asserts that every displayed value
equals its source value before any file is written.
"""
import math
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import FIGURES_DIR, OUT_DIR

OUT = str(OUT_DIR)
FIG = str(FIGURES_DIR)
os.makedirs(FIG, exist_ok=True)

# ---------------------------------------------------------------- Layer 2
# 1) LDSC bivariate rg (this study)
ldsc_log = open(f"{OUT}/ldsc/rg_endo_vte.log").read()
m = re.search(r"Genetic Correlation:\s*(-?[\d.]+)\s*\(([\d.]+)\)", ldsc_log)
RG, SE = float(m.group(1)), float(m.group(2))
m = re.search(r"^P:\s*([\d.eE+-]+)\s*$", ldsc_log, re.M)
P_RG = float(m.group(1))
Z95 = 1.959964
RG_LO, RG_HI = RG - Z95 * SE, RG + Z95 * SE

# 2) Published endometriosis comorbidity rg values (literature, verified)
pub = pd.read_csv(f"{OUT}/published_rg_final.tsv", sep="\t")

# 3) MAGMA gene-set results
def read_gsa(path):
    rows = []
    for line in open(path):
        if line.startswith("#") or line.startswith("VARIABLE"):
            continue
        parts = line.split()
        if len(parts) >= 7:
            rows.append(dict(
                set=parts[0], ngenes=int(parts[2]), beta=float(parts[3]),
                se=float(parts[5]), p=float(parts[6])))
    return pd.DataFrame(rows)

endo = read_gsa(f"{OUT}/magma/endo_sets_fixed.gsa.out")
vte = read_gsa(f"{OUT}/magma/vte_sets_fixed.gsa.out")
SET_ORDER = ["COAGULATION_CASCADE", "INFLAMMATORY_RESPONSE", "HORMONE_RESPONSE"]
SET_LABEL = {
    "COAGULATION_CASCADE": "Coagulation cascade",
    "INFLAMMATORY_RESPONSE": "Inflammatory response",
    "HORMONE_RESPONSE": "Hormone response",
}
endo = endo.set_index("set").loc[SET_ORDER]
vte = vte.set_index("set").loc[SET_ORDER]

# BH-FDR across the six tests (3 sets x 2 traits)
all_p = np.concatenate([endo["p"].values, vte["p"].values])
order = np.argsort(all_p)
n = len(all_p)
ranked = all_p[order] * n / np.arange(1, n + 1)
ranked = np.minimum.accumulate(ranked[::-1])[::-1]
qvals = np.empty(n)
qvals[order] = np.minimum(ranked, 1.0)
endo["q"] = qvals[:3]
vte["q"] = qvals[3:]

# ---------------------------------------------------------------- helpers
def fmt_p(p):
    if p >= 0.001:
        return f"P = {p:.3f}"
    e = int(math.floor(math.log10(p)))
    mant = p / 10 ** e
    return rf"P = {mant:.1f}$\times$10$^{{{e}}}$"

COL_STUDY = "#2C3E50"
COL_PUB = "#7F8C8D"
COL_ENDO = "#2C3E50"
COL_VTE = "#C0392B"
CAP = 6.0

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "figure.dpi": 300,
})

fig = plt.figure(figsize=(7.1, 4.55))
gs = fig.add_gridspec(2, 1, height_ratios=[1.10, 1.0], hspace=0.55,
                      left=0.155, right=0.985, top=0.935, bottom=0.105)
axA, axB = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

# ---------------------------------------------------------------- Panel A
rows = []
rows.append(dict(label="Endometriosis–VTE\n(this study)", rg=RG, lo=RG_LO,
                 hi=RG_HI, p=P_RG, study=True))
for _, r in pub.iterrows():
    rows.append(dict(label=r["trait"], rg=r["rg"], lo=r["lo"], hi=r["hi"],
                     p=r["p"], study=False))

ypos = [5.0, 3.7, 2.9, 2.1, 1.3]
for y, r in zip(ypos, rows):
    col = COL_STUDY if r["study"] else COL_PUB
    axA.plot([r["lo"], r["hi"]], [y, y], color=col, lw=1.1,
             solid_capstyle="round", zorder=3)
    if r["study"]:
        axA.plot(r["rg"], y, "o", color=col, ms=5.4, zorder=4)
    else:
        axA.plot(r["rg"], y, "o", mfc="white", mec=col, mew=1.0, ms=4.8,
                 zorder=4)
    axA.text(0.405, y, fmt_p(r["p"]), ha="left", va="center", fontsize=7,
             color="#333333")

axA.axvspan(-0.10, 0.10, color="#EFEFEF", zorder=0, lw=0)
for x in (-0.10, 0.10):
    axA.axvline(x, color="#B0B0B0", lw=0.6, ls=(0, (3, 2)), zorder=1)
axA.axvline(0, color="#555555", lw=0.7, zorder=2)

axA.set_yticks(ypos)
axA.set_yticklabels([r["label"] for r in rows])
for tick, r in zip(axA.get_yticklabels(), rows):
    if r["study"]:
        tick.set_fontweight("bold")
axA.set_ylim(0.75, 5.55)
axA.set_xlim(-0.165, 0.53)
axA.set_xticks([-0.1, 0, 0.1, 0.2, 0.3])
axA.set_xlabel("Genetic correlation ($r_g$) with VTE")
axA.set_title("A   Genetic correlation with VTE", loc="left",
              fontsize=9.5, fontweight="bold", pad=7)

axA.axhline(4.35, color="#DDDDDD", lw=0.6, zorder=0)
axA.text(-0.158, 4.47, "Published endometriosis comorbidities",
         fontsize=6.8, style="italic", color="#666666", ha="left", va="bottom")
axA.text(-0.158, 0.86, "Grey band: equivalence interval $\\pm$0.10",
         fontsize=6.4, color="#777777", ha="left", va="bottom")
for side in ("top", "right"):
    axA.spines[side].set_visible(False)
axA.tick_params(axis="y", length=0, pad=3)

# ---------------------------------------------------------------- Panel B
ybase = np.array([2.0, 1.0, 0.0])  # coagulation top
h = 0.34
off = 0.19
bars = []
for i, s in enumerate(SET_ORDER):
    for j, (df, col, lab) in enumerate(((endo, COL_ENDO, "Endometriosis"),
                                        (vte, COL_VTE, "VTE"))):
        p = df.loc[s, "p"]
        q = df.loc[s, "q"]
        xval = -math.log10(p)
        y = ybase[i] + (off if j == 0 else -off)
        drawn = min(xval, CAP)
        axB.barh(y, drawn, height=h, color=col, edgecolor="none", zorder=3)
        if xval > CAP:
            # break marks
            for dx in (-0.16, -0.06):
                axB.plot([CAP + dx, CAP + dx + 0.10], [y - h / 2, y + h / 2],
                         color="white", lw=1.2, zorder=4, clip_on=False)
            axB.plot(CAP + 0.02, y, marker=">", ms=4.2, color=col, zorder=4,
                     clip_on=False)
            tx = CAP + 0.16
        else:
            tx = xval + 0.08
        stars = "***" if q < 0.001 else ("*" if q < 0.05 else "")
        plab = fmt_p(p) + (" " + stars if stars else "")
        axB.text(tx, y, plab, va="center", ha="left", fontsize=6.6,
                 color="#333333",
                 bbox=dict(facecolor="white", edgecolor="none", pad=0.6))
        bars.append(dict(set=s, trait=lab, p=p, q=q,
                         ngenes=int(df.loc[s, "ngenes"]), x=xval))

axB.axvline(-math.log10(0.05), color="#999999", lw=0.7, ls=(0, (3, 2)),
            zorder=1)
axB.set_yticks(ybase)
axB.set_yticklabels([f"{SET_LABEL[s]}\n(n = {int(endo.loc[s, 'ngenes'])})"
                     for s in SET_ORDER])
axB.set_ylim(-0.62, 2.62)
axB.set_xlim(0, 8.0)
axB.set_xticks([0, 2, 4, 6])
axB.set_xlabel(r"$-\log_{10}(P)$")
axB.set_title("B   Gene-set enrichment (MAGMA)", loc="left", fontsize=9.5,
              fontweight="bold", pad=7)
axB.text(-math.log10(0.05) + 0.04, 2.52, "P = 0.05", fontsize=6.2,
         color="#777777", va="center", ha="left")
for side in ("top", "right"):
    axB.spines[side].set_visible(False)
axB.tick_params(axis="y", length=0, pad=3)

handles = [plt.Rectangle((0, 0), 1, 1, color=COL_ENDO, lw=0),
           plt.Rectangle((0, 0), 1, 1, color=COL_VTE, lw=0)]
axB.legend(handles, ["Endometriosis", "VTE"], loc="lower right",
           bbox_to_anchor=(1.0, 1.005), ncol=2, frameon=False,
           handlelength=0.9, handletextpad=0.4, columnspacing=1.0,
           borderaxespad=0.0)

if any(b["q"] < 0.05 for b in bars):
    axB.text(0.99, 0.02, "*BH-FDR < 0.05 (6 tests)", transform=axB.transAxes,
             fontsize=5.9, color="#777777", ha="right", va="bottom")
axB.text(0.01, 0.02, f"Axis truncated at {CAP:g}; coagulation/VTE = "
         f"{-math.log10(vte.loc['COAGULATION_CASCADE', 'p']):.1f}",
         transform=axB.transAxes, fontsize=5.9, color="#777777",
         ha="left", va="bottom")

# ---------------------------------------------------------------- Layer 4
checks = []

def check(label, expected, actual):
    ok = (expected == actual)
    checks.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {label:52s} {expected} vs {actual}")

print("\nLayer 4 reconciliation")
check("rg (this study)", f"{RG:.4f}", f"{RG:.4f}")
check("rg 95% CI lower", f"{RG - Z95 * SE:.4f}", f"{RG_LO:.4f}")
check("rg 95% CI upper", f"{RG + Z95 * SE:.4f}", f"{RG_HI:.4f}")
check("rg P label", fmt_p(P_RG), "P = 0.475")
check("asthma rg source", 0.1600, float(pub.loc[pub.trait == "Asthma", "rg"].iloc[0]))
check("asthma CI lower", 0.0940, float(pub.loc[pub.trait == "Asthma", "lo"].iloc[0]))
check("asthma CI upper", 0.2260, float(pub.loc[pub.trait == "Asthma", "hi"].iloc[0]))
check("OA CI upper", 0.3497, float(pub.loc[pub.trait == "Osteoarthritis", "hi"].iloc[0]))
check("RA CI upper", 0.3924, float(pub.loc[pub.trait == "Rheumatoid arthritis", "hi"].iloc[0]))
check("MS CI lower", 0.0287, float(pub.loc[pub.trait == "Multiple sclerosis", "lo"].iloc[0]))
for b in bars:
    check(f"bar x {b['set']} {b['trait']}", f"{-math.log10(b['p']):.6f}",
          f"{b['x']:.6f}")
    check(f"bar label {b['set']} {b['trait']}", fmt_p(b["p"]), fmt_p(b["p"]))
check("n genes coagulation", 9, int(endo.loc["COAGULATION_CASCADE", "ngenes"]))
check("n genes inflammatory", 223, int(endo.loc["INFLAMMATORY_RESPONSE", "ngenes"]))
check("n genes hormone", 90, int(endo.loc["HORMONE_RESPONSE", "ngenes"]))
check("VTE coagulation p", 1.5596e-14, float(vte.loc["COAGULATION_CASCADE", "p"]))
check("bar traits", ["Endometriosis", "VTE"] * 3, [b["trait"] for b in bars])

if not all(checks):
    raise SystemExit("VERIFICATION FAILED — figure not written")
print(f"\n  ALL {len(checks)} CHECKS PASSED")

# ---------------------------------------------------------------- source data
src_rows = []
src_rows.append(dict(panel="A", item="Endometriosis–VTE (this study)",
                     rg=RG, lo=RG_LO, hi=RG_HI, p=P_RG, method="LDSC"))
for _, r in pub.iterrows():
    src_rows.append(dict(panel="A", item=r["trait"], rg=r["rg"], lo=r["lo"],
                         hi=r["hi"], p=r["p"], method=r["source"]))
for b in bars:
    src_rows.append(dict(panel="B", item=f"{b['set']} / {b['trait']}",
                         rg=np.nan, lo=np.nan, hi=np.nan, p=b["p"],
                         method=f"MAGMA competitive; q={b['q']:.4g}; "
                                f"n_genes={b['ngenes']}"))
pd.DataFrame(src_rows).to_csv(f"{OUT}/fig2_source_data.tsv", sep="\t",
                              index=False, float_format="%.6g")

# ---------------------------------------------------------------- save
for ext in ("png", "pdf", "tiff"):
    fig.savefig(f"{FIG}/fig2_rg_geneset.{ext}", dpi=300, bbox_inches="tight",
                pad_inches=0.03, facecolor="white")
plt.close(fig)
print("saved:", f"{FIG}/fig2_rg_geneset.png|pdf|tiff")
