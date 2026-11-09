#!/usr/bin/env python3
"""Fig 3 (V6.2, revised): antagonistic pleiotropy at F5 rs6025 and sentinel-SNP ORs.

Panel A  forest plot of rs6025 across 4 thrombotic + 4 gynaecological FinnGen R13
         phenotypes, ordered by OR inside each spectrum (protocol Table 41).
Panel B  sentinel-SNP ORs for the 9 coagulation genes. The sentinel is the SNP with
         the lowest VTE P inside the MAGMA gene window (+/-10 kb, 1000G EUR); the
         endometriosis OR of that same SNP is shown unselected. This replaces the
         SNPwise-mean MAGMA scatter, which diluted single strong variants.
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

a = pd.read_csv(os.path.join(OUT, "fig3a_f5_spectrum.tsv"), sep="\t")
b = pd.read_csv(os.path.join(OUT, "fig3b_sentinel_or.tsv"), sep="\t")
b = b.sort_values("OR_vte", ascending=False).reset_index(drop=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7.5,
    "axes.linewidth": 0.7,
    "axes.labelsize": 8.2,
    "xtick.labelsize": 7.2,
    "ytick.labelsize": 7.2,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "dejavusans",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

COL = {"thrombotic": "#B2182B", "gynaecological": "#2166AC"}
MARK = {"thrombotic": "s", "gynaecological": "o"}
Y_POS = [8.2, 7.3, 6.4, 5.5, 3.7, 2.8, 1.9, 1.0]

fig = plt.figure(figsize=(7.1, 6.8))

# ------------------------------------------------------------------ Panel A
axA = fig.add_axes([0.205, 0.585, 0.59, 0.355])
for i, row in a.iterrows():
    y = Y_POS[i]
    col = COL[row["lineage"]]
    axA.plot([row["OR_lo"], row["OR_hi"]], [y, y], color=col, lw=1.15,
             solid_capstyle="round", zorder=3)
    axA.plot(row["OR"], y, marker=MARK[row["lineage"]], ms=4.7, color=col,
             mec="white", mew=0.6, zorder=4)
    axA.text(1.018, y, f"{row['OR']:.2f} ({row['OR_lo']:.2f}\u2013{row['OR_hi']:.2f})",
             transform=axA.get_yaxis_transform(), va="center", ha="left",
             fontsize=6.9, color="#222222")

axA.axvline(1.0, color="#555555", lw=0.8, ls=(0, (4, 3)), zorder=1)
axA.axhline(4.55, color="#D9D9D9", lw=0.6, zorder=0)
axA.set_xscale("log")
axA.set_xlim(0.62, 3.2)
axA.set_xticks([0.7, 0.8, 1.0, 1.5, 2.0, 3.0])
axA.set_xticklabels(["0.7", "0.8", "1.0", "1.5", "2.0", "3.0"])
axA.xaxis.set_minor_locator(NullLocator())
axA.set_yticks(Y_POS)
axA.set_yticklabels(a["phenotype"])
axA.set_ylim(0.55, 9.05)
axA.set_xlabel("Odds ratio per rs6025 T allele (95% CI)")
axA.text(0.628, 8.88, "Thrombotic spectrum", color=COL["thrombotic"],
         fontsize=7.2, fontstyle="italic", fontweight="bold", va="center", ha="left",
         bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
axA.text(0.628, 4.28, "Gynaecological spectrum", color=COL["gynaecological"],
         fontsize=7.2, fontstyle="italic", fontweight="bold", va="center", ha="left",
         bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
axA.text(1.018, 1.004, "OR (95% CI)", transform=axA.transAxes,
         fontsize=6.8, fontweight="bold", va="bottom", ha="left", color="#333333")
axA.set_title("A   F5 rs6025: two phenotypic spectra", loc="left",
              fontsize=9.5, fontweight="bold", pad=7)
for side in ("top", "right"):
    axA.spines[side].set_visible(False)

# ------------------------------------------------------------------ Panel B
axB = fig.add_axes([0.245, 0.085, 0.585, 0.375])
Y_B = list(range(len(b), 0, -1))
YOFF = 0.14
for i, row in b.iterrows():
    y = Y_B[i]
    opp = bool(row["opposite"])
    axB.plot([row["OR_endo"], row["OR_vte"]], [y, y],
             color="#E69F00" if opp else "#CCCCCC", lw=1.0,
             alpha=0.95 if opp else 0.8, zorder=2)
    axB.plot([row["OR_endo_lo"], row["OR_endo_hi"]], [y - YOFF, y - YOFF],
             color="#2166AC", lw=1.05, solid_capstyle="round", zorder=3)
    axB.plot(row["OR_endo"], y - YOFF, marker="o", ms=4.2, color="#2166AC",
             mec="white", mew=0.5, zorder=4)
    axB.plot([row["OR_vte_lo"], row["OR_vte_hi"]], [y + YOFF, y + YOFF],
             color="#B2182B", lw=1.05, solid_capstyle="round", zorder=3)
    axB.plot(row["OR_vte"], y + YOFF, marker="s", ms=4.2, color="#B2182B",
             mec="white", mew=0.5, zorder=4)

axB.axvline(1.0, color="#555555", lw=0.8, ls=(0, (4, 3)), zorder=1)
axB.set_xscale("log")
axB.set_xlim(0.70, 2.60)
axB.set_xticks([0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5])
axB.set_xticklabels(["0.7", "0.8", "0.9", "1.0", "1.2", "1.5", "2.0", "2.5"])
axB.xaxis.set_minor_locator(NullLocator())
axB.set_yticks(Y_B)
axB.set_yticklabels([f"{g} ({s})" for g, s in zip(b["gene"], b["SNP"])])
axB.set_ylim(0.35, 9.65)
axB.set_xlabel("OR per sentinel allele (95% CI)")
axB.set_title("B   Coagulation genes: sentinel-SNP ORs", loc="left",
              fontsize=9.5, fontweight="bold", pad=7)
for side in ("top", "right"):
    axB.spines[side].set_visible(False)

for lbl in axB.get_yticklabels():
    if lbl.get_text().startswith("F5 "):
        lbl.set_fontweight("bold")

handles = [
    Line2D([], [], marker="o", color="none", mfc="#2166AC", mec="white", ms=4.2,
           label="Endometriosis"),
    Line2D([], [], marker="s", color="none", mfc="#B2182B", mec="white", ms=4.2,
           label="VTE"),
    Line2D([], [], color="#E69F00", lw=1.0, label="Opposite sign"),
]
axB.legend(handles=handles, loc="lower right", frameon=False, fontsize=6.4,
           handlelength=1.5, labelspacing=0.35, borderaxespad=0.4)

# ------------------------------------------------------------------ Layer 4 checks
checks = []


def check(name, expected, observed):
    ok = (math.isclose(float(expected), float(observed), rel_tol=1e-9, abs_tol=1e-12)
          if isinstance(expected, (int, float)) and not isinstance(expected, bool)
          else expected == observed)
    checks.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {name:48s} {expected} vs {observed}")


print("\nLayer 4 reconciliation")
check("rows in Fig 3A", 8, len(a))
check("lineage order", ["thrombotic"] * 4 + ["gynaecological"] * 4, a["lineage"].tolist())
check("OR > 1 count", 4, int((a["OR"] > 1).sum()))
check("OR < 1 count", 4, int((a["OR"] < 1).sum()))
check("rounded ORs vs protocol Table 41", [3.28, 2.75, 2.35, 1.69, 0.92, 0.87, 0.79, 0.69],
      [round(v, 2) for v in a["OR"]])
check("all 8 associations P < 0.05", 8, int((a["p"] < 0.05).sum()))
for _, r in a.iterrows():
    check(f"OR = exp(beta) {r['phenocode']}", math.exp(r["beta"]), r["OR"])
    check(f"CI lower {r['phenocode']}",
          math.exp(r["beta"] - 1.959963984540054 * r["se"]), r["OR_lo"])
    check(f"CI upper {r['phenocode']}",
          math.exp(r["beta"] + 1.959963984540054 * r["se"]), r["OR_hi"])

check("coagulation genes", 9, len(b))
check("gene set", sorted(["F5", "F2", "ABO", "F11", "FGG", "PROCR", "KNG1", "THBD", "SERPINC1"]),
      sorted(b["gene"].tolist()))
for _, r in b.iterrows():
    check(f"OR_endo = exp(beta) {r['gene']}", math.exp(r["beta_endo"]), r["OR_endo"])
    check(f"OR_vte = exp(beta) {r['gene']}", math.exp(r["beta_vte"]), r["OR_vte"])
    check(f"endo CI lower {r['gene']}",
          math.exp(r["beta_endo"] - 1.96 * r["se_endo"]), r["OR_endo_lo"])
    check(f"vte CI lower {r['gene']}",
          math.exp(r["beta_vte"] - 1.96 * r["se_vte"]), r["OR_vte_lo"])

f5 = b[b["gene"] == "F5"].iloc[0]
check("F5 sentinel is rs6025", "rs6025", f5["SNP"])
f5a = a[a["phenocode"] == "N14_ENDOMETRIOSIS"].iloc[0]
check("F5 endo OR matches Fig 3A", f5a["OR"], f5["OR_endo"])
check("F5 VTE OR", 2.345822297, f5["OR_vte"])
check("opposite-sign genes", ["F5", "FGG", "PROCR", "THBD"],
      sorted(b.loc[b["opposite"], "gene"].tolist()))
sig_opp = b[(b["opposite"]) & (b["p_endo"] < 0.05) & (b["p_vte"] < 0.05)]
check("opposite + P<0.05 both sides", ["F5"], sig_opp["gene"].tolist())
check("no missing values", 0, int(a.isna().sum().sum() + b.isna().sum().sum()))

if not all(checks):
    raise SystemExit("VERIFICATION FAILED - figure not written")
print(f"\n  ALL {len(checks)} CHECKS PASSED")

# ------------------------------------------------------------------ source data
rows = []
for _, r in a.iterrows():
    rows.append(dict(panel="A", item=r["phenotype"], code=r["phenocode"],
                     n_case=r["n_case"], n_control=r["n_control"],
                     beta=r["beta"], se=r["se"], OR=r["OR"], OR_lo=r["OR_lo"],
                     OR_hi=r["OR_hi"], p=r["p"],
                     detail=r["finngen_phenostring"]))
for _, r in b.iterrows():
    rows.append(dict(panel="B", item=f"{r['gene']} ({r['SNP']})", code=r["SNP"],
                     n_case="", n_control="",
                     beta=f"{r['beta_endo']:.6g} / {r['beta_vte']:.6g}",
                     se=f"{r['se_endo']:.6g} / {r['se_vte']:.6g}",
                     OR=f"{r['OR_endo']:.6g} / {r['OR_vte']:.6g}",
                     OR_lo=f"{r['OR_endo_lo']:.6g} / {r['OR_vte_lo']:.6g}",
                     OR_hi=f"{r['OR_endo_hi']:.6g} / {r['OR_vte_hi']:.6g}",
                     p=f"{r['p_endo']:.6g} / {r['p_vte']:.6g}",
                     detail=f"opposite={r['opposite']}; n_variants={r['n_variants']}"))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "fig3_source_data.tsv"), sep="\t", index=False)

for ext in ("png", "pdf", "tiff"):
    path = os.path.join(FIG, f"fig3_f5_coag.{ext}")
    if ext == "tiff":
        fig.savefig(path, dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    else:
        fig.savefig(path, dpi=300)
print("saved: figures/fig3_f5_coag.png|pdf|tiff")
