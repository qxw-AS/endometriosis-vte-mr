#!/usr/bin/env python3
"""Figure S4: mirrored LocusZoom plots for the two shared coagulation loci.

Columns: F5 (R03) and ABO (R28).  Top row = endometriosis, bottom row = VTE.
Points are coloured by LD r2 with the shared lead variant computed in 1000G
EUR Phase 3 (PLINK 1.9), because FinnGen R13 is GRCh38 and 1000G is GRCh37;
LD is therefore merged on rsID, never on position.
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter, NullFormatter
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.8,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "figure.dpi": 300,
})

COL = {
    "shared": "#D73027",
    "endo": "#1A1A1A",
    "gene": "#2C3E50",
    "grid": "#E6E6E6",
    "sex": "#D9D9D9",
}

LD_BINS = [
    (0.8, 1.0001, "#D73027"),
    (0.6, 0.8, "#FDAE61"),
    (0.4, 0.6, "#ABD9E9"),
    (0.2, 0.4, "#74ADD1"),
    (-0.0001, 0.2, "#4575B4"),
    (None, None, "#BDBDBD"),   # no LD information
]

LOCI = [
    dict(name="F5", region="R03", chrom=1, start=168117123, end=170238202,
         sentinel="rs6025", ld="f5_ld.ld", endo_lead="rs1209731",
         gene_start=169514166, gene_end=169586588, ymax_vte=260,
         xlab="Position on chr1 (Mb)", coloc="PP.H4 = 0.994"),
    dict(name="ABO", region="R28", chrom=9, start=132493162, end=133779427,
         sentinel="rs529565", ld="abo_ld.ld", endo_lead="rs635634",
         gene_start=133250375, gene_end=133275204, ymax_vte=180,
         xlab="Position on chr9 (Mb)", coloc="PP.H3 = 0.953"),
]

ENDO_FILE = os.path.join(OUT, "m3", "region_snps_endo.tsv.gz")
VTE_FILE = os.path.join(OUT, "m3", "region_snps_vte.tsv.gz")


def verify(checks):
    ok = True
    for label, exp, act in checks:
        good = (exp == act) if isinstance(exp, str) else np.isclose(exp, act, rtol=1e-4, atol=1e-8)
        ok &= bool(good)
        print(f"    {'PASS' if good else 'FAIL'} {label:46s} expected={exp} actual={act}")
    if not ok:
        raise ValueError("Layer 4 verification failed")


def load_region(fn, chrom, start, end):
    d = pd.read_csv(fn, sep="\t")
    d = d[(d["#chrom"] == chrom) & (d["pos"].between(start, end))].copy()
    d["rsid"] = d["rsids"].str.split(",").str[0]
    z = d["beta"] / d["sebeta"]
    # FinnGen stores p=0 when the Wald p underflows; rebuild from beta/se.
    logp = np.where(d["pval"] > 0, -np.log10(d["pval"]),
                    -np.log10(np.maximum(2 * stats.norm.sf(np.abs(z)), 1e-320)))
    d["logp"] = logp
    return d


def ld_colors(snps, ld, anchor):
    ldd = ld.set_index("SNP_B")["R2"]
    r2 = snps["rsid"].map(ldd)
    colors = np.full(len(snps), LD_BINS[-1][2], dtype=object)
    for lo, hi, c in LD_BINS[:-1]:
        colors[(r2 > lo) & (r2 <= hi)] = c
    return r2.values, colors


def panel(ax, d, r2, colors, ymax, ylabel, sentinel, endo_lead,
          locus_start, locus_end, show_x, threshold=True):
    ax.scatter(d["pos"] / 1e6, d["logp"], s=5, c=list(colors),
               linewidths=0, rasterized=True, zorder=2)
    if threshold:
        ax.axhline(-np.log10(5e-8), color="#7F8C8D", lw=0.6, ls=(0, (4, 3)), zorder=1)
    # shared lead
    srow = d[d["rsid"] == sentinel]
    if len(srow):
        ax.scatter(srow["pos"] / 1e6, srow["logp"], marker="D", s=34,
                   facecolor=COL["shared"], edgecolor="black", linewidth=0.6, zorder=5)
    erow = d[d["rsid"] == endo_lead]
    if len(erow) and endo_lead != sentinel:
        ax.scatter(erow["pos"] / 1e6, erow["logp"], marker="^", s=40,
                   facecolor="#FFFFFF", edgecolor=COL["endo"], linewidth=1.0, zorder=5)
    ax.set_xlim(locus_start / 1e6, locus_end / 1e6)
    ax.set_ylim(0, ymax)
    ax.set_ylabel(ylabel, labelpad=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=2.5, width=0.6)
    if not show_x:
        ax.tick_params(labelbottom=False)
    else:
        fmt = ScalarFormatter(useOffset=False)
        fmt.set_scientific(False)
        ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.set_xlabel("Position (Mb)")


def main():
    print("=" * 78)
    print("Figure S4 - mirrored LocusZoom (F5 and ABO)")
    print("=" * 78)

    endo = pd.read_csv(ENDO_FILE, sep="\t")
    vte = pd.read_csv(VTE_FILE, sep="\t")
    coloc = pd.read_csv(os.path.join(OUT, "m3", "coag_loci_coloc.tsv"), sep="\t")
    for c in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4",
              "placo_min_p_gene_100kb"):
        coloc[c] = pd.to_numeric(coloc[c], errors="coerce")
    coloc = coloc.set_index("gene")

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(3, 2, height_ratios=[3.0, 3.0, 0.70],
                          hspace=0.10, wspace=0.24,
                          left=0.075, right=0.985, top=0.915, bottom=0.155)

    # drop duplicated header rows (region files concatenate phenotypes)
    endo = endo[endo["#chrom"] != "#chrom"].copy()
    vte = vte[vte["#chrom"] != "#chrom"].copy()
    for c in ("#chrom", "pos", "pval", "beta", "sebeta"):
        endo[c] = pd.to_numeric(endo[c], errors="coerce")
        vte[c] = pd.to_numeric(vte[c], errors="coerce")
    endo = endo.dropna(subset=["#chrom", "pos", "pval", "beta", "sebeta"])
    vte = vte.dropna(subset=["#chrom", "pos", "pval", "beta", "sebeta"])
    for c in ("#chrom", "pos"):
        endo[c] = endo[c].astype(int)
        vte[c] = vte[c].astype(int)

    checks = []
    panel_handles = []

    for j, loc in enumerate(LOCI):
        de = load_region(ENDO_FILE, loc["chrom"], loc["start"], loc["end"])
        dv = load_region(VTE_FILE, loc["chrom"], loc["start"], loc["end"])
        ld = pd.read_csv(os.path.join(FIG, loc["ld"]), sep=r"\s+")

        r2e, ce = ld_colors(de, ld, loc["sentinel"])
        r2v, cv = ld_colors(dv, ld, loc["sentinel"])

        ax_e = fig.add_subplot(gs[0, j])
        ax_v = fig.add_subplot(gs[1, j], sharex=ax_e)
        ax_g = fig.add_subplot(gs[2, j], sharex=ax_e)

        panel(ax_e, de, r2e, ce, 10.2, "Endometriosis  −log$_{10}$(P)",
              loc["sentinel"], loc["endo_lead"], loc["start"], loc["end"], show_x=False)
        panel(ax_v, dv, r2v, cv, loc["ymax_vte"], "VTE  −log$_{10}$(P)",
              loc["sentinel"], loc["endo_lead"], loc["start"], loc["end"], show_x=False)

        # locus titles and coloc annotation
        ax_e.set_title(f"{loc['name']} locus  (chr{loc['chrom']}:"
                       f"{loc['start']/1e6:.1f}–{loc['end']/1e6:.1f} Mb)",
                       fontsize=9, fontweight="bold", pad=5, loc="left")
        ax_e.text(0.985, 0.92, f"coloc {loc['coloc']}", transform=ax_e.transAxes,
                  ha="right", va="top", fontsize=7.2, color="#333333")

        # gene track
        ax_g.set_ylim(0, 1)
        ax_g.set_xlim(loc["start"] / 1e6, loc["end"] / 1e6)
        ax_g.add_patch(Rectangle((loc["gene_start"] / 1e6, 0.30),
                                 (loc["gene_end"] - loc["gene_start"]) / 1e6, 0.34,
                                 facecolor=COL["gene"], edgecolor="none", zorder=2))
        ax_g.text((loc["gene_end"] / 1e6) + 0.012, 0.47, loc["name"],
                  ha="left", va="center", fontsize=7, fontstyle="italic")
        ax_g.set_yticks([])
        ax_g.set_xlabel(loc["xlab"])
        ax_g.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax_g.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax_g.spines["left"].set_visible(False)
        ax_g.spines["right"].set_visible(False)
        ax_g.spines["top"].set_visible(False)
        ax_g.spines["bottom"].set_linewidth(0.6)
        ax_g.tick_params(axis="x", direction="out", length=2.5, width=0.6)

        # Layer 4 numbers per locus
        exp = {
            "F5": dict(n_e=16704, n_v=16742, n_ld=22903, max_e=8.65, max_v=249.50,
                       s_e=8.51, s_v=249.50),
            "ABO": dict(n_e=10730, n_v=10744, n_ld=21963, max_e=7.81, max_v=161.87,
                        s_e=4.16, s_v=161.87),
        }[loc["name"]]
        checks += [
            (f"{loc['name']}: endo region SNPs", exp["n_e"], len(de)),
            (f"{loc['name']}: VTE region SNPs", exp["n_v"], len(dv)),
            (f"{loc['name']}: LD pairs", exp["n_ld"], len(ld)),
            (f"{loc['name']}: LD anchor max r2", 1.0, float(ld["R2"].max())),
            (f"{loc['name']}: sentinel log10P endo", exp["s_e"],
             round(float(de.loc[de["rsid"] == loc["sentinel"], "logp"].iloc[0]), 2)),
            (f"{loc['name']}: sentinel log10P VTE", exp["s_v"],
             round(float(dv.loc[dv["rsid"] == loc["sentinel"], "logp"].iloc[0]), 2)),
            (f"{loc['name']}: max endo -log10P", exp["max_e"], round(float(de["logp"].max()), 2)),
            (f"{loc['name']}: max VTE -log10P", exp["max_v"], round(float(dv["logp"].max()), 2)),
        ]
        print(f"\n[{loc['name']}] endo SNPs={len(de)}, VTE SNPs={len(dv)}, LD pairs={len(ld)}")
        print(f"  endo lead: {de.loc[de['logp'].idxmax(), 'rsid']} "
              f"log10P={de['logp'].max():.2f}; VTE lead: "
              f"{dv.loc[dv['logp'].idxmax(), 'rsid']} log10P={dv['logp'].max():.2f}")
        print(f"  coloc {loc['name']}: PP.H4={coloc.loc[loc['name'], 'PP.H4']:.6f}, "
              f"PP.H3={coloc.loc[loc['name'], 'PP.H3']:.6f}")

    # coloc assertions against the source table
    checks += [
        ("coloc F5 PP.H4", 0.993662, round(float(coloc.loc["F5", "PP.H4"]), 6)),
        ("coloc ABO PP.H3", 0.952859, round(float(coloc.loc["ABO", "PP.H3"]), 6)),
        ("coloc F5 PLACO min p", 2.35e-82, float(coloc.loc["F5", "placo_min_p_gene_100kb"])),
    ]
    print("\nLayer 4 verification:")
    verify(checks)

    # legend
    handles = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=c,
               markeredgecolor="none", markersize=6, label=lab)
        for lab, c in [
            ("LD r$^2$ 0.8–1.0", "#D73027"),
            ("0.6–0.8", "#FDAE61"),
            ("0.4–0.6", "#ABD9E9"),
            ("0.2–0.4", "#74ADD1"),
            ("<0.2", "#4575B4"),
            ("no LD data", "#BDBDBD"),
        ]
    ]
    handles += [
        Line2D([0], [0], marker="D", color="none", markerfacecolor=COL["shared"],
               markeredgecolor="black", markersize=6, label="Shared lead"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor="white",
               markeredgecolor=COL["endo"], markersize=7, label="Endometriosis lead"),
        Line2D([0], [0], color="#7F8C8D", lw=0.7, ls=(0, (4, 3)),
               label="P = 5×10$^{-8}$"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=9, frameon=False,
               handletextpad=0.35, columnspacing=1.1, bbox_to_anchor=(0.53, 0.015))

    base = os.path.join(FIG, "FigS4_LocusZoom")
    for ext in ("png", "pdf", "tiff"):
        fig.savefig(f"{base}.{ext}", dpi=300, bbox_inches="tight",
                    pad_inches=0.12, facecolor="white", edgecolor="none")
    plt.close(fig)
    print("\nSaved: FigS4_LocusZoom.{png,pdf,tiff}")


if __name__ == "__main__":
    main()
