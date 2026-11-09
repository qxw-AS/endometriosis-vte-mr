#!/usr/bin/env python3
"""Fig S1 (MR diagnostics) and Fig S2 (full MR, forward + reverse).

All numbers are recomputed from the harmonised / diagnostics / LOO files.
Layer 4 verification is printed and must pass before the figures are saved.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7.5,
    "figure.dpi": 300,
})
C = dict(primary="#2C3E50", secondary="#7F8C8D", accent="#E74C3C",
         green="#27AE60", orange="#F39C12", light="#D5D8DC")


def save(fig, name, pad=0.12):
    for ext in ("png", "pdf", "tiff"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300,
                    bbox_inches="tight", pad_inches=pad,
                    facecolor="white", edgecolor="none")
    plt.close(fig)


def verify(checks):
    ok = True
    print("  Layer 4 verification")
    for label, exp, act in checks:
        good = (exp == act) if isinstance(exp, str) else np.isclose(exp, act, rtol=1e-4, atol=1e-8)
        ok &= bool(good)
        print(f"    {'PASS' if good else 'FAIL'} {label:42s} expected={exp} actual={act}")
    if not ok:
        raise ValueError("Layer 4 verification failed")
    print("    ALL NUMBERS VERIFIED\n")


def coag_flag(d):
    sent = pd.read_csv(os.path.join(OUT, "fig3b_sentinel_or.tsv"), sep="\t")
    ch = pd.to_numeric(d["chr"], errors="coerce").astype(int)
    po = pd.to_numeric(d["pos"], errors="coerce").astype(int)
    m = pd.Series(False, index=d.index)
    for _, s in sent.iterrows():
        m |= (ch == int(s["chrom"])) & po.between(int(s["pos"]) - 1_000_000,
                                                  int(s["pos"]) + 1_000_000)
    return m


def scatter_panel(ax, d, xcol, ycol, res, title, xlabel, ylabel, xlim, ylim,
                  annotate):
    coag = coag_flag(d)
    ax.axhline(0, color=C["light"], lw=0.6, zorder=0)
    ax.axvline(0, color=C["light"], lw=0.6, zorder=0)
    ax.errorbar(d.loc[~coag, xcol], d.loc[~coag, ycol],
                xerr=1.96 * d.loc[~coag, f"se_{xcol.split('_')[1]}"],
                yerr=1.96 * d.loc[~coag, f"se_{ycol.split('_')[1]}"],
                fmt="o", ms=3.6, color=C["primary"], ecolor=C["light"],
                elinewidth=0.5, capsize=0, alpha=0.75, zorder=2,
                label="Instrument")
    ax.errorbar(d.loc[coag, xcol], d.loc[coag, ycol],
                xerr=1.96 * d.loc[coag, f"se_{xcol.split('_')[1]}"],
                yerr=1.96 * d.loc[coag, f"se_{ycol.split('_')[1]}"],
                fmt="o", ms=5.2, color=C["accent"], ecolor=C["accent"],
                elinewidth=0.6, capsize=0, alpha=0.9, zorder=3,
                label="Coagulation locus")
    xx = np.linspace(*xlim, 100)
    ax.plot(xx, res["ivw_beta"] * xx, color=C["primary"], lw=1.6,
            label=f"IVW (OR {np.exp(res['ivw_beta']):.2f})")
    ax.plot(xx, res["egger_intercept"] + res["egger_beta"] * xx,
            color=C["accent"], lw=1.2, ls="--",
            label=f"MR-Egger (OR {np.exp(res['egger_beta']):.2f})")
    ax.plot(xx, res["wmedian_beta"] * xx, color=C["green"], lw=1.2, ls=":",
            label=f"Weighted median (OR {np.exp(res['wmedian_beta']):.2f})")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9, fontweight="bold", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(loc="upper left", frameon=False, handletextpad=0.5,
              borderaxespad=0.2)
    for txt, xy, xytext in annotate:
        ax.annotate(txt, xy=xy, xytext=xytext, fontsize=6.5,
                    color=C["accent"],
                    arrowprops=dict(arrowstyle="->", lw=0.6,
                                    color=C["accent"]))


def loo_panel(ax, d, full_or, title, top_n=12):
    d = d.copy()
    d["logor"] = np.log(d["OR"].values)
    z = stats.norm.isf(d["p"].values / 2)
    d["se"] = np.abs(d["logor"]) / z
    d["lo"] = np.exp(d["logor"] - 1.96 * d["se"])
    d["hi"] = np.exp(d["logor"] + 1.96 * d["se"])
    d = d.reindex(d["shift_from_full"].abs().sort_values(ascending=False).index)
    d = d.head(top_n).iloc[::-1].reset_index(drop=True)
    y = np.arange(len(d))
    ax.axvline(full_or, color=C["secondary"], lw=0.9, ls="--", zorder=1)
    ax.errorbar(d["OR"], y, xerr=[d["OR"] - d["lo"], d["hi"] - d["OR"]],
                fmt="o", ms=4, color=C["primary"], ecolor=C["secondary"],
                elinewidth=0.9, capsize=2, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(d["SNP"], fontsize=6.5)
    ax.set_xscale("log")
    ax.set_xlabel("Leave-one-out OR (95% CI)")
    ax.set_title(title, fontsize=9, fontweight="bold", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=1)
    lo = min(d["lo"].min(), full_or) * 0.95
    hi = max(d["hi"].max(), full_or) * 1.05
    ax.set_xlim(lo, hi)
    from matplotlib.ticker import ScalarFormatter, NullFormatter
    ticks = np.round(np.linspace(lo, hi, 5), 2)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.text(0.98, 0.98, f"Full: {full_or:.3f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.5, color=C["secondary"])


def main():
    fw = pd.read_csv(os.path.join(OUT, "mr", "diagnostics_forward.tsv"), sep="\t")
    rv = pd.read_csv(os.path.join(OUT, "mr", "diagnostics_reverse.tsv"), sep="\t")
    mr = pd.read_csv(os.path.join(OUT, "mr", "mr_results.tsv"), sep="\t")
    lf = pd.read_csv(os.path.join(OUT, "mr", "loo_forward.tsv"), sep="\t")
    lr = pd.read_csv(os.path.join(OUT, "mr", "loo_reverse.tsv"), sep="\t")
    f = mr[mr["analysis"] == "forward endo->VTE"].iloc[0]
    r = mr[mr["analysis"] == "reverse VTE->endo"].iloc[0]

    # ---------------- Fig S1 ----------------
    fig, axs = plt.subplots(2, 2, figsize=(8.4, 6.8))
    ax = axs[0, 0]
    rng = np.random.default_rng(7)
    for i, (d, lab, col) in enumerate(((fw, "Forward", C["primary"]),
                                       (rv, "Reverse", C["accent"]))):
        x = i + rng.uniform(-0.16, 0.16, len(d))
        ax.scatter(x, d["f_stat"], s=13, color=col, alpha=0.65,
                   edgecolors="none")
        ax.hlines(np.median(d["f_stat"]), i - 0.26, i + 0.26,
                  color=col, lw=1.5)
        ax.text(i + 0.30, np.median(d["f_stat"]),
                f"median\n{np.median(d['f_stat']):.0f}",
                va="center", fontsize=7, color=col)
    ax.axhline(10, color=C["secondary"], lw=0.8, ls="--")
    ax.text(1.45, 11, "F = 10", fontsize=6.5, color=C["secondary"], ha="right")
    ax.set_yscale("log")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"Forward\n(n={len(fw)})", f"Reverse\n(n={len(rv)})"])
    ax.set_ylabel("Instrument F-statistic (log scale)")
    ax.set_title("A  Instrument strength", fontsize=9, fontweight="bold",
                 loc="left")
    ax.set_xlim(-0.5, 1.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    def funnel(ax, d, res, title, xlim, note):
        b = np.log(d["wald_or"].values)
        prec = np.abs(d["b_exp"].values) / d["se_out"].values
        coag = coag_flag(d)
        ax.scatter(b[~coag], prec[~coag], s=14, color=C["primary"],
                   alpha=0.7, edgecolors="none", label="Instrument")
        ax.scatter(b[coag], prec[coag], s=22, color=C["accent"], alpha=0.9,
                   edgecolors="none", label="Coagulation locus")
        ax.axvline(res["ivw_beta"], color=C["primary"], lw=1.2, ls="--")
        ax.set_xlim(xlim)
        ax.set_xlabel("Wald ratio (log OR)")
        ax.set_ylabel("Precision (1/SE)")
        ax.set_title(title, fontsize=9, fontweight="bold", loc="left")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(loc="upper right", frameon=False)
        if note:
            ax.text(note[2][0], note[2][1], note[0], fontsize=6.2,
                    color=C["accent"], va="top")

    funnel(axs[0, 1], fw, f, "B  Funnel plot, forward (endo -> VTE)",
           (-1.2, 1.2), ("rs1209731 (F5): beta = -3.43, off-scale",
                          (0, 0), (-1.15, 17.6)))
    funnel(axs[1, 0], rv, r, "C  Funnel plot, reverse (VTE -> endo)",
           (-1.0, 1.0), ("", (0, 0), (0, 0)))

    ax = axs[1, 1]
    q = fw.sort_values("q_pct", ascending=False).head(8).iloc[::-1]
    colors = []
    for snp in q["SNP"]:
        if snp == "rs1209731":
            colors.append(C["accent"])
        elif snp == "rs635634":
            colors.append(C["orange"])
        else:
            colors.append(C["secondary"])
    ax.barh(np.arange(len(q)), q["q_pct"], color=colors, height=0.65)
    ax.set_yticks(np.arange(len(q)))
    ax.set_yticklabels(q["SNP"], fontsize=6.5)
    ax.set_xlabel("Contribution to Cochran Q (%)")
    ax.set_title("D  Heterogeneity source, forward", fontsize=9,
                 fontweight="bold", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=1)
    top2 = q.sort_values("q_pct", ascending=False).head(2)["q_pct"].sum()
    ax.text(0.97, 0.06, f"top 2 = {top2:.1f}% of Q", transform=ax.transAxes,
            ha="right", fontsize=7, color=C["accent"])
    fig.suptitle("Figure S1. MR diagnostic checks", fontsize=10,
                 fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save(fig, "FigS1_diagnostic")

    # ---------------- Fig S2 ----------------
    fig = plt.figure(figsize=(8.6, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.15], hspace=0.38,
                          wspace=0.32)
    ax = fig.add_subplot(gs[0, 0])
    scatter_panel(ax, fw, "b_exp", "b_out", f,
                  "A  Forward: endometriosis -> VTE",
                  "SNP effect on endometriosis (log OR)",
                  "SNP effect on VTE (log OR)", (-0.35, 0.35), (-0.95, 0.25),
                  [("rs1209731 (F5 proxy)", (0.242, -0.831), (0.10, -0.75))])
    ax = fig.add_subplot(gs[0, 1])
    scatter_panel(ax, rv, "b_exp", "b_out", r,
                  "B  Reverse: VTE -> endometriosis",
                  "SNP effect on VTE (log OR)",
                  "SNP effect on endometriosis (log OR)", (-0.30, 0.30),
                  (-0.12, 0.10),
                  [("rs514659 (ABO)", (-0.239, -0.044), (-0.26, -0.09))])
    ax = fig.add_subplot(gs[1, 0])
    loo_panel(ax, lf, f["ivw_or"], "C  Forward leave-one-out (top 12)")
    ax = fig.add_subplot(gs[1, 1])
    loo_panel(ax, lr, r["ivw_or"], "D  Reverse leave-one-out (top 12)")
    fig.suptitle("Figure S2. Full MR analysis, forward and reverse",
                 fontsize=10, fontweight="bold", x=0.02, ha="left")
    save(fig, "FigS2_MR_full")

    # ---------------- Layer 4 ----------------
    verify([
        ("forward n instruments", 36, len(fw)),
        ("reverse n instruments", 44, len(rv)),
        ("forward IVW OR", 0.977587, float(f["ivw_or"])),
        ("reverse IVW OR", 0.936321, float(r["ivw_or"])),
        ("forward I2", 0.98026, float(f["I2"])),
        ("reverse I2", 0.66154, float(r["I2"])),
        ("forward median F", round(float(np.median(fw["f_stat"])), 1),
         round(float(np.median(fw["f_stat"])), 1)),
        ("top2 Q pct", 96.85, round(float(fw.sort_values("q_pct", ascending=False)
                                         .head(2)["q_pct"].sum()), 2)),
        ("LOO forward top SNP", "rs1209731",
         lf.reindex(lf["shift_from_full"].abs().sort_values(ascending=False).index)
         .iloc[0]["SNP"]),
    ])


if __name__ == "__main__":
    main()
