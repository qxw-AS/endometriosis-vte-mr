#!/usr/bin/env python3
"""M3 stage 5: integrate coloc + PLACO into the protocol deliverables.

High-confidence shared locus = coloc PP.H4 >= 0.8 AND PLACO p < 5e-8 in the
same candidate region. Suggestive tier = 0.5 <= PP.H4 < 0.8 AND PLACO p < 5e-8.
MHC is reported separately and excluded from the primary count. Direction
consistency is evaluated on the coloc top shared SNP (sign of beta).
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m3")

regions = pd.read_csv(os.path.join(OUT, "candidate_regions.tsv"), sep="\t",
                      dtype={"chr": str})
coloc = pd.read_csv(os.path.join(OUT, "coloc_results.tsv"), sep="\t",
                    dtype={"chr": str})
loci = pd.read_csv(os.path.join(OUT, "placo_loci.tsv"), sep="\t",
                   dtype={"chr38": "Int64"})

# region-level PLACO aggregation
agg = (loci.dropna(subset=["region"]).groupby("region")
       .agg(n_placo_loci=("SNP", "size"), placo_min_p=("P", "min"),
            placo_lead=("SNP", lambda s: loci.loc[s.index]
                        .sort_values("P").iloc[0]["SNP"]))
       .reset_index())
summ = regions.merge(coloc, on=["region", "chr", "start", "end"],
                     how="left", suffixes=("", "_coloc"))
summ = summ.merge(agg, on="region", how="left")
summ["n_placo_loci"] = summ["n_placo_loci"].fillna(0).astype(int)
summ["placo_min_p"] = summ["placo_min_p"].fillna(1.0)


def classify(r):
    has_coloc = r["PP.H4"] >= 0.8
    has_sugg = r["PP.H4"] >= 0.5
    has_placo = r["placo_min_p"] < 5e-8
    if has_coloc and has_placo:
        return "HIGH_CONFIDENCE"
    if has_sugg and has_placo:
        return "SUGGESTIVE_SHARED"
    if has_coloc:
        return "COLOC_ONLY"
    if has_placo:
        return "PLACO_ONLY"
    return "NO_EVIDENCE"


summ["class"] = summ.apply(classify, axis=1)
if "mhc" in summ.columns:
    summ["mhc"] = summ["mhc"].fillna(False).astype(bool)


def direction(r):
    if pd.isna(r["beta_endo"]) or pd.isna(r["beta_vte"]):
        return "NA"
    return "concordant" if np.sign(r["beta_endo"]) == np.sign(r["beta_vte"]) \
        else "discordant"


summ["direction"] = summ.apply(direction, axis=1)
summ.to_csv(os.path.join(OUT, "m3_region_summary.tsv"), sep="\t", index=False,
            float_format="%.6g")

non_mhc = summ[~summ["mhc"]]
high = non_mhc[non_mhc["class"] == "HIGH_CONFIDENCE"]
sugg = non_mhc[non_mhc["class"] == "SUGGESTIVE_SHARED"]
high.to_csv(os.path.join(OUT, "m3_high_confidence.tsv"), sep="\t", index=False,
            float_format="%.6g")


def binom_line(df, label):
    n = len(df)
    k = int((df["direction"] == "concordant").sum())
    if n == 0:
        return f"{label}: n=0, NOT_ASSESSED"
    p2 = stats.binomtest(k, n, 0.5, alternative="two-sided").pvalue
    p1_less = stats.binomtest(k, n, 0.5, alternative="less").pvalue
    return (f"{label}: n={n}, concordant={k}/{n} = {k/n:.1%}, "
            f"binomial two-sided p={p2:.4g} (one-sided depletion p={p1_less:.4g})")


# PLACO-only loci outside candidate regions
placo_only = loci[loci["region"].isna()].copy()
placo_only_dir = placo_only.dropna(subset=["beta_endo", "beta_vte"])
placo_only_dir = placo_only_dir.assign(
    direction=np.where(np.sign(placo_only_dir["beta_endo"]) ==
                       np.sign(placo_only_dir["beta_vte"]),
                       "concordant", "discordant"))

# coagulation gene panel (M6.1): coloc PP from the containing region, PLACO min p
# within sentinel +/- 100 kb (GRCh38). The PLACO hit table stores float64 exact
# p-values; the all-score npz is float32 and underflows for extreme p, so it is
# only used as a fallback.
sent = pd.read_csv(os.path.join(ROOT, "out", "fig3b_sentinel_or.tsv"), sep="\t")
hits = pd.read_csv(os.path.join(OUT, "placo_hits_p1e5.tsv.gz"), sep="\t")

rows = []
for _, g in sent.iterrows():
    c = str(int(g["chrom"]))
    pos = int(float(g["pos"]))
    h = hits[(hits["chr"].astype(str) == c) &
             (hits["pos"] >= pos - 100_000) & (hits["pos"] <= pos + 100_000)]
    pstr = f"{h['placo_p'].min():.3g}" if len(h) else ">1e-5"
    reg = regions[(regions["chr"] == c) & (regions["start"] <= pos) &
                  (regions["end"] >= pos)]
    row = dict(gene=g["gene"], sentinel=g["SNP"], chr=c, pos=pos,
               region=reg.iloc[0]["region"] if not reg.empty else "NOT_IN_CANDIDATE",
               placo_min_p_gene_100kb=pstr,
               sentinel_beta_endo=float(g["beta_endo"]),
               sentinel_beta_vte=float(g["beta_vte"]),
               sentinel_direction=("concordant" if np.sign(g["beta_endo"]) ==
                                   np.sign(g["beta_vte"]) else "discordant"))
    if not reg.empty:
        cr = coloc[coloc["region"] == reg.iloc[0]["region"]].iloc[0]
        for i in range(5):
            row[f"PP.H{i}"] = cr[f"PP.H{i}"]
        row["top_shared"] = cr["top_shared"]
        row["beta_endo"] = cr["beta_endo"]
        row["beta_vte"] = cr["beta_vte"]
        row["direction"] = "concordant" if np.sign(cr["beta_endo"]) == \
            np.sign(cr["beta_vte"]) else "discordant"
    else:
        for i in range(5):
            row[f"PP.H{i}"] = "NOT_ASSESSED"
        row["top_shared"] = "NOT_ASSESSED"
        row["beta_endo"] = row["beta_vte"] = np.nan
        row["direction"] = "NA"
    rows.append(row)
coag = pd.DataFrame(rows)
coag.to_csv(os.path.join(OUT, "coag_loci_coloc.tsv"), sep="\t", index=False,
            float_format="%.6g")

# sdY sensitivity range for headline regions
sens = pd.read_csv(os.path.join(OUT, "coloc_sensitivity.tsv"), sep="\t")
rng = sens.groupby("region")["PP_H4"].agg(["min", "max"]).reset_index()
summ = summ.merge(rng.rename(columns={"min": "PP.H4_min_sdy",
                                      "max": "PP.H4_max_sdy"}),
                  on="region", how="left")
non_mhc = summ[~summ["mhc"]]
high = non_mhc[non_mhc["class"] == "HIGH_CONFIDENCE"]
sugg = non_mhc[non_mhc["class"] == "SUGGESTIVE_SHARED"]
high.to_csv(os.path.join(OUT, "m3_high_confidence.tsv"), sep="\t", index=False,
            float_format="%.6g")

lines = []
lines.append("M3 共定位与共享因果变异（V6.2 方案）")
lines.append("=" * 60)
lines.append(f"候选区域：{len(regions)} 个（37 个内异症 lead + 46 个 VTE lead，"
             f"±500 kb 并集；MHC {int(regions['mhc'].sum())} 个单列）")
lines.append(f"LAVA：NOT_RUN（环境无 R/LAVA）；按方案阶段 2 回退到 lead ±500 kb 窗口")
lines.append(f"coloc.abf：sdY=0.15，p1=p2=1e-4，p12=1e-5；敏感性 sdY=0.10/0.20/0.30")
params = pd.read_csv(os.path.join(OUT, "placo_params.txt"), sep="\t",
                     header=None).set_index(0)[1]
lines.append(f"PLACO：全基因组 {int(params['n_variants']):,} 个常见变异；"
             f"rho={float(params['rho']):.4f}；"
             f"p<5e-8 的 SNP {int(params['n_p_lt_5e8']):,} 个")
lines.append("")
lines.append(f"PP.H4 >= 0.8 的区域：{int((non_mhc['PP.H4'] >= 0.8).sum())} 个"
             f"（MHC 另计：{int((summ['mhc'] & (summ['PP.H4'] >= 0.8)).sum())}）")
lines.append(f"PP.H4 0.5-0.8：{int(((non_mhc['PP.H4'] >= 0.5) & (non_mhc['PP.H4'] < 0.8)).sum())} 个")
lines.append("")
lines.append("高置信共享位点（PP.H4>=0.8 且 PLACO p<5e-8，MHC 除外）：")
lines.append(high[["region", "chr", "start", "end", "n_snp", "PP.H3", "PP.H4",
                   "top_shared", "beta_endo", "beta_vte", "direction",
                   "placo_min_p", "PP.H4_min_sdy", "PP.H4_max_sdy"]]
             .to_string(index=False))
lines.append("")
lines.append("提示级共享位点（0.5<=PP.H4<0.8 且 PLACO p<5e-8）：")
lines.append(sugg[["region", "chr", "start", "end", "n_snp", "PP.H3", "PP.H4",
                   "top_shared", "beta_endo", "beta_vte", "direction",
                   "placo_min_p", "PP.H4_min_sdy", "PP.H4_max_sdy"]]
             .to_string(index=False))
lines.append("  注：R55(PROCR) 的 PP.H4 随 sdY 在 0.58（0.10）至 0.36（0.30）间移动，"
             "跨过 0.5 阈值，只能算提示级；R03(F5) 在全部 sdY 下 PP.H4>=0.993。")
lines.append("")
lines.append("方向一致率：")
lines.append("  " + binom_line(high, "高置信位点"))
two = pd.concat([high, sugg])
lines.append("  " + binom_line(two, "高置信+提示级"))
lines.append("  " + binom_line(summ[(summ["PP.H4"] >= 0.5) & (~summ["mhc"])],
                                "所有 PP.H4>=0.5 区域"))
lines.append("  " + binom_line(placo_only_dir,
                                "候选区域外的 PLACO 位点（探索性）"))
lines.append("")
lines.append("PLACO 位点总数：%d（区域内 %d，区域外 %d；MHC %d）"
             % (len(loci), loci["region"].notna().sum(), len(placo_only),
                int(loci["mhc"].fillna(False).sum())))
lines.append("")
lines.append("凝血基因面板（M6.1；PP 来自所在候选区域，PLACO 为 sentinel ±100 kb）：")
lines.append(coag[["gene", "sentinel", "region", "PP.H0", "PP.H1", "PP.H2",
                   "PP.H3", "PP.H4", "top_shared", "beta_endo", "beta_vte",
                   "direction", "sentinel_beta_endo", "sentinel_beta_vte",
                   "sentinel_direction", "placo_min_p_gene_100kb"]]
             .to_string(index=False))
lines.append("")
lines.append("MHC 区域（单列，不进主计数）：")
mhc = summ[summ["mhc"]]
lines.append(mhc[["region", "chr", "start", "end", "PP.H3", "PP.H4",
                  "placo_min_p", "n_placo_loci"]].to_string(index=False))
lines.append("")
lines.append("方法披露：")
lines.append("- coloc.abf 用 Python 复现（Wakefield ABF + 标准 H0-H4 公式），"
             "已用模拟数据自检（共享/独立/单性状场景）。")
lines.append("- PLACO 按 Ray & Chatterjee 2020 S1 Appendix Eq.3 的近似 p 值实现，"
             "Z 相关由 FinnGen 共享对照的空变异估计并做 R^-1/2 去相关；"
             "40 万对独立正态 Z 自检：P(p<0.05)=0.0494，P(p<1e-4)=1.35e-4"
             "（后者为 40/40 万期望值附近的 2.2 SD 波动）。")
lines.append("- 两病均用 FinnGen R13 双侧汇总统计，对照高度重叠；coloc 的独立样本假设"
             "不完全满足，PP.H4 可能偏乐观，结论上限 C2。")
lines.append("- coloc 假设每个性状在区域内只有一个因果变异；ABO 等强信号区域"
             "PP.H3 占优，可能反映两个独立信号，也可能是单因果假设的局限。")
lines.append("- 未做 LAVA（无 R）；未做 SuSiE 精细定位。")

report = "\n".join(lines)
with open(os.path.join(OUT, "m3_report.txt"), "w") as fh:
    fh.write(report + "\n")
print(report)

# refresh summary with sensitivity columns
summ.to_csv(os.path.join(OUT, "m3_region_summary.tsv"), sep="\t", index=False,
            float_format="%.6g")
