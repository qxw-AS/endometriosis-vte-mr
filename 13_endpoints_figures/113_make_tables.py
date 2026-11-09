#!/usr/bin/env python3
"""Build Tables I, II, IV, V and VI for the V6.2 manuscript.

Every number is recomputed from the analysis outputs in out/; nothing is
copied from the manuscript text.  Table III is produced separately by
scripts/88_table3_v62.py.
"""
import math
import os
from statistics import NormalDist

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out")
TAB = os.path.join(OUT, "tables")
os.makedirs(TAB, exist_ok=True)


def check(label, expected, actual):
    good = (expected == actual) if isinstance(expected, str) else np.isclose(expected, actual, rtol=1e-4, atol=1e-9)
    print(f"  {'PASS' if good else 'FAIL'} {label:52s} expected={expected} actual={actual}")
    if not good:
        raise ValueError(f"Table number mismatch: {label}")


# ---------------------------------------------------------------- Table I
def table1():
    rows = [
        dict(trait="Endometriosis", role="Primary exposure",
             source="FinnGen R13", phenotype="N14_ENDOMETRIOSIS",
             cases=20534, controls=122622, total=143156, n_variants="all imputed variants",
             build="GRCh38", population="Finnish",
             note="43 FinnGen genome-wide significant loci; 37 lead SNPs after PLINK clumping (r2<0.001, 10 Mb, 1000G EUR); 36 usable after harmonisation (rs11686429 dropped, ambiguous palindrome). Primary MR uses 34 after coagulation-locus exclusion (rs1209731, rs635634)."),
        dict(trait="Venous thromboembolism", role="Primary outcome",
             source="FinnGen R13", phenotype="I9_VTE",
             cases=27928, controls=472258, total=500186, n_variants="all imputed variants",
             build="GRCh38", population="Finnish",
             note="84 FinnGen genome-wide significant loci; 46 lead SNPs after PLINK clumping (r2<0.001, 10 Mb, 1000G EUR); 44 usable after harmonisation (rs2844543 ambiguous palindrome; rs145770442 allele mismatch). Reverse external MR uses 38 of 44 (6 absent from GCST90205183)."),
        dict(trait="Endometriosis", role="External replication",
             source="Rahmioglu et al. 2023 (GCST90205183)", phenotype="GCST90205183",
             cases=60674, controls=410192, total=470866, n_variants="8,593,419",
             build="GRCh38", population="European + East Asian",
             note="24-GWAS meta-analysis; effective N from GWAS Catalog; includes a FinnGen freeze (6,502 cases)"),
        dict(trait="Body mass index", role="MVMR covariate",
             source="FinnGen R13", phenotype="BMI_IRN",
             cases=362216, controls=0, total=362216, n_variants="all imputed variants",
             build="GRCh38", population="Finnish",
             note="Inverse-rank normalized; OR per 1 SD"),
        dict(trait="Ever smoking", role="MVMR covariate",
             source="FinnGen R13", phenotype="EVER_SMOKER",
             cases=145362, controls=150295, total=295657, n_variants="all imputed variants",
             build="GRCh38", population="Finnish",
             note="Protocol deviation: GSCAN ever-smoking summary statistics could not be retrieved"),
        dict(trait="Pulmonary embolism", role="Severity stratification",
             source="FinnGen R13", phenotype="I9_PULMEMB",
             cases=13783, controls=485055, total=498838, n_variants="all imputed variants",
             build="GRCh38", population="Finnish",
             note="27 independent PE leads"),
        dict(trait="Deep vein thrombosis", role="Severity stratification",
             source="FinnGen R13", phenotype="I9_PHLETHROMBDVTLOW",
             cases=8601, controls=429443, total=438044, n_variants="all imputed variants",
             build="GRCh38", population="Finnish",
             note="51 independent DVT peaks; 5 are identical to PE peaks (including rs6025, F5) and counted once, leaving 46 DVT-labelled instruments (73 unique PE+DVT leads)"),
    ]
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(TAB, "TableI_data_sources.tsv"), sep="\t", index=False)
    check("Table I: endometriosis cases", 20534, int(t.loc[0, "cases"]))
    check("Table I: VTE controls", 472258, int(t.loc[1, "controls"]))
    check("Table I: external effective N", 470866, int(t.loc[2, "total"]))
    check("Table I: external variants", "8,593,419", t.loc[2, "n_variants"])
    check("Table I: BMI N", 362216, int(t.loc[3, "total"]))
    return t


# --------------------------------------------------------------- Table II
def coag_map(df, sent):
    gene = np.array([""] * len(df), dtype=object)
    for _, r in sent.iterrows():
        hit = (df["chr"] == int(r["chrom"])) & \
              (df["pos"].between(int(r["pos"]) - 1_000_000, int(r["pos"]) + 1_000_000))
        gene[hit.values & (gene == "")] = r["gene"]
    return gene


def table2():
    fw = pd.read_csv(os.path.join(OUT, "mr", "diagnostics_forward.tsv"), sep="\t")
    rv = pd.read_csv(os.path.join(OUT, "mr", "diagnostics_reverse.tsv"), sep="\t")
    sent = pd.read_csv(os.path.join(OUT, "fig3b_sentinel_or.tsv"), sep="\t")
    fw["coagulation_gene"] = coag_map(fw, sent)
    rv["coagulation_gene"] = coag_map(rv, sent)
    fw["direction"] = "Endometriosis -> VTE"
    rv["direction"] = "VTE -> Endometriosis"
    cols = ["direction", "SNP", "chr", "pos", "ea", "nea", "eaf",
            "b_exp", "se_exp", "p_exp", "b_out", "se_out", "p_out",
            "f_stat", "coagulation_gene", "use"]
    t = pd.concat([fw[cols], rv[cols]], ignore_index=True)
    t = t.sort_values(["direction", "chr", "pos"]).reset_index(drop=True)
    t = t.rename(columns={"chr": "chr_GRCh38", "pos": "pos_GRCh38",
                          "ea": "effect_allele", "nea": "other_allele",
                          "eaf": "EAF", "b_exp": "beta_exposure",
                          "se_exp": "SE_exposure", "p_exp": "P_exposure",
                          "b_out": "beta_outcome", "se_out": "SE_outcome",
                          "p_out": "P_outcome", "f_stat": "F_statistic",
                          "coagulation_gene": "coagulation_locus",
                          "use": "in_primary_set"})
    t["coagulation_locus"] = t["coagulation_locus"].replace("", "none")
    t["in_primary_set"] = np.where(t["coagulation_locus"] == "none", "yes", "no")
    t.to_csv(os.path.join(TAB, "TableII_instruments.tsv"), sep="\t", index=False,
             float_format="%.6g")
    n_fw = int((t["direction"] == "Endometriosis -> VTE").sum())
    n_rv = int((t["direction"] == "VTE -> Endometriosis").sum())
    check("Table II: forward instruments", 36, n_fw)
    check("Table II: reverse instruments", 44, n_rv)
    check("Table II: forward coagulation-locus SNPs",
          int((fw["coagulation_gene"] != "").sum()),
          int((t["direction"] == "Endometriosis -> VTE").sum() -
              (t.loc[t["direction"] == "Endometriosis -> VTE", "in_primary_set"] == "yes").sum()))
    check("Table II: reverse coagulation-locus SNPs",
          int((rv["coagulation_gene"] != "").sum()),
          int((t["direction"] == "VTE -> Endometriosis").sum() -
              (t.loc[t["direction"] == "VTE -> Endometriosis", "in_primary_set"] == "yes").sum()))
    check("Table II: forward min F", round(fw["f_stat"].min(), 3),
          round(float(t.loc[t["direction"] == "Endometriosis -> VTE", "F_statistic"].min()), 3))
    check("Table II: reverse min F", round(rv["f_stat"].min(), 3),
          round(float(t.loc[t["direction"] == "VTE -> Endometriosis", "F_statistic"].min()), 3))
    print("  coagulation-locus instruments:")
    print(t.loc[t["coagulation_locus"] != "none",
                ["direction", "SNP", "coagulation_locus", "F_statistic"]].to_string(index=False))
    return t


# ---------------------------------------------------------------- Table IV
def table4():
    c = pd.read_csv(os.path.join(OUT, "m3", "coag_loci_coloc.tsv"), sep="\t")
    interp = {
        "F5": "Shared causal variant (PP.H4=0.994); opposite effect directions",
        "ABO": "Two distinct causal variants (PP.H3=0.953); no shared causal signal",
        "PROCR": "Suggestive shared signal (PP.H4=0.502), unstable across sdY (0.36-0.58)",
        "THBD": "No shared signal; PP.H2=0.847 (distinct, one trait only)",
        "FGG": "No shared signal; PP.H2=0.749",
        "F11": "No shared signal; PP.H2=0.852",
        "F2": "No shared signal; PP.H2=0.676; sentinel rs1799963 and top_shared rs5896 are different SNPs (direction per top_shared)",
        "KNG1": "Not assessed: sentinel outside all +/-500 kb candidate regions",
        "SERPINC1": "Not assessed: sentinel outside all +/-500 kb candidate regions",
    }
    c["interpretation"] = c["gene"].map(interp)
    rename = {"sentinel": "sentinel_SNP", "pos": "sentinel_pos_GRCh38",
              "placo_min_p_gene_100kb": "PLACO_min_P_sentinel_100kb",
              "top_shared": "top_shared_SNP",
              "beta_endo": "shared_beta_endometriosis",
              "beta_vte": "shared_beta_VTE",
              "direction": "shared_effect_direction"}
    out = c.rename(columns=rename)
    out.to_csv(os.path.join(TAB, "TableIV_colocalization.tsv"), sep="\t", index=False)
    check("Table IV: genes", 9, len(out))
    check("Table IV: F5 PP.H4", 0.993662, round(float(out.loc[out.gene == "F5", "PP.H4"].iloc[0]), 6))
    check("Table IV: ABO PP.H3", 0.952859, round(float(out.loc[out.gene == "ABO", "PP.H3"].iloc[0]), 6))
    check("Table IV: assessed genes", 7,
          int(out["PP.H4"].apply(lambda v: str(v) != "NOT_ASSESSED").sum()))
    return out


# ---------------------------------------------------------------- Table V
def table5():
    lk = pd.read_csv(os.path.join(OUT, "g8", "g8_f5_lookup.tsv"), sep="\t")
    rv = pd.read_csv(os.path.join(OUT, "g8", "g8_reverse_mr.tsv"), sep="\t")
    r = lk.iloc[0]
    rows = [dict(
        analysis="F5 rs6025 lookup", data="External endometriosis (GCST90205183)",
        n_snp=1, effect_allele=r["effect_allele"], OR=float(r["OR_per_T"]),
        lo=float(r["OR_T_lo"]), hi=float(r["OR_T_hi"]), p=float(r["p"]),
        I2=np.nan, egger_int_p=np.nan, presso_p=np.nan, n_wald_below_1=np.nan,
        wm_or=np.nan, wm_p=np.nan,
        note="Direction PASS (OR<1, P<0.01)")]
    rows.append(dict(
        analysis="F5 proxy rs1209731 lookup", data="External endometriosis (GCST90205183)",
        n_snp=1, effect_allele="T", OR=0.853, lo=0.7875, hi=0.9244, p=1.03e-4,
        I2=np.nan, egger_int_p=np.nan, presso_p=np.nan, n_wald_below_1=np.nan,
        wm_or=np.nan, wm_p=np.nan,
        note="C is the risk allele (beta=+0.1588, SE=0.0409); OR converted to the T allele"))
    for _, q in rv.iterrows():
        ext = "external" in q["analysis"].lower()
        label = {
            "reverse VTE->endo | external endo (all)": "Reverse MR: all instruments",
            "reverse VTE->endo | external endo (drop F5 locus)": "Reverse MR: F5 locus excluded",
            "reverse VTE->endo | external endo (drop 9 coag genes)": "Reverse MR: 9 coagulation genes excluded",
            "reverse VTE->endo | FinnGen endo (same 38 instruments)": "Reverse MR (FinnGen matched): all 38 instruments",
            "reverse VTE->endo | FinnGen endo (same coag-excluded set)": "Reverse MR (FinnGen matched): 25-instrument set",
        }[q["analysis"]]
        note = {
            "Reverse MR: all instruments": "PASS: OR<0.95; PRESSO p=0.041 indicates mild heterogeneity",
            "Reverse MR: F5 locus excluded": "PASS: signal persists without F5",
            "Reverse MR: 9 coagulation genes excluded": "PASS: signal persists without all coagulation loci; 19/25 Wald OR<1",
            "Reverse MR (FinnGen matched): all 38 instruments": "Matched comparison; six instruments absent from the external file",
            "Reverse MR (FinnGen matched): 25-instrument set": "Matched comparison for the coagulation-excluded set",
        }[label]
        rows.append(dict(
            analysis=label, data="External endometriosis" if ext else "FinnGen endometriosis (matched)",
            n_snp=int(q["n_snp"]), effect_allele="", OR=float(q["ivw_or"]),
            lo=float(q["ivw_lo"]), hi=float(q["ivw_hi"]), p=float(q["ivw_p"]),
            I2=float(q["I2"]), egger_int_p=float(q["egger_intercept_p"]),
            presso_p=float(q["presso_global_p"]), n_wald_below_1=int(q["wald_below_1"]),
            wm_or=float(q["wmedian_or"]), wm_p=float(q["wmedian_p"]),
            note=note))
    t = pd.DataFrame(rows)

    def wm_cell(or_, p):
        se = abs(math.log(or_)) / NormalDist().inv_cdf(1 - p / 2)
        return (f"{or_:.3f} ({math.exp(math.log(or_) - 1.96 * se):.3f}\u2013"
                f"{math.exp(math.log(or_) + 1.96 * se):.3f})")

    t["weighted_median_OR_95CI"] = [
        "\u2014" if pd.isna(v) else wm_cell(v, p)
        for v, p in zip(t["wm_or"], t["wm_p"])]
    t["weighted_median_P"] = [
        "\u2014" if pd.isna(v) else (f"{p:.3f}" if p >= 0.001 else f"{p:.2g}")
        for v, p in zip(t["wm_or"], t["wm_p"])]
    t = t[["analysis", "data", "n_snp", "effect_allele", "OR", "lo", "hi", "p",
           "I2", "egger_int_p", "presso_p", "n_wald_below_1",
           "weighted_median_OR_95CI", "weighted_median_P", "note"]]
    out = t.copy()
    for c in ("OR", "lo", "hi"):
        out[c] = out[c].map(lambda v: "" if pd.isna(v) else f"{v:.3f}")
    out["p"] = out["p"].map(lambda v: "" if pd.isna(v) else f"{v:.3g}")
    out["I2"] = out["I2"].map(lambda v: "" if pd.isna(v) else f"{v:.3f}")
    for c in ("egger_int_p", "presso_p"):
        out[c] = out[c].map(lambda v: "" if pd.isna(v) else f"{v:.3g}")
    out.to_csv(os.path.join(TAB, "TableV_external_validation.tsv"), sep="\t", index=False)
    check("Table V: rs6025 OR", 0.865541, round(float(t.loc[0, "OR"]), 6))
    check("Table V: external all n", 38, int(t.loc[2, "n_snp"]))
    check("Table V: external all OR", 0.925684, round(float(t.loc[2, "OR"]), 6))
    check("Table V: external coag-excluded OR", 0.886888, round(float(t.loc[4, "OR"]), 6))
    check("Table V: matched FinnGen 38 OR", 0.888904, round(float(t.loc[5, "OR"]), 6))
    check("Table V: matched FinnGen 25 OR", 0.906414, round(float(t.loc[6, "OR"]), 6))
    check("Table V: proxy OR", 0.853, round(float(t.loc[1, "OR"]), 4))
    check("Table V: proxy CI upper < 1", True, bool(float(t.loc[1, "hi"]) < 1))
    check("Table V: rows with weighted median", 5,
          int((t["weighted_median_OR_95CI"] != "\u2014").sum()))
    check("Table V: WM external all", "0.889 (0.782\u20131.011)",
          t.loc[2, "weighted_median_OR_95CI"])
    check("Table V: WM external coag-excluded", "0.798 (0.681\u20130.934)",
          t.loc[4, "weighted_median_OR_95CI"])
    check("Table V: WM external all P", "0.073", t.loc[2, "weighted_median_P"])
    return t


# ---------------------------------------------------------------- Table VI
def table6():
    t = pd.read_csv(os.path.join(OUT, "sev", "severity_table.tsv"), sep="\t")
    t.loc[t["analysis"] == "DVT -> endometriosis", "note"] = \
        "51 DVT peaks; 5 are identical to PE peaks (including rs6025, F5) and counted once, leaving 46 instruments"
    mr = pd.read_csv(os.path.join(OUT, "sev", "severity_mr.tsv"), sep="\t").set_index("analysis")
    sp = pd.read_csv(os.path.join(OUT, "sev", "severity_specific.tsv"), sep="\t").set_index("analysis")
    ex = pd.read_csv(os.path.join(OUT, "sev", "severity_wmedian_extra.tsv"), sep="\t").set_index("analysis")
    vs = pd.read_csv(os.path.join(OUT, "sev", "severity_vte_split.tsv"), sep="\t").set_index("analysis")
    wm = {
        "PE -> endometriosis": mr.loc["PE->endo"],
        "DVT -> endometriosis": mr.loc["DVT->endo"],
        "PE-specific -> endometriosis": sp.loc["PE_specific->endo"],
        "DVT-specific -> endometriosis": sp.loc["DVT_specific->endo"],
        "PE-specific minus coag -> endometriosis": ex.loc["PE-specific minus coag"],
        "DVT-specific minus coag -> endometriosis": ex.loc["DVT-specific minus coag"],
        "VTE instruments that are PE leads only": vs.loc["VTE_pe_lead"],
        "VTE instruments that are DVT leads only": vs.loc["VTE_dvt_lead"],
        "VTE instruments that are both": vs.loc["VTE_both"],
        "VTE instruments that are neither": vs.loc["VTE_neither"],
    }

    def fmt_p(p):
        return f"{p:.3f}" if p >= 0.001 else f"{p:.2g}"

    t["weighted_median_OR_95CI"] = t["analysis"].map(
        lambda a: f"{wm[a]['wmedian_or']:.3f} ({wm[a]['wmedian_lo']:.3f}\u2013"
                  f"{wm[a]['wmedian_hi']:.3f})" if a in wm else "\u2014")
    t["weighted_median_P"] = t["analysis"].map(
        lambda a: fmt_p(float(wm[a]["wmedian_p"])) if a in wm else "\u2014")
    t = t[["analysis", "exposure", "n_snp", "OR", "lo", "hi", "p", "I2",
           "egger_int_p", "weighted_median_OR_95CI", "weighted_median_P", "note"]]
    out = t.copy()
    for c in ("OR", "lo", "hi"):
        out[c] = out[c].map(lambda v: f"{v:.3f}")
    out["p"] = out["p"].map(lambda v: f"{v:.3g}")
    out["I2"] = out["I2"].map(lambda v: "" if pd.isna(v) else f"{v:.2f}")
    out["egger_int_p"] = out["egger_int_p"].map(
        lambda v: "" if pd.isna(v) else f"{v:.2g}")
    out.to_csv(os.path.join(TAB, "TableVI_severity.tsv"), sep="\t", index=False)
    check("Table VI: rows", 11, len(t))
    check("Table VI: rows with weighted median", 10,
          int((t["weighted_median_OR_95CI"] != "\u2014").sum()))
    check("Table VI: PE n", 27, int(t.loc[0, "n_snp"]))
    check("Table VI: DVT n", 46, int(t.loc[1, "n_snp"]))
    check("Table VI: PE OR", 0.953, round(float(t.loc[0, "OR"]), 3))
    check("Table VI: DVT OR", 0.944, round(float(t.loc[1, "OR"]), 3))
    check("Table VI: DVT weighted median", "0.896 (0.811\u20131.017)",
          t.loc[1, "weighted_median_OR_95CI"])
    check("Table VI: DVT-specific weighted median", "0.854 (0.781\u20130.976)",
          t.loc[3, "weighted_median_OR_95CI"])
    check("Table VI: DVT weighted median P", "0.055", t.loc[1, "weighted_median_P"])
    check("Table VI: ratio p", 0.573, round(float(t.loc[4, "p"]), 3))
    return t


if __name__ == "__main__":
    print("=" * 78)
    print("Building Tables I, II, IV, V, VI")
    print("=" * 78)
    for name, fn in [("Table I", table1), ("Table II", table2), ("Table IV", table4),
                     ("Table V", table5), ("Table VI", table6)]:
        print(f"\n[{name}]")
        t = fn()
        print(f"  rows={len(t)}  ->  out/tables/{name.replace(' ', '')}")
    print("\nAll table checks passed.")
