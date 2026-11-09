#!/usr/bin/env python3
"""M4 step 2: pair exposure/outcome statistics for the clumped instruments.

Both traits are FinnGen R13, so alleles are on the same reference and the main
task is (i) confirming the same variant, (ii) allele alignment, (iii) palindromic
rule per M0 (drop ambiguous palindromes with MAF > 0.42, otherwise align by EAF).

Sample sizes (FinnGen R13 pheweb endpoint metadata):
  N14_ENDOMETRIOSIS  20,534 cases / 122,622 controls
  I9_VTE             27,928 cases / 472,258 controls
"""
import os
import re

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "out", "mr")

FILES = {
    "endo": os.path.join(DATA, "finngen_R13_N14_ENDOMETRIOSIS.gz"),
    "vte": os.path.join(DATA, "finngen_R13_I9_VTE.gz"),
}
N = {"endo": 20534 + 122622, "vte": 27928 + 472258}
CASE = {"endo": 20534, "vte": 27928}


def read_index(trait):
    path = os.path.join(OUT, f"{trait}_clump.clumped")
    df = pd.read_csv(path, sep=r"\s+")
    df = df[["SNP", "P"]].dropna()
    df["SNP"] = df["SNP"].astype(str)
    return df["SNP"].tolist()


def fetch(path, targets):
    """Stream the GWAS once, returning rows whose rsids field contains a target."""
    pat = re.compile(r"(?:^|,)(?:" + "|".join(map(re.escape, targets)) + r")(?:,|$)")
    keep = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=2_000_000,
                             dtype={"#chrom": str, "rsids": str}, low_memory=False):
        m = chunk["rsids"].fillna("").map(lambda s: bool(pat.search(s)))
        if m.any():
            keep.append(chunk[m])
    df = pd.concat(keep, ignore_index=True)
    df = df.rename(columns={"#chrom": "chr", "pval": "p", "sebeta": "se", "af_alt": "af"})
    df["rsids"] = df["rsids"].astype(str)
    df["SNP"] = df["rsids"].str.split(",").str[0]
    return df


def harmonise(exp, out, exp_name, out_name):
    e = exp.copy()
    o = out.copy()
    o = o.rename(columns={"ref": "o_ref", "alt": "o_alt", "af": "o_af", "beta": "b_out",
                          "se": "se_out", "p": "p_out", "pos": "o_pos"})
    m = e.merge(o[["SNP", "o_ref", "o_alt", "o_af", "b_out", "se_out", "p_out"]],
                on="SNP", how="inner", validate="many_to_one")
    rows = []
    for _, r in m.iterrows():
        ea, nea, eaf = r["alt"].upper(), r["ref"].upper(), float(r["af"])
        oa, ona = r["o_alt"].upper(), r["o_ref"].upper()
        b = float(r["b_out"])
        pat = {ea, nea}
        palindromic = pat in ({"A", "T"}, {"C", "G"})
        if (ea, nea) == (oa, ona):
            action = "aligned"
        elif (ea, nea) == (ona, oa):
            b = -b
            action = "swapped"
        elif (ea, nea) == ({"A": "T", "T": "A", "C": "G", "G": "C"}[oa],
                           {"A": "T", "T": "A", "C": "G", "G": "C"}[ona]):
            b = -b
            action = "strand_flip"
        else:
            action = "allele_mismatch"
        if palindromic:
            maf = min(eaf, 1 - eaf)
            if maf > 0.42:
                action = "ambiguous_palindrome_drop"
            elif abs(eaf - float(r["o_af"])) > 0.2 and not abs(eaf + float(r["o_af"]) - 1) > 0.2:
                action = action + "_eaf_warn"
        rows.append(dict(SNP=r["SNP"], chr=r["chr"], pos=r["pos"],
                         ea=ea, nea=nea, eaf=eaf,
                         b_exp=float(r["beta"]), se_exp=float(r["se"]), p_exp=float(r["p"]),
                         b_out=b, se_out=float(r["se_out"]), p_out=float(r["p_out"]),
                         palindromic=palindromic, action=action,
                         exposure=exp_name, outcome=out_name))
    h = pd.DataFrame(rows)
    h["f_stat"] = (h["b_exp"] / h["se_exp"]) ** 2
    h["r2"] = h["b_exp"] ** 2 / (h["b_exp"] ** 2 + h["se_exp"] ** 2 * N[exp_name])
    good = h["action"].isin(["aligned", "swapped", "strand_flip"])
    h["use"] = good & (h["f_stat"] > 10)
    return h


endo_idx = read_index("endo")
vte_idx = read_index("vte")
targets = sorted(set(endo_idx) | set(vte_idx))
print(f"index SNPs: endo {len(endo_idx)}, vte {len(vte_idx)}, union {len(targets)}")

stats = {}
for trait in ("endo", "vte"):
    d = fetch(FILES[trait], targets)
    d = d.drop_duplicates("SNP")
    stats[trait] = d
    print(f"{trait}: matched {len(d)} of {len(targets)} requested variants")

# exposure tables come from the significant-SNP extract (has exact beta/se/ref/alt/af)
sig = {t: pd.read_csv(os.path.join(OUT, f"{t}_sig.tsv"), sep="\t").drop_duplicates("SNP")
       for t in ("endo", "vte")}

forward = harmonise(sig["endo"][sig["endo"]["SNP"].isin(endo_idx)],
                    stats["vte"], "endo", "vte")
reverse = harmonise(sig["vte"][sig["vte"]["SNP"].isin(vte_idx)],
                    stats["endo"], "vte", "endo")

for name, h in (("forward", forward), ("reverse", reverse)):
    h.to_csv(os.path.join(OUT, f"harmonised_{name}.tsv"), sep="\t", index=False,
             float_format="%.10g")
    print(f"\n{name}: {len(h)} pairs, usable {int(h['use'].sum())}")
    print(h["action"].value_counts().to_string())
    print(f"  mean F (usable) = {h.loc[h['use'], 'f_stat'].mean():.1f}; "
          f"min F = {h.loc[h['use'], 'f_stat'].min():.1f}")
print("\nSaved harmonised_forward.tsv / harmonised_reverse.tsv")
