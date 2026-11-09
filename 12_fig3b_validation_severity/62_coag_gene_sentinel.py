#!/usr/bin/env python3
"""Fig 3B revision: sentinel-SNP ORs for the 9 coagulation genes.

For each gene, take all SNPs annotated by MAGMA (gene body +/-10 kb, 1000G EUR),
find the sentinel SNP with the strongest evidence across the two traits
(max of -log10 P), then report its OR in endometriosis and in VTE.

This replaces the SNPwise-mean MAGMA scatter (which dilutes single strong
variants such as rs6025) with a directional comparison on a common scale.
"""
import os
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "out")
GENES = ["F5", "F2", "ABO", "F11", "FGG", "PROCR", "KNG1", "THBD", "SERPINC1"]
CHROMS = {"1", "3", "4", "9", "11", "20"}
CACHE = os.path.join(OUT, "fig3b_sentinel_variants.tsv")

# ---------------------------------------------------------------- gene SNP lists
annot = os.path.join(OUT, "magma", "gene_annot_fixed.genes.annot")
gene_snps, gene_region = {}, {}
with open(annot) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if parts[0] in GENES:
            gene_region[parts[0]] = parts[1]
            gene_snps[parts[0]] = [s for s in parts[2:] if s.startswith("rs")]
snp2gene = {s: g for g, ss in gene_snps.items() for s in ss}
print("annotated SNPs:", {g: len(v) for g, v in gene_snps.items()})
print("union SNPs:", len(snp2gene))
assert set(gene_snps) == set(GENES)


def read_trait(path, label):
    t0 = time.time()
    keep = []
    for chunk in pd.read_csv(path, sep="\t", chunksize=1_000_000,
                             dtype={"#chrom": str, "rsids": str},
                             low_memory=False):
        chunk = chunk[chunk["#chrom"].isin(CHROMS)]
        if chunk.empty:
            continue
        mask = [any(s in snp2gene for s in str(v).split(",")) for v in chunk["rsids"]]
        sub = chunk.loc[mask, ["#chrom", "pos", "ref", "alt", "rsids",
                               "pval", "beta", "sebeta", "af_alt"]].copy()
        if not sub.empty:
            keep.append(sub)
    df = pd.concat(keep, ignore_index=True)
    # one row per variant; keep the first rsID that is in our gene map
    df["SNP"] = df["rsids"].apply(lambda s: next(x for x in str(s).split(",") if x in snp2gene))
    df = df.drop_duplicates("SNP")
    df["trait"] = label
    df = df.rename(columns={"pval": f"p_{label}", "beta": f"beta_{label}",
                            "sebeta": f"se_{label}", "af_alt": f"af_{label}"})
    print(f"{label}: {len(df)} variants retained ({time.time()-t0:.0f}s)")
    return df


if os.path.exists(CACHE):
    both = pd.read_csv(CACHE, sep="\t")
    print(f"reuse {CACHE} ({len(both)} variants)")
else:
    e = read_trait(os.path.join(DATA, "finngen_R13_N14_ENDOMETRIOSIS.gz"), "endo")
    v = read_trait(os.path.join(DATA, "finngen_R13_I9_VTE.gz"), "vte")
    both = e.merge(v[["SNP", "p_vte", "beta_vte", "se_vte", "af_vte"]], on="SNP", how="inner")
    both.to_csv(CACHE, sep="\t", index=False, float_format="%.10g")
    print("merged variants:", len(both))

both["gene"] = both["SNP"].map(snp2gene)
both = both[(both["af_endo"].between(0.001, 0.999)) & (both["af_vte"].between(0.001, 0.999))]
both = both[np.isfinite(both["beta_endo"]) & np.isfinite(both["beta_vte"])]

rows = []
for g in GENES:
    sub = both[both["gene"] == g].sort_values(["p_vte", "p_endo"])
    if sub.empty:
        raise SystemExit(f"no variant for {g}")
    r = sub.iloc[0]
    rows.append(dict(
        gene=g, SNP=r["SNP"], chrom=r["#chrom"], pos=int(r["pos"]),
        ref=r["ref"], alt=r["alt"],
        beta_endo=float(r["beta_endo"]), se_endo=float(r["se_endo"]),
        OR_endo=float(np.exp(r["beta_endo"])),
        OR_endo_lo=float(np.exp(r["beta_endo"] - 1.96 * r["se_endo"])),
        OR_endo_hi=float(np.exp(r["beta_endo"] + 1.96 * r["se_endo"])),
        p_endo=float(r["p_endo"]), af_endo=float(r["af_endo"]),
        beta_vte=float(r["beta_vte"]), se_vte=float(r["se_vte"]),
        OR_vte=float(np.exp(r["beta_vte"])),
        OR_vte_lo=float(np.exp(r["beta_vte"] - 1.96 * r["se_vte"])),
        OR_vte_hi=float(np.exp(r["beta_vte"] + 1.96 * r["se_vte"])),
        p_vte=float(r["p_vte"]), af_vte=float(r["af_vte"]),
        n_variants=len(sub), gene_region=gene_region[g],
        sentinel_criterion="lowest VTE P within gene +/-10 kb",
    ))
out = pd.DataFrame(rows)
out["opposite"] = np.sign(out["beta_endo"]) != np.sign(out["beta_vte"])
out = out.sort_values("OR_vte", ascending=True).reset_index(drop=True)
out.to_csv(os.path.join(OUT, "fig3b_sentinel_or.tsv"), sep="\t", index=False, float_format="%.10g")
show = out[["gene", "SNP", "OR_endo", "OR_endo_lo", "OR_endo_hi", "p_endo",
            "OR_vte", "OR_vte_lo", "OR_vte_hi", "p_vte", "opposite", "n_variants"]]
print(show.to_string(index=False))
print("opposite-direction genes:", int(out["opposite"].sum()))
