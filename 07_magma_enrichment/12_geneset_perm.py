#!/usr/bin/env python3
"""Gene-level competitive enrichment for coagulation / inflammation / hormone sets.

Data: FinnGen R13 GWAS p-values on HapMap3 SNPs (GRCh38), Ensembl GRCh38 gene
locations. Gene statistic: min SNP P with Bonferroni correction over the number
of HapMap3 SNPs assigned to the gene (gene body +/- 10 kb).
Competitive test: 10,000 label permutations stratified by SNP-count decile,
seed 20261108. Real data only; no simulation of GWAS values.
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm
from intervaltree import IntervalTree

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DATA_DIR, OUT_DIR

SEED = 20261108
NPERM = 10000
DATA = str(DATA_DIR)
OUT = str(OUT_DIR)

genes = pd.read_csv(f"{DATA}/magma_gene_loc.txt", sep="\t", header=None,
                    names=["gene", "chr", "start", "end"])
genes = genes.drop_duplicates("gene")
snps = pd.read_csv(f"{DATA}/hm3_snp_p.tsv.gz", sep="\t")
print("genes", len(genes), "snps", len(snps))

# assign SNPs to genes (gene body +/- 10 kb)
recs = []
for chrom, gsub in genes.groupby("chr"):
    ssub = snps[snps.CHR == chrom]
    if ssub.empty:
        continue
    tree = IntervalTree()
    for r in gsub.itertuples():
        tree.addi(r.start - 10000, r.end + 10000 + 1, r.gene)
    for r in ssub.itertuples():
        hits = tree.at(r.BP)
        for h in hits:
            recs.append((h.data, r.SNP, r.P_endo, r.P_vte))
assign = pd.DataFrame(recs, columns=["gene", "SNP", "P_endo", "P_vte"])
print("assignments", len(assign), "genes with SNP", assign.gene.nunique())

gene_stats = assign.groupby("gene").agg(
    n_snp=("SNP", "nunique"),
    min_p_endo=("P_endo", "min"),
    min_p_vte=("P_vte", "min"),
).reset_index()
# Bonferroni over SNPs in the gene (conservative, LD-agnostic)
gene_stats["p_gene_endo"] = np.minimum(1.0, gene_stats.min_p_endo * gene_stats.n_snp)
gene_stats["p_gene_vte"] = np.minimum(1.0, gene_stats.min_p_vte * gene_stats.n_snp)
for t in ["endo", "vte"]:
    # bounded one-sided association score: -log10(gene P)
    gene_stats[f"z_{t}"] = -np.log10(gene_stats[f"p_gene_{t}"].clip(lower=1e-300))
gene_stats.to_csv(f"{OUT}/gene_stats.tsv.gz", sep="\t", index=False)

# exclude MHC 25-34 Mb (chr6) from background and gene sets
mhc = genes[(genes.chr == 6) & (genes.end >= 25_000_000) & (genes.start <= 34_000_000)].gene
gene_stats = gene_stats[~gene_stats.gene.isin(mhc)].copy()
gene_stats["decile"] = pd.qcut(gene_stats.n_snp, 10, labels=False, duplicates="drop")
print("background genes after MHC exclusion", len(gene_stats), "MHC genes", len(mhc))

sets = {}
for line in open(f"{DATA}/magma_gene_sets.txt"):
    p = line.rstrip("\n").split("\t")
    sets[p[0]] = p[1:]

rng = np.random.default_rng(SEED)
bg = gene_stats[["gene", "z_endo", "z_vte", "decile"]].reset_index(drop=True)
pool = {d: bg[bg.decile == d].index.to_numpy() for d in sorted(bg.decile.unique())}

rows = []
for name, glist in sets.items():
    sub = bg[bg.gene.isin(glist)].reset_index(drop=True)
    if len(sub) < 2:
        continue
    for trait in ["endo", "vte"]:
        obs = sub[f"z_{trait}"].mean()
        # stratified permutation: draw one gene per observed gene, same SNP decile
        null = np.empty(NPERM)
        for b in range(NPERM):
            idx = np.concatenate([rng.choice(pool[d], size=(sub.decile == d).sum(), replace=True)
                                  for d in sorted(sub.decile.unique())])
            null[b] = bg.loc[idx, f"z_{trait}"].mean()
        p_emp = (1 + np.sum(null >= obs)) / (1 + NPERM)
        z_emp = (obs - null.mean()) / (null.std(ddof=1) + 1e-12)
        rows.append(dict(set=name, trait=trait, n_genes=len(sub), mean_z=obs,
                         null_mean=null.mean(), null_sd=null.std(ddof=1),
                         enrichment_z=z_emp, p_perm=p_emp,
                         n_perm=NPERM, seed=SEED))
res = pd.DataFrame(rows)
res.to_csv(f"{OUT}/geneset_enrichment.tsv", sep="\t", index=False)
print(res.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
