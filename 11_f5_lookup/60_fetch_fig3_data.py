#!/usr/bin/env python3
"""Fetch and assemble the two Fig 3 source tables (V6.2).

Fig 3A : F5 rs6025 (1:169549811:C:T) across 8 FinnGen R13 phenotypes,
         downloaded live from the FinnGen R13 pheweb variant API.
Fig 3B : MAGMA gene-level P values for the 9 coagulation genes,
         read from the same MAGMA runs used for Fig 2B.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "out")
os.makedirs(OUT, exist_ok=True)

API = "https://r13.finngen.fi/api/variant/1-169549811-C-T"
RAW = os.path.join(DATA, "finngen_r13_rs6025_api.json")

# phenotype order follows the protocol (Table 41): thrombotic then gynaecological
PHENOS = [
    ("thrombotic", "I9_PHLETHROMBDVTLOW", "Deep vein thrombosis"),
    ("thrombotic", "I9_THROMBOTH", "Other embolism/thrombosis"),
    ("thrombotic", "I9_VTE", "Venous thromboembolism"),
    ("thrombotic", "I9_PULMEMB", "Pulmonary embolism"),
    ("gynaecological", "N14_OVARYCYST", "Ovarian cyst"),
    ("gynaecological", "N14_FEMALEINFERT", "Female infertility"),
    ("gynaecological", "N14_ENDOMETRIOSIS", "Endometriosis"),
    ("gynaecological", "N14_MENORRHAGIA", "Menorrhagia"),
]


def download_api():
    if os.path.exists(RAW) and os.path.getsize(RAW) > 500_000:
        with open(RAW) as fh:
            d = json.load(fh)
        if len(d.get("results", [])) > 2000:
            print(f"reuse cached API response ({os.path.getsize(RAW)} bytes)")
            return d
    for attempt in range(4):
        try:
            r = requests.get(API, timeout=90)
            r.raise_for_status()
            d = r.json()
            if len(d.get("results", [])) < 2000:
                raise ValueError("unexpected API payload")
            with open(RAW, "w") as fh:
                json.dump(d, fh)
            print(f"downloaded API response ({len(r.content)} bytes)")
            return d
        except Exception as exc:  # pragma: no cover
            print(f"attempt {attempt + 1} failed: {exc}")
            time.sleep(5 * (attempt + 1))
    raise SystemExit("could not download FinnGen R13 API response")


def build_fig3a(d):
    by_code = {x["phenocode"]: x for x in d["results"]}
    rows = []
    for lineage, code, label in PHENOS:
        x = by_code[code]
        beta, se = float(x["beta"]), float(x["sebeta"])
        rows.append(dict(
            lineage=lineage, phenotype=label, phenocode=code,
            n_case=int(x["n_case"]), n_control=int(x["n_control"]),
            beta=beta, se=se,
            OR=float(np.exp(beta)),
            OR_lo=float(np.exp(beta - 1.959963984540054 * se)),
            OR_hi=float(np.exp(beta + 1.959963984540054 * se)),
            p=float(x["pval"]),
            finngen_phenostring=x["phenostring"],
            source="FinnGen R13 pheweb variant API 1:169549811:C:T",
        ))
    df = pd.DataFrame(rows)
    # deterministic display order: protocol lineage (thrombotic first), then OR descending
    df["lineage_order"] = df["lineage"].map({"thrombotic": 0, "gynaecological": 1})
    df = df.sort_values(["lineage_order", "OR"], ascending=[True, False]).reset_index(drop=True)
    return df


def build_fig3b():
    genes = ["F5", "F2", "ABO", "F11", "FGG", "PROCR", "KNG1", "THBD", "SERPINC1"]
    e = pd.read_csv(os.path.join(OUT, "magma", "endo_gene_fixed.genes.out"), sep=r"\s+")
    v = pd.read_csv(os.path.join(OUT, "magma", "vte_gene_fixed.genes.out"), sep=r"\s+")
    m = (e[["GENE", "CHR", "START", "STOP", "NSNPS", "ZSTAT", "P"]]
         .merge(v[["GENE", "ZSTAT", "P"]], on="GENE", suffixes=("_endo", "_vte")))
    m = m[m["GENE"].isin(genes)].copy()
    m["x_endo_mlogp"] = -np.log10(m["P_endo"])
    m["y_vte_mlogp"] = -np.log10(m["P_vte"])
    m = m.sort_values("y_vte_mlogp", ascending=False).reset_index(drop=True)
    missing = sorted(set(genes) - set(m["GENE"]))
    print(f"Fig 3B: {len(m)}/9 coagulation genes with MAGMA P; missing={missing}")
    if len(missing) > 3:
        raise SystemExit("MAGMA gene-level P missing for >3 genes: panel must be a table")
    return m


if __name__ == "__main__":
    d = download_api()
    a = build_fig3a(d)
    b = build_fig3b()
    a.to_csv(os.path.join(OUT, "fig3a_f5_spectrum.tsv"), sep="\t", index=False,
             float_format="%.10g")
    b.to_csv(os.path.join(OUT, "fig3b_coag_genes.tsv"), sep="\t", index=False,
             float_format="%.10g")
    print(a[["phenotype", "phenocode", "n_case", "n_control", "OR", "OR_lo",
             "OR_hi", "p"]].to_string(index=False))
    print(b[["GENE", "CHR", "START", "STOP", "NSNPS", "ZSTAT_endo", "P_endo",
             "ZSTAT_vte", "P_vte", "x_endo_mlogp", "y_vte_mlogp"]].to_string(index=False))
    print("wrote out/fig3a_f5_spectrum.tsv, out/fig3b_coag_genes.tsv")
