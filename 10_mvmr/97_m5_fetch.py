#!/usr/bin/env python3
"""M5 step 2: fetch FinnGen R13 effects for MVMR instruments.

Uses the per-phenotype pheweb endpoint because the full /api/variant/<id>
payload returns beta=None for many phenotype-variant pairs (verified for
BMI_IRN), while /api/variant/<id>/<pheno> returns the correct beta.

SE is reconstructed from beta and the Wald p-value (or mlogp when p
underflows), which reproduces the stored sebeta to <1% for checked SNPs.
"""
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import log, sqrt

import numpy as np
import pandas as pd
import requests
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m5")
CACHE = os.path.join(OUT, "api_cache")
os.makedirs(CACHE, exist_ok=True)

PHENOS = ["N14_ENDOMETRIOSIS", "I9_VTE", "BMI_IRN", "EVER_SMOKER",
          "I9_PULMEMB", "I9_PHLETHROMBDVTLOW"]
BASE = "https://r13.finngen.fi/api/variant"


def get_one(vid, pheno):
    fn = os.path.join(CACHE, f"{vid}__{pheno}.json")
    if os.path.exists(fn) and os.path.getsize(fn) > 0:
        try:
            return vid, pheno, json.load(open(fn))
        except Exception:
            pass
    url = f"{BASE}/{vid}/{pheno}"
    for attempt in range(4):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 200:
                j = r.json()
                json.dump(j, open(fn, "w"))
                return vid, pheno, j
            if r.status_code == 404:
                json.dump({"results": None}, open(fn, "w"))
                return vid, pheno, {"results": None}
        except Exception:
            time.sleep(0.5 * (attempt + 1))
    return vid, pheno, {"results": None}


def parse(j):
    r = (j or {}).get("results")
    if not r:
        return None
    beta = r.get("beta")
    p = r.get("pval")
    mlogp = r.get("mlogp")
    if beta is None or p is None:
        return None
    beta = float(beta)
    p = float(p)
    if p > 0:
        z = stats.norm.isf(p / 2)
    elif mlogp:
        z = sqrt(2 * float(mlogp) * log(10))
    else:
        return None
    if not np.isfinite(z) or z <= 0 or beta == 0:
        return None
    return dict(beta=beta, se=abs(beta) / z, p=p, maf=r.get("maf"),
                maf_case=r.get("maf_case"), maf_control=r.get("maf_control"),
                n_case=r.get("n_case"), n_control=r.get("n_control"))


def main():
    u = pd.read_csv(os.path.join(OUT, "union_all.tsv"), sep="\t")
    u["chr"] = u["chr"].astype(int)
    u["pos"] = u["pos"].astype(int)
    u["vid"] = (u["chr"].astype(str) + "-" + u["pos"].astype(str) + "-" +
                u["ref"].astype(str) + "-" + u["alt"].astype(str))
    u["group"] = "MVMR_" + u["source"].astype(str)

    rv = pd.read_csv(os.path.join(ROOT, "out", "mr", "harmonised_reverse.tsv"),
                     sep="\t")
    rv = rv[rv["use"]].drop_duplicates("SNP").copy()
    rv["chr"] = rv["chr"].astype(int)
    rv["pos"] = rv["pos"].astype(int)
    rv["vid"] = (rv["chr"].astype(str) + "-" + rv["pos"].astype(str) + "-" +
                 rv["ea"].astype(str) + "-" + rv["nea"].astype(str))
    rv["group"] = "VTE"
    rv = rv.rename(columns={"ea": "alt", "nea": "ref", "b_exp": "beta",
                            "se_exp": "sebeta"})
    rv = rv[["SNP", "chr", "pos", "ref", "alt", "beta", "sebeta", "vid",
             "group"]]
    allv = pd.concat([u[["SNP", "chr", "pos", "ref", "alt", "beta", "sebeta",
                         "vid", "group"]], rv], ignore_index=True)
    # a variant can be both an MVMR instrument and a VTE instrument
    allv = allv.drop_duplicates("vid").reset_index(drop=True)

    # probe phenotype per group: use the trait the variant was selected for
    probe = {"MVMR_BMI": "BMI_IRN", "MVMR_SMOKE": "EVER_SMOKER",
             "MVMR_ENDO": "N14_ENDOMETRIOSIS", "VTE": "I9_VTE"}

    def alt_of(v):
        p = v.split("-")
        return "-".join(p[:2] + [p[3], p[2]])

    # first pass: resolve allele order with the group-specific probe
    chosen = {}
    probe_tasks = []
    for _, r in allv.iterrows():
        v = r["vid"]
        ph = probe[r["group"]]
        probe_tasks.append((v, ph))
        probe_tasks.append((alt_of(v), ph))
    print(f"variants={len(allv)}, probe queries={len(probe_tasks)}")
    pres = {}
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs = [ex.submit(get_one, v, ph) for v, ph in probe_tasks]
        for f in as_completed(futs):
            vid, ph, j = f.result()
            pres[(vid, ph)] = j
    for _, r in allv.iterrows():
        v = r["vid"]
        ph = probe[r["group"]]
        if (pres.get((v, ph)) or {}).get("results"):
            chosen[v] = (v, False)
        elif (pres.get((alt_of(v), ph)) or {}).get("results"):
            chosen[v] = (alt_of(v), True)
    print(f"resolved variants={len(chosen)}/{len(allv)}")

    # second pass: remaining phenotypes for the resolved variant ID
    rest = [ph for ph in PHENOS]
    tasks = []
    for v, (cand, _) in chosen.items():
        for ph in rest:
            tasks.append((cand, ph))
    print(f"remaining queries={len(tasks)}")
    res = dict(pres)
    done = 0
    with ThreadPoolExecutor(max_workers=10) as ex:
        futs = [ex.submit(get_one, v, ph) for v, ph in tasks]
        for f in as_completed(futs):
            vid, ph, j = f.result()
            res[(vid, ph)] = j
            done += 1
            if done % 500 == 0:
                print(f"  {done}/{len(tasks)}", flush=True)

    rows = []
    for _, r in allv.iterrows():
        if r["vid"] not in chosen:
            continue
        chosen_vid, flip = chosen[r["vid"]]
        for ph in PHENOS:
            pr = parse(res.get((chosen_vid, ph)))
            if pr is None:
                continue
            rows.append(dict(SNP=r["SNP"], vid=chosen_vid, group=r["group"],
                             pheno=ph, beta=pr["beta"] * (-1 if flip else 1),
                             se=pr["se"], p=pr["p"], maf=pr["maf"],
                             n_case=pr["n_case"], n_control=pr["n_control"],
                             flipped=flip))
    long = pd.DataFrame(rows)
    long.to_csv(os.path.join(OUT, "instrument_effects_long.tsv"), sep="\t",
                index=False, float_format="%.6g")
    wide = long.pivot_table(index=["SNP", "vid", "group"], columns="pheno",
                            values=["beta", "se", "p"], aggfunc="first")
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    wide.to_csv(os.path.join(OUT, "instrument_matrix.tsv"), sep="\t",
                index=False, float_format="%.6g")
    print("long rows", len(long), "wide rows", len(wide))
    print(wide.groupby("group").size().to_string())
    print("\ncoverage (non-missing beta per phenotype):")
    print(long.groupby("pheno").size().to_string())


if __name__ == "__main__":
    main()
