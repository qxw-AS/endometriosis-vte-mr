#!/usr/bin/env python3
"""M3 stage 1: candidate regions from both traits' clumped lead SNPs.

Protocol M3 stage 1: union of lead SNP +/- 500 kb from both traits; MHC listed
separately and excluded from the primary count. LAVA (M2) was not run in this
environment, so stage 2 falls back to the protocol-defined lead +/- 500 kb
windows.
"""
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out", "m3")
os.makedirs(OUT, exist_ok=True)

WINDOW = 500_000
MHC = ("6", 25_000_000, 35_000_000)

frames = []
for trait in ("endo", "vte"):
    d = pd.read_csv(os.path.join(ROOT, "out", "mr", f"{trait}_leads_hg38.tsv"),
                    sep="\t")
    d["trait"] = trait
    frames.append(d)
leads = pd.concat(frames, ignore_index=True)
leads["chr"] = leads["chr"].astype(str)
leads["start"] = (leads["pos"] - WINDOW).clip(lower=0)
leads["end"] = leads["pos"] + WINDOW

rows = []
for chrom, g in leads.groupby("chr", sort=False):
    g = g.sort_values("start")
    cur_s, cur_e = None, None
    members = []
    for _, r in g.iterrows():
        if cur_s is None:
            cur_s, cur_e, members = r["start"], r["end"], [r]
        elif r["start"] <= cur_e:
            cur_e = max(cur_e, r["end"])
            members.append(r)
        else:
            rows.append((chrom, cur_s, cur_e, members))
            cur_s, cur_e, members = r["start"], r["end"], [r]
    if cur_s is not None:
        rows.append((chrom, cur_s, cur_e, members))

recs = []
for i, (chrom, s, e, members) in enumerate(
        sorted(rows, key=lambda x: (int(x[0]), x[1])), start=1):
    m = pd.DataFrame(members)
    mhc = chrom == MHC[0] and not (e < MHC[1] or s > MHC[2])
    recs.append(dict(
        region=f"R{i:02d}", chr=chrom, start=int(s), end=int(e),
        length_bp=int(e - s),
        n_leads=len(m),
        n_endo=int((m["trait"] == "endo").sum()),
        n_vte=int((m["trait"] == "vte").sum()),
        leads_endo=",".join(m.loc[m["trait"] == "endo", "SNP"].astype(str)),
        leads_vte=",".join(m.loc[m["trait"] == "vte", "SNP"].astype(str)),
        mhc=mhc,
    ))
regions = pd.DataFrame(recs)
regions["lava_status"] = "NOT_RUN"
regions.to_csv(os.path.join(OUT, "candidate_regions.tsv"), sep="\t", index=False)
regions[["chr", "start", "end", "region"]].to_csv(
    os.path.join(OUT, "candidate_regions.bed"), sep="\t", index=False, header=False)

print(f"{len(leads)} lead SNPs -> {len(regions)} merged regions "
      f"({int(regions['mhc'].sum())} MHC), total {regions['length_bp'].sum()/1e6:.1f} Mb")
print(regions[["region", "chr", "start", "end", "n_endo", "n_vte", "mhc"]]
      .to_string(index=False))
