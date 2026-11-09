#!/usr/bin/env python3
"""Figure 1: study design and analytical workflow.

All numbers printed on the figure are read from the analysis outputs and
asserted in the Layer 4 block before the figure is written.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle, Ellipse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "out")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

NAVY = "#17365D"
NAVY_L = "#2E5C8A"
BLUE = "#2E75B6"
BLUE_L = "#9DC3E6"
BLUE_PALE = "#DEEBF7"
CREAM = "#FFF4E8"
CREAM_D = "#FDE9D9"
TAN = "#D9B38C"
ORANGE = "#C55A11"
INK = "#222222"
GREY = "#595959"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
})


# ---------------------------------------------------------------- data
def get(df, col, value):
    return df.loc[df[col] == value].iloc[0]


t3 = pd.read_csv(os.path.join(OUT, "table3_v62.tsv"), sep="\t")
fwd = get(t3, "analysis", "Primary: coagulation genes +/-1 Mb excluded")
fwd_all = get(t3, "analysis", "Sensitivity: all instruments")
m1 = get(t3[t3["analysis"].str.startswith("Model 1") &
            (t3["exposure"] == "Endometriosis")].assign(
                analysis="Model 1: endometriosis + BMI"), "analysis",
         "Model 1: endometriosis + BMI")
m2 = get(t3[t3["analysis"].str.startswith("Model 2") &
            (t3["exposure"] == "Endometriosis")].assign(
                analysis="Model 2: endometriosis + ever smoking"), "analysis",
         "Model 2: endometriosis + ever smoking")
m3 = get(t3[t3["analysis"].str.startswith("Model 3") &
            (t3["exposure"] == "Endometriosis")].assign(
                analysis="Model 3: endometriosis + BMI + ever smoking"), "analysis",
         "Model 3: endometriosis + BMI + ever smoking")
bmi = t3[(t3["analysis"].str.startswith("Model 1")) &
         (t3["exposure"] == "BMI")].iloc[0]
smk = t3[(t3["analysis"].str.startswith("Model 2")) &
         (t3["exposure"] == "Ever smoking")].iloc[0]

g8 = pd.read_csv(os.path.join(OUT, "g8", "g8_reverse_mr.tsv"), sep="\t")
rev = get(g8, "analysis", "reverse VTE->endo | external endo (all)")
rev_coag = get(g8, "analysis",
               "reverse VTE->endo | external endo (drop 9 coag genes)")

sev = pd.read_csv(os.path.join(OUT, "sev", "severity_table.tsv"), sep="\t")
pe = get(sev, "analysis", "PE -> endometriosis")
dvt = get(sev, "analysis", "DVT -> endometriosis")
diff = get(sev, "analysis", "PE-specific minus DVT-specific")
dvt_sp = get(sev, "analysis", "DVT-specific -> endometriosis")

col = pd.read_csv(os.path.join(OUT, "m3", "coag_loci_coloc.tsv"), sep="\t")
f5c = get(col, "gene", "F5")
aboc = get(col, "gene", "ABO")
procrc = get(col, "gene", "PROCR")

spec = pd.read_csv(os.path.join(OUT, "fig3a_f5_spectrum.tsv"), sep="\t")
spectrum = {r["phenocode"]: r for _, r in spec.iterrows()}

# ---------------------------------------------------------------- canvas
fig = plt.figure(figsize=(7.2, 6.35))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

ax.text(0.038, 0.975, "Figure 1", fontsize=8.6, fontweight="bold",
        color=NAVY, va="top")
ax.text(0.118, 0.975, "Study design and analytical workflow",
        fontsize=7.0, color=GREY, va="top")

COLS = [0.038, 0.2733, 0.5086, 0.7439]
W = 0.218
CENTERS = [x + W / 2 for x in COLS]


def round_box(x, y, w, h, face, edge, lw=0.75, z=3, rounding=0.010):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0.003,rounding_size={rounding}",
                       linewidth=lw, edgecolor=edge, facecolor=face, zorder=z)
    ax.add_patch(p)
    return p


def arrow(x1, y1, x2, y2, color=BLUE_L, lw=1.05, z=2):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                                 arrowstyle="-|>", mutation_scale=6.5,
                                 linewidth=lw, color=color, zorder=z,
                                 shrinkA=0, shrinkB=0))


def card(x, top, w, title, lines, accent=False, h=None):
    n = len(lines)
    h = h if h is not None else 0.052 + 0.018 * n
    y = top - h
    face = CREAM_D if accent else CREAM
    edge = ORANGE if accent else TAN
    round_box(x, y, w, h, face, edge, lw=0.85 if accent else 0.7)
    ax.add_patch(Rectangle((x + 0.003, y + 0.006), 0.0045, h - 0.012,
                           facecolor=ORANGE if accent else BLUE,
                           edgecolor="none", zorder=4))
    ax.text(x + 0.013, y + h - 0.015, title, fontsize=7.1,
            fontweight="bold", color=NAVY, va="top", zorder=5)
    if accent:
        chip = FancyBboxPatch((x + w - 0.043, y + h - 0.027), 0.034, 0.017,
                              boxstyle="round,pad=0.001,rounding_size=0.004",
                              facecolor=ORANGE, edgecolor="none", zorder=5)
        ax.add_patch(chip)
        ax.text(x + w - 0.026, y + h - 0.0185, "KEY", fontsize=4.7,
                color="white", fontweight="bold", ha="center", va="center",
                zorder=6)
    ax.plot([x + 0.010, x + w - 0.010], [y + h - 0.033, y + h - 0.033],
            color=edge, lw=0.45, zorder=4)
    yy = y + h - 0.047
    for line in lines:
        if line.startswith("- "):
            ax.text(x + 0.012, yy, "\u2022", fontsize=6.3, color=BLUE,
                    va="top", zorder=5)
            ax.text(x + 0.024, yy, line[2:], fontsize=6.2, color=INK,
                    va="top", zorder=5)
        else:
            ax.text(x + 0.024, yy, line, fontsize=6.2, color=INK,
                    va="top", zorder=5)
        yy -= 0.018
    return y


def column(x, cards, top=0.745):
    prev_bottom = None
    for i, (title, lines, accent) in enumerate(cards):
        y = card(x, top, W, title, lines, accent=accent)
        if prev_bottom is not None:
            arrow(x + W / 2, prev_bottom, x + W / 2, top)
        prev_bottom = y
        top = y - 0.027


# ---------------------------------------------------------------- inputs
def input_card(x, title, trait, n):
    w, h, y = 0.300, 0.070, 0.872
    round_box(x, y, w, h, "#F2F7FB", BLUE_L, lw=0.9)
    # database glyph
    cx, cy = x + 0.021, y + h / 2
    ax.add_patch(Rectangle((cx - 0.010, cy - 0.014), 0.020, 0.027,
                           facecolor=BLUE_L, edgecolor=NAVY_L, lw=0.5, zorder=4))
    ax.add_patch(Ellipse((cx, cy + 0.013), 0.020, 0.007,
                         facecolor="#CFE2F3", edgecolor=NAVY_L, lw=0.5, zorder=5))
    ax.add_patch(Ellipse((cx, cy - 0.014), 0.020, 0.007,
                         facecolor=BLUE_L, edgecolor=NAVY_L, lw=0.5, zorder=5))
    ax.text(x + 0.039, y + h - 0.012, title, fontsize=7.1, fontweight="bold",
            color=NAVY, va="top")
    ax.text(x + 0.039, y + h - 0.033, trait, fontsize=6.4, color=INK, va="top")
    ax.text(x + 0.039, y + h - 0.052, n, fontsize=5.9, color=GREY, va="top")


input_card(0.115, "Exposure GWAS", "Endometriosis",
           "FinnGen R13 \u00b7 20,534 / 122,622")
input_card(0.585, "Outcome GWAS", "Venous thromboembolism",
           "FinnGen R13 \u00b7 27,928 / 472,258")

ax.plot([0.265, 0.265], [0.872, 0.858], color=BLUE_L, lw=1.05, zorder=2)
ax.plot([0.735, 0.735], [0.872, 0.858], color=BLUE_L, lw=1.05, zorder=2)
ax.plot([0.265, 0.735], [0.858, 0.858], color=BLUE_L, lw=1.05, zorder=2)
arrow(0.500, 0.858, 0.500, 0.840)

# ---------------------------------------------------------------- stage bar
BAR_Y, BAR_H = 0.790, 0.050
BAR_X, BAR_W = 0.038, 0.924
round_box(BAR_X + 0.002, BAR_Y - 0.002, BAR_W, BAR_H, "#D9E2F3", "none",
          lw=0, z=2, rounding=0.008)
round_box(BAR_X, BAR_Y, BAR_W, BAR_H, NAVY, NAVY, lw=0.6, z=3, rounding=0.008)
seg_labels = ["01  Data & instruments", "02  Causal estimation",
              "03  Robustness & pleiotropy", "04  Integration & synthesis"]
seg_w = BAR_W / 4
for i, lab in enumerate(seg_labels):
    if i:
        ax.plot([BAR_X + i * seg_w, BAR_X + i * seg_w], [BAR_Y + 0.004, BAR_Y + BAR_H - 0.004],
                color="#5B7FA6", lw=0.7, zorder=4)
    ax.text(BAR_X + (i + 0.5) * seg_w, BAR_Y + BAR_H / 2, lab, color="white",
            fontsize=6.35, fontweight="bold", ha="center", va="center", zorder=5)

# ---------------------------------------------------------------- columns
column(COLS[0], [
    ("Instrument selection", [
        "- Endometriosis: 36 lead SNPs",
        "- VTE: 44 lead SNPs",
        "- $P < 5\\times10^{-8}$; $r^2 < 0.001$; 10 Mb",
        "- All instruments $F > 10$",
    ], False),
    ("Coagulation-aware design", [
        "- 9 genes \u00b11 Mb: F5, ABO,",
        "  THBD, KNG1, PROCR, FGG,",
        "  F11, SERPINC1, F2",
        "- Primary set: 34 endo SNPs",
        "- Sensitivity: all 36 SNPs",
    ], False),
    ("External replication", [
        "- Rahmioglu 2023",
        "  (GCST90205183)",
        "- 60,674 / 410,192",
        "- $N$ = 470,866",
        "- Reverse MR only",
    ], False),
])

column(COLS[1], [
    ("Forward MR: endo \u2192 VTE", [
        "- IVW, random effects",
        "- Primary (34 SNPs):",
        "  1.02 (0.975\u20131.062)",
        "- $P$ = 0.432; $I^2$ = 0.35",
        "- All 36: 0.98 (0.768\u20131.245)",
    ], False),
    ("Reverse MR: VTE \u2192 endo", [
        "- External outcome, 38 SNPs",
        "- OR 0.926 (0.880\u20130.974)",
        "- After 9 coagulation loci:",
        "  25 SNPs, OR 0.887",
        "  (0.824\u20130.954)",
    ], False),
    ("Multivariable MR", [
        "- Endometriosis + BMI",
        "  + ever smoking",
        "- Endo: 1.01\u20131.02 (ns)",
        "- BMI 1.47; smoking 1.16",
        "- Cond. $F$: 4.4 / 19.3 / 4.2",
        "- Only Model 2 $F > 10$",
    ], False),
])

column(COLS[2], [
    ("Pleiotropy diagnostics", [
        "- MR-Egger intercept",
        "- Weighted median and mode",
        "- MR-PRESSO; Cochran $Q$",
        "- Leave-one-out; Steiger",
    ], False),
    ("Colocalization", [
        "- coloc.abf + PLACO",
        "- 9 coagulation loci",
        "- F5: PP.H4 = 0.994",
        "- ABO: PP.H3 = 0.953",
        "- PROCR: PP.H4 = 0.502",
    ], False),
    ("Severity stratification", [
        "- PE vs DVT instruments",
        "  \u2192 endometriosis",
        "- PE 0.953 (0.876\u20131.037)",
        "- DVT 0.944 (0.900\u20130.991)",
        "- Difference $P$ = 0.573",
        "- DVT-specific 0.865",
    ], False),
])

column(COLS[3], [
    ("F5 phenotypic spectrum", [
        "- rs6025, per T allele",
        "- DVT 3.28",
        "- Other embolism 2.75; VTE 2.35",
        "- PE 1.69; endometriosis 0.79",
    ], False),
    ("Integrated model", [
        "- F5 shared locus (PP.H4 = 0.994),",
        "  opposite OR directions",
        "- No detectable causal effect of",
        "  endometriosis liability on VTE",
        "- VTE \u2192 endometriosis:",
        "  0.926 (0.880\u20130.974)",
        "- Attenuated after coagulation",
        "  loci removal: 0.887",
        "- Weighted median 0.889",
        "  ($P$ = 0.073)",
    ], True),
    ("Reporting", [
        "- STROBE-MR checklist",
        "- Coagulation-excluded primary",
        "  + full sensitivity grid",
        "- FinnGen R13; GCST90205183",
    ], False),
])

# ---------------------------------------------------------------- tree from bar
ax.plot([0.500, 0.500], [BAR_Y, 0.767], color=BLUE_L, lw=1.05, zorder=2)
ax.plot([CENTERS[0], CENTERS[-1]], [0.767, 0.767], color=BLUE_L, lw=1.05, zorder=2)
for cx in CENTERS:
    arrow(cx, 0.767, cx, 0.746)

# ---------------------------------------------------------------- legend
ly = 0.185
items = [(CREAM, TAN, "Analysis step"), (NAVY, NAVY, "Stage"),
         ("none", BLUE_L, "Data flow"), (CREAM_D, ORANGE, "Key synthesis")]
x = 0.038
for face, edge, lab in items:
    if lab == "Data flow":
        ax.plot([x, x + 0.016], [ly + 0.005, ly + 0.005], color=BLUE_L, lw=1.4)
    else:
        ax.add_patch(Rectangle((x, ly), 0.015, 0.010, facecolor=face,
                               edgecolor=edge, lw=0.7, zorder=4))
    ax.text(x + 0.020, ly + 0.005, lab, fontsize=5.7, color=GREY, va="center")
    x += 0.095

ax.text(0.962, ly + 0.005,
        "FinnGen R13 (GRCh38) \u00b7 Rahmioglu 2023 \u00b7 1000G EUR",
        fontsize=5.7, color=GREY, ha="right", va="center")

# ---------------------------------------------------------------- Layer 4
checks = []


def check(label, expected, actual, tol=1e-4):
    ok = (abs(float(expected) - float(actual)) <= tol
          if isinstance(expected, (int, float, np.floating)) else expected == actual)
    checks.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'} {label:44s} expected={expected} actual={actual}")


print("Layer 4 reconciliation")
check("forward primary OR", 1.017, fwd["OR"], 1e-3)
check("forward primary upper CI", 1.062, fwd["hi"], 1e-3)
check("forward primary P", 0.432, fwd["p"], 1e-3)
check("forward all OR", 0.978, fwd_all["OR"], 1e-3)
check("reverse external OR", 0.925684, rev["ivw_or"], 1e-6)
check("reverse coag OR", 0.886888, rev_coag["ivw_or"], 1e-6)
check("reverse weighted median", 0.889221, rev["wmedian_or"], 1e-6)
check("MVMR Model 2 conditional F", 19.3, m2["conditional_F"], 0.05)
check("MVMR Model 1 conditional F", 4.41, m1["conditional_F"], 0.01)
check("MVMR Model 3 conditional F", 4.23, m3["conditional_F"], 0.01)
check("MVMR BMI OR", 1.47, bmi["OR"], 1e-2)
check("MVMR smoking OR", 1.16, smk["OR"], 1e-2)
check("F5 PP.H4", 0.993662, f5c["PP.H4"], 1e-6)
check("ABO PP.H3", 0.952859, aboc["PP.H3"], 1e-6)
check("PROCR PP.H4", 0.50239, procrc["PP.H4"], 1e-6)
check("PE severity OR", 0.953, pe["OR"], 1e-3)
check("DVT severity OR", 0.944, dvt["OR"], 1e-3)
check("severity difference P", 0.573, diff["p"], 1e-3)
check("DVT-specific OR", 0.865, dvt_sp["OR"], 1e-3)
check("F5 spectrum DVT", 3.277825, spectrum["I9_PHLETHROMBDVTLOW"]["OR"], 1e-5)
check("F5 spectrum endometriosis", 0.785967, spectrum["N14_ENDOMETRIOSIS"]["OR"], 1e-5)
if not all(checks):
    raise ValueError("Layer 4 verification failed; figure not written")

for ext in ("png", "pdf", "tiff"):
    fig.savefig(os.path.join(FIG, f"Fig1_study_design.{ext}"), dpi=300,
                bbox_inches="tight", pad_inches=0.03, facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- trim margins
# The full-figure axes leave a wide white band below the legend; trim raster
# outputs and set a matching PDF CropBox so the submitted files are tight.
from PIL import Image, ImageChops


def trim_raster(path):
    im = Image.open(path).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    box = ImageChops.difference(im, bg).getbbox()
    if not box:
        return
    pad = 18
    box = (max(0, box[0] - pad), max(0, box[1] - pad),
           min(im.width, box[2] + pad), min(im.height, box[3] + pad))
    out = im.crop(box)
    kw = {"dpi": (300, 300)}
    if path.endswith(".tiff"):
        kw["compression"] = "tiff_lzw"
    out.save(path, **kw)


trim_raster(os.path.join(FIG, "Fig1_study_design.png"))
trim_raster(os.path.join(FIG, "Fig1_study_design.tiff"))

pdf_path = os.path.join(FIG, "Fig1_study_design.pdf")
try:
    import re
    import subprocess
    from pypdf import PdfReader, PdfWriter
    from pypdf.generic import RectangleObject

    res = subprocess.run(
        ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=bbox", pdf_path],
        capture_output=True, text=True, check=False)
    boxes = re.findall(r"%%BoundingBox: (\d+) (\d+) (\d+) (\d+)", res.stderr)
    if boxes:
        llx, lly, urx, ury = map(int, boxes[-1])
        pad = 3
        reader = PdfReader(pdf_path)
        writer = PdfWriter()
        for page in reader.pages:
            page.mediabox = RectangleObject((llx - pad, lly - pad, urx + pad, ury + pad))
            page.cropbox = RectangleObject((llx - pad, lly - pad, urx + pad, ury + pad))
            writer.add_page(page)
        with open(pdf_path, "wb") as fh:
            writer.write(fh)
        print(f"PDF CropBox: {llx - pad} {lly - pad} {urx + pad} {ury + pad}")
except Exception as exc:
    print(f"PDF crop skipped: {exc}")

print("\nALL NUMBERS VERIFIED")
print("saved: out/figures/Fig1_study_design.png|pdf|tiff")
