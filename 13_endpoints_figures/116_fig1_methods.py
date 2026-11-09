#!/usr/bin/env python3
"""Figure 1 (methods-only version): study design and analytical workflow.

The figure labels methods and data types only. It deliberately contains no
sample sizes, thresholds, effect estimates, instrument counts, or other
numbers; Layer 4 asserts that no digit appears in any displayed string.
"""
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "out", "figures")
os.makedirs(FIG, exist_ok=True)

NAVY = "#17365D"
NAVY_L = "#2E5C8A"
BLUE = "#2E75B6"
BLUE_L = "#9DC3E6"
CREAM = "#FFF4E8"
CREAM_D = "#FDE9D9"
TAN = "#D9B38C"
ORANGE = "#C55A11"
INK = "#222222"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
})

fig = plt.figure(figsize=(7.2, 6.0))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

COLS = [0.038, 0.2733, 0.5086, 0.7439]
W = 0.218
CENTERS = [x + W / 2 for x in COLS]
TEXTS = []


def say(x, y, s, size=6.2, color=INK, weight="normal", ha="left", va="top", z=5):
    TEXTS.append(s)
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
            ha=ha, va=va, zorder=z)


def round_box(x, y, w, h, face, edge, lw=0.75, z=3, rounding=0.010):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0.003,rounding_size={rounding}",
        linewidth=lw, edgecolor=edge, facecolor=face, zorder=z))


def arrow(x1, y1, x2, y2, color=BLUE_L, lw=1.05, z=2):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=6.5,
        linewidth=lw, color=color, zorder=z, shrinkA=0, shrinkB=0))


def card(x, top, title, lines, accent=False):
    h = 0.112 if len(lines) == 2 else 0.132
    y = top - h
    face = CREAM_D if accent else CREAM
    edge = ORANGE if accent else TAN
    round_box(x, y, W, h, face, edge, lw=0.85 if accent else 0.7)
    ax.add_patch(Rectangle((x + 0.003, y + 0.006), 0.0045, h - 0.012,
                           facecolor=ORANGE if accent else BLUE,
                           edgecolor="none", zorder=4))
    say(x + 0.013, y + h - 0.015, title, size=7.1, color=NAVY,
        weight="bold")
    ax.plot([x + 0.010, x + W - 0.010], [y + h - 0.033, y + h - 0.033],
            color=edge, lw=0.45, zorder=4)
    yy = y + h - 0.047
    for line in lines:
        ax.text(x + 0.012, yy, "\u2022", fontsize=6.3, color=BLUE,
                va="top", zorder=5)
        say(x + 0.024, yy, line, size=6.2)
        yy -= 0.020
    return y


def column(x, cards, top=0.735):
    prev = None
    for title, lines, accent in cards:
        y = card(x, top, title, lines, accent=accent)
        if prev is not None:
            arrow(x + W / 2, prev, x + W / 2, top)
        prev = y
        top = y - (0.058 if len(lines) == 2 else 0.038)


# ---------------------------------------------------------------- inputs
def input_card(x, title, trait):
    w, h, y = 0.300, 0.075, 0.862
    round_box(x, y, w, h, "#F2F7FB", BLUE_L, lw=0.9)
    cx, cy = x + 0.021, y + h / 2
    ax.add_patch(Rectangle((cx - 0.010, cy - 0.014), 0.020, 0.027,
                           facecolor=BLUE_L, edgecolor=NAVY_L, lw=0.5, zorder=4))
    ax.add_patch(Ellipse((cx, cy + 0.013), 0.020, 0.007,
                         facecolor="#CFE2F3", edgecolor=NAVY_L, lw=0.5, zorder=5))
    ax.add_patch(Ellipse((cx, cy - 0.014), 0.020, 0.007,
                         facecolor=BLUE_L, edgecolor=NAVY_L, lw=0.5, zorder=5))
    say(x + 0.039, y + h - 0.012, title, size=7.1, color=NAVY, weight="bold")
    say(x + 0.039, y + h - 0.033, trait, size=6.4)
    say(x + 0.039, y + h - 0.052, "Summary statistics", size=5.9, color="#595959")


input_card(0.115, "Exposure GWAS", "Endometriosis")
input_card(0.585, "Outcome GWAS", "Venous thromboembolism")

ax.plot([0.265, 0.265], [0.862, 0.848], color=BLUE_L, lw=1.05, zorder=2)
ax.plot([0.735, 0.735], [0.862, 0.848], color=BLUE_L, lw=1.05, zorder=2)
ax.plot([0.265, 0.735], [0.848, 0.848], color=BLUE_L, lw=1.05, zorder=2)
arrow(0.500, 0.848, 0.500, 0.832)

# ---------------------------------------------------------------- stage bar
BAR_Y, BAR_H = 0.782, 0.050
BAR_X, BAR_W = 0.038, 0.924
round_box(BAR_X + 0.002, BAR_Y - 0.002, BAR_W, BAR_H, "#D9E2F3", "none",
          lw=0, z=2, rounding=0.008)
round_box(BAR_X, BAR_Y, BAR_W, BAR_H, NAVY, NAVY, lw=0.6, z=3, rounding=0.008)
seg_labels = ["Instrument selection", "Causal estimation",
              "Robustness & pleiotropy", "Integration & synthesis"]
seg_w = BAR_W / 4
for i, lab in enumerate(seg_labels):
    if i:
        ax.plot([BAR_X + i * seg_w, BAR_X + i * seg_w],
                [BAR_Y + 0.004, BAR_Y + BAR_H - 0.004],
                color="#5B7FA6", lw=0.7, zorder=4)
    say(BAR_X + (i + 0.5) * seg_w, BAR_Y + BAR_H / 2, lab, size=6.35,
        color="white", weight="bold", ha="center", va="center")

# ---------------------------------------------------------------- methods
column(COLS[0], [
    ("Instrument selection", [
        "Lead variant selection",
        "LD pruning",
        "Instrument strength screening",
    ], False),
    ("Coagulation-aware filtration", [
        "Coagulation-gene exclusion",
        "Sensitivity instrument set",
        "Locus-aware MR",
    ], False),
    ("External replication", [
        "Independent endometriosis GWAS",
        "Reverse-direction analysis",
        "Replication of MR estimates",
    ], False),
])

column(COLS[1], [
    ("Univariable MR", [
        "IVW and MR-Egger",
        "Weighted median and mode",
        "Heterogeneity assessment",
    ], False),
    ("Reverse MR", [
        "VTE instruments \u2192 endometriosis",
        "External endometriosis outcome",
        "Coagulation-locus sensitivity",
    ], False),
    ("Multivariable MR", [
        "Endometriosis, BMI, smoking",
        "Conditional instrument diagnostics",
        "Direct-effect estimation",
    ], False),
])

column(COLS[2], [
    ("Pleiotropy diagnostics", [
        "MR-Egger intercept",
        "MR-PRESSO and Cochran's Q",
        "Leave-one-out and Steiger",
    ], False),
    ("Colocalization", [
        "coloc.abf",
        "PLACO",
    ], False),
    ("Severity-stratified MR", [
        "PE-instrumented analysis",
        "DVT-instrumented analysis",
        "Instrument-subset comparison",
    ], False),
])

column(COLS[3], [
    ("Phenome-wide association", [
        "Coagulation-factor variant",
        "Thrombotic traits",
        "Gynaecological traits",
    ], False),
    ("Evidence synthesis", [
        "MR and colocalization",
        "MVMR and severity stratification",
        "Cross-trait consistency",
    ], True),
    ("Reporting", [
        "STROBE-MR",
        "Sensitivity analyses",
        "Reproducible workflow",
    ], False),
])

# ---------------------------------------------------------------- tree
ax.plot([0.500, 0.500], [BAR_Y, 0.754], color=BLUE_L, lw=1.05, zorder=2)
ax.plot([CENTERS[0], CENTERS[-1]], [0.754, 0.754], color=BLUE_L, lw=1.05, zorder=2)
for cx in CENTERS:
    arrow(cx, 0.754, cx, 0.736)

# ---------------------------------------------------------------- Layer 4
print("Layer 4 reconciliation: methods-only constraint")
bad = sorted({s for s in TEXTS if re.search(r"\d", s)})
for s in bad:
    print(f"  FAIL digit in displayed text: {s!r}")
if bad:
    raise ValueError("numeric content found; figure not written")

forbidden = [r"\bOR\b", r"\bCI\b", r"\bSNP\b", r"\bMb\b", r"\bPP\b",
             r"\bSE\b", r"\bI2\b", r"\bcases?\b", r"\bcontrols?\b",
             r"\bModel\b", r"\bFigure\b", r"\bStage\b"]
for pat in forbidden:
    hits = sorted({s for s in TEXTS if re.search(pat, s, flags=re.I)})
    for s in hits:
        print(f"  FAIL forbidden token {pat}: {s!r}")
    if hits:
        raise ValueError("non-method annotation found; figure not written")

print(f"  PASS no digits in {len(TEXTS)} displayed strings")
print("  PASS no sample sizes, estimates, thresholds or stage numbers")
print("  PASS methods-only labels")

for ext in ("png", "pdf", "tiff"):
    fig.savefig(os.path.join(FIG, f"Fig1_study_design.{ext}"), dpi=300,
                bbox_inches="tight", pad_inches=0.03, facecolor="white")
plt.close(fig)

# ---------------------------------------------------------------- trim margins
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
except Exception as exc:
    print(f"PDF crop skipped: {exc}")

print("saved: out/figures/Fig1_study_design.png|pdf|tiff")
