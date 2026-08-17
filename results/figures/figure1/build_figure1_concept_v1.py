"""Figure 1 concept rebuild — editable vector skeleton (v1).

Rebuilds the AI-generated Figure 1 mock as a true SVG: every label is live text,
every shape a vector object, grouped per panel (a-d) so the file can be edited
in Illustrator/Inkscape/Affinity. Data-driven values are read from the
ncbi-phylum-2026-08-04-v1 release outputs; unresolved items are rendered as
clearly-marked [TBD] slots.

Output: outputs/analysis/ncbi-phylum-2026-08-04-v1/figure1_redesign_2026-08-11_v1/
        Figure1_concept_skeleton.svg

Canvas: 180 x 120 mm (Nature Comms double column), viewBox 1080 x 720
(6 units per mm; 6 pt text = 2.117 mm = 12.7 units).
"""

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REL = ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
OUT_DIR = REL / "figure1_redesign_2026-08-11_v1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- data inputs
tier_counts = pd.read_csv(REL / "annotation" / "tier_counts.csv")
pos_tiers = (tier_counts[tier_counts["mode"] == "POS"]
             .groupby("annotation_tier")["n_features"].sum())
TIER_ORDER = ["Gold", "Silver", "Bronze", "Unidentified"]
tier_pct = {t: 100.0 * pos_tiers.get(t, 0) / pos_tiers.sum() for t in TIER_ORDER}

# fc-weighted + rules = the PRIMARY estimate in the final house-style render
# (r_render/fig6_climgrass_v2.R uses this composition; marker_panel not shown)
fcw = pd.read_csv(REL / "climgrass" / "figure5_redesign_2026-08-08_v2_archlips"
                  / "composition_fcweighted_kingdom_ci.csv")
fc = {r["kingdom"]: r for _, r in fcw.iterrows()}

simper = json.load(open(REL / "simper" / "summary.json"))
bio = json.load(open(REL / "biomarker_discovery" / "summary.json"))
N_BIOMARKERS = bio["validation"]["atlas_rows"]          # 11371
N_BIO_PHYLA = bio["validation"]["observed_phyla"]       # 16

arch = json.load(open(REL / "climgrass" / "strict16_archlips_extended_2026-08-08_v1"
                      / "RUN_SUMMARY.json"))
N_MATCHES = arch["spectral_qc"]["rows"]                 # 2313
N_SUBSTRATE = arch["simper_features"] + arch["archlips_added"]  # 736

sub = json.load(open(REL / "figure3" / "substrate" / "substrate_summary.json"))
N_POS_CORE = sub["population"]["collection_core_samples"]["POS"]   # 168
N_NEG_CORE = sub["population"]["collection_core_samples"]["NEG"]   # 195

N_FILTERED_POS = 45525  # figure3/substrate/pos.csv row count (verified 2026-08-11)

# literature-expected ranges (README figure5_redesign v2; development anchors)
EXPECTED = {"Bacteria": (35, 50), "Fungi": (15, 30), "Plantae": (15, 30),
            "Animalia": (1, 5), "Protozoa": (0, 1), "Archaea": (2, 8)}

# ------------------------------------------------------------------- palette
KING = {
    "Bacteria": "#7B52AB", "Archaea": "#1B9E8F", "Fungi": "#D9A420",
    "Plantae": "#3E9C35", "Animalia": "#D64541", "Protozoa": "#3B6FD4",
}
TIER_COL = {"Gold": "#E8B931", "Silver": "#B8B8B8", "Bronze": "#C98A4B",
            "Unidentified": "#E3E3E3"}
INK = "#1A1A1A"
GREY = "#6E6E6E"
LINE = "#BFBFBF"
BOX = "#FAFAFA"
PURPLE = "#7B52AB"
BLUE = "#2E5FA3"
FONT = "Helvetica, Arial, sans-serif"

S = []  # svg fragments


def esc(t):
    return (str(t).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def text(x, y, t, size=12, fill=INK, weight="normal", anchor="start",
         style="normal", spacing=None):
    extra = f' font-style="{style}"' if style != "normal" else ""
    if spacing:
        extra += f' letter-spacing="{spacing}"'
    S.append(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" '
             f'font-size="{size}" fill="{fill}" font-weight="{weight}" '
             f'text-anchor="{anchor}"{extra}>{esc(t)}</text>')


def rect(x, y, w, h, fill="none", stroke="none", sw=1, rx=0, dash=None, opac=1):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    o = f' opacity="{opac}"' if opac != 1 else ""
    S.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
             f'rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}{o}/>')


def line(x1, y1, x2, y2, stroke=LINE, sw=1, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    S.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
             f'stroke="{stroke}" stroke-width="{sw}"{d}/>')


def arrow(x1, y1, x2, y2, stroke=GREY, sw=1.2):
    S.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
             f'stroke="{stroke}" stroke-width="{sw}" marker-end="url(#arr)"/>')


def circle(cx, cy, r, fill, stroke="none", sw=1):
    S.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}" '
             f'stroke="{stroke}" stroke-width="{sw}"/>')


def panel_letter(x, y, letter, title):
    text(x, y, letter, size=19, weight="bold")
    text(x + 24, y, title, size=14.5, weight="bold")


# ================================================================== document
S.append(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1080 720" '
    'width="180mm" height="120mm" font-family="Helvetica, Arial, sans-serif">')
S.append('<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" '
         'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
         f'<path d="M 0 1 L 9 5 L 0 9 z" fill="{GREY}"/></marker></defs>')
S.append('<rect x="0" y="0" width="1080" height="720" fill="white"/>')

# panel divider lines (mock style)
line(545, 12, 545, 708, stroke="#E4E4E4", sw=1)
line(12, 372, 1068, 372, stroke="#E4E4E4", sw=1)

# ============================================================ PANEL A =======
S.append('<g id="panel_a">')
panel_letter(16, 34, "a", "Reference lipidome atlas from soil organisms")

# six organism groups with minimal icons
groups = list(KING.items())
gx0, gy = 30, 66
text(gx0, gy, "Six organism groups", size=11.5, weight="bold")
slot_w = 84
for i, (name, col) in enumerate(groups):
    cx = gx0 + 12 + i * slot_w
    text(cx + 18, gy + 20, name, size=11, fill=col, weight="bold", anchor="middle")
    icy = gy + 44
    S.append(f'<g id="icon_{name.lower()}" stroke="{col}" fill="none" stroke-width="1.6">')
    if name == "Bacteria":     # three rods
        for j, (dx, dy, rot) in enumerate([(-10, -4, -20), (4, -8, 15), (-2, 6, 5)]):
            S.append(f'<rect x="{cx+dx:.0f}" y="{icy+dy:.0f}" width="20" height="8" '
                     f'rx="4" transform="rotate({rot} {cx+dx+10} {icy+dy+4})"/>')
    elif name == "Archaea":    # cluster of cocci
        for dx, dy in [(-8, -6), (4, -8), (10, 2), (-2, 2), (-12, 5), (2, 10)]:
            S.append(f'<circle cx="{cx+dx+8:.0f}" cy="{icy+dy:.0f}" r="4.6"/>')
    elif name == "Fungi":      # mushroom
        S.append(f'<path d="M {cx-6} {icy} q 14 -18 28 0 z"/>')
        S.append(f'<path d="M {cx+5} {icy} l 0 12 m 4 -12 l 0 12"/>')
    elif name == "Plantae":    # grass blades
        S.append(f'<path d="M {cx+8} {icy+12} q -10 -8 -12 -22 M {cx+8} {icy+12} '
                 f'q 0 -14 2 -24 M {cx+8} {icy+12} q 10 -6 14 -20"/>')
    elif name == "Animalia":   # mite: body + legs
        S.append(f'<ellipse cx="{cx+8:.0f}" cy="{icy:.0f}" rx="10" ry="7"/>')
        for sgn in (-1, 1):
            for k in range(3):
                S.append(f'<path d="M {cx+8+sgn*9} {icy-4+k*4} l {sgn*8} {-2+k*2}"/>')
    else:                      # Protozoa: amoeba blob + flagellate
        S.append(f'<path d="M {cx-4} {icy} q -4 -10 6 -10 q 4 -6 10 -2 q 8 -2 8 6 '
                 f'q 6 6 -2 9 q -2 6 -9 4 q -8 4 -11 -2 q -6 -1 -2 -5 z"/>')
        S.append(f'<circle cx="{cx+7:.0f}" cy="{icy:.0f}" r="2" fill="{col}" stroke="none"/>')
    S.append('</g>')

# stats row (organism count = taxid-resolved analysis units, Figure 3b SSU
# curated freeze v3: 105 distinct NCBI taxids across the 16 analysis phyla)
units = pd.read_csv(REL / "figure3" / "evolutionary_tree_reference"
                    / "curated_freeze_v3"
                    / "ssu_curated_unit_representations_v3.csv")
N_ANALYSIS_ORG = units["analysis_unit_taxid"].nunique()   # 105
# + 4 collection-only organisms in the below-threshold phyla (Anabaena torulosa,
# lichen composite, Linnemannia gamsii, Rhogostoma minus) -> 109 collected
N_COLLECTED_ORG = N_ANALYSIS_ORG + 4
sy = gy + 78
line(30, sy - 16, 526, sy - 16, stroke="#DDDDDD")
stats = (f"{N_COLLECTED_ORG} organisms ({N_ANALYSIS_ORG} in statistics)   |   "
         "19 collection phyla (16 analysis)   |   "
         f"{N_POS_CORE} POS / {N_NEG_CORE} NEG samples   |   6 batches")
text(278, sy, stats, size=11, anchor="middle")
line(30, sy + 10, 526, sy + 10, stroke="#DDDDDD")

# pipeline row
py = sy + 52
steps = [("Biomass", "(sample material)"),
         ("Lipid extraction", "(modified Bligh–Dyer)"),
         ("C30 reversed-phase", "LC separation"),
         ("Orbitrap MS/MS", "(positive & negative)"),
         ("Per-batch processing", "(MZmine/GNPS2)")]
step_w, gap = 88, 14
px0 = 30
for i, (t1, t2) in enumerate(steps):
    x = px0 + i * (step_w + gap)
    rect(x, py - 26, step_w, 42, fill=BOX, stroke=LINE, sw=1, rx=5)
    text(x + step_w / 2, py - 9, t1, size=9.5, weight="bold", anchor="middle")
    text(x + step_w / 2, py + 4, t2, size=8.5, fill=GREY, anchor="middle")
    if i < 4:
        arrow(x + step_w + 1, py - 5, x + step_w + gap - 1, py - 5)

# cross-batch alignment box (dashed, mock style)
ay = py + 40
rect(30, ay, 496, 84, fill="white", stroke=GREY, sw=1, rx=6, dash="5 4")
text(278, ay + 17, "Cross-batch alignment (positive mode)", size=11,
     weight="bold", anchor="middle")
al_items = [(">15,000 anchor", "features (5 ppm / 2.0 min)"),
            ("LOESS retention", "time correction"),
            ("5 ppm m/z, 0.5 min RT", "consensus alignment")]
ax0 = 44
for i, (t1, t2) in enumerate(al_items):
    x = ax0 + i * 122
    text(x + 50, ay + 42, t1, size=9, anchor="middle")
    text(x + 50, ay + 54, t2, size=9, anchor="middle")
    if i < 2:
        arrow(x + 104, ay + 44, x + 120, ay + 44)
arrow(ax0 + 348, ay + 44, ax0 + 364, ay + 44)
rect(ax0 + 366, ay + 22, 116, 42, fill="#F1F1F1", stroke=INK, sw=1.2, rx=4)
text(ax0 + 424, ay + 36, "Consensus atlas", size=9.5, weight="bold", anchor="middle")
text(ax0 + 424, ay + 48, "273,248 features", size=8, anchor="middle")
text(ax0 + 424, ay + 58, "(1.5–25.0 min RT)", size=8, anchor="middle")
text(278, ay + 76, "Negative mode processed in parallel", size=9, fill=GREY,
     anchor="middle", style="italic")
S.append('</g>')

# ============================================================ PANEL B =======
S.append('<g id="panel_b">')
panel_letter(562, 34, "b", "Two complementary analytical layers")

bx1, bx2, by, bw, bh = 562, 822, 50, 246, 268
# --- layer 1: biomarker atlas
rect(bx1, by, bw, bh, fill="white", stroke=LINE, sw=1.2, rx=8)
text(bx1 + bw / 2, by + 20, "1. Phylum-enriched biomarker atlas", size=11,
     weight="bold", fill=PURPLE, anchor="middle")
text(bx1 + bw / 2, by + 40, "Cross-batch composite score", size=9.5,
     weight="bold", anchor="middle")
text(bx1 + bw / 2, by + 52,
     "(w₁=0.30, w₂=0.25, w₃=0.20, w₄=0.15, w₅=0.10)",
     size=8.5, anchor="middle")
text(bx1 + bw / 2, by + 64, "+ cross-batch IndVal (P < 0.05, FDR)", size=8.5,
     anchor="middle")
arrow(bx1 + bw / 2, by + 70, bx1 + bw / 2, by + 82)
text(bx1 + bw / 2, by + 96, f"{N_BIOMARKERS:,} biomarkers", size=11.5,
     weight="bold", fill=PURPLE, anchor="middle")
text(bx1 + bw / 2, by + 108, f"across {N_BIO_PHYLA} phyla (positive mode)",
     size=8.5, anchor="middle")
arrow(bx1 + bw / 2, by + 114, bx1 + bw / 2, by + 126)
text(bx1 + bw / 2, by + 140, "Annotation confidence (positive mode)", size=9.5,
     weight="bold", fill=PURPLE, anchor="middle")
# stacked annotation bar (data-driven)
bar_x, bar_y, bar_w, bar_h = bx1 + 22, by + 150, bw - 44, 14
cx = bar_x
for t in TIER_ORDER:
    w = bar_w * tier_pct[t] / 100.0
    rect(cx, bar_y, w, bar_h, fill=TIER_COL[t], stroke="white", sw=0.5)
    if tier_pct[t] >= 8:
        text(cx + w / 2, bar_y + 10.5, f"{tier_pct[t]:.1f}%", size=7.5,
             anchor="middle",
             fill="white" if t in ("Bronze",) else INK)
    cx += w
# tier legend (2 x 2)
lg_y = bar_y + 26
for i, t in enumerate(TIER_ORDER):
    lx = bar_x + (i % 2) * (bar_w / 2)
    ly = lg_y + (i // 2) * 14
    rect(lx, ly - 8, 9, 9, fill=TIER_COL[t], stroke=LINE, sw=0.5)
    lab = {"Gold": "Gold (molecular species)", "Silver": "Silver (partial)",
           "Bronze": "Bronze (lipid class)", "Unidentified": "Unidentified"}[t]
    text(lx + 13, ly, lab, size=7.5)
arrow(bx1 + bw / 2, lg_y + 24, bx1 + bw / 2, lg_y + 36)
text(bx1 + bw / 2, lg_y + 50, "External validation", size=9.5, weight="bold",
     fill=GREY, anchor="middle")
text(bx1 + bw / 2, lg_y + 62,
     "fastMASST | Pan-ReDU  [rerun pending for this atlas]", size=8,
     fill=GREY, anchor="middle", style="italic")

# --- layer 2: distributed fingerprints
rect(bx2, by, bw, bh, fill="white", stroke=LINE, sw=1.2, rx=8)
text(bx2 + bw / 2, by + 20, "2. Distributed fingerprint analysis", size=11,
     weight="bold", fill=BLUE, anchor="middle")
text(bx2 + bw / 2, by + 40, "Quality-filtered cross-batch feature space",
     size=9.5, weight="bold", anchor="middle")
text(bx2 + bw / 2, by + 52,
     f"({N_FILTERED_POS:,} positive-mode features", size=8.5, anchor="middle")
text(bx2 + bw / 2, by + 64,
     f"across {simper['n_phyla']} phyla for fingerprinting)", size=8.5,
     anchor="middle")
arrow(bx2 + bw / 2, by + 70, bx2 + bw / 2, by + 82)
text(bx2 + bw / 2, by + 96, "Bray–Curtis dissimilarity", size=9.5,
     weight="bold", anchor="middle")
text(bx2 + bw / 2, by + 108, "Ordination & UPGMA dendrograms", size=8.5,
     anchor="middle")
text(bx2 + bw / 2, by + 119, "(chemotaxonomic structure)", size=8.5,
     anchor="middle", fill=GREY)
arrow(bx2 + bw / 2, by + 125, bx2 + bw / 2, by + 137)
text(bx2 + bw / 2, by + 151, "SIMPER per phylum", size=9.5, weight="bold",
     fill=BLUE, anchor="middle")
text(bx2 + bw / 2, by + 163, "(top 3,000 features ranked)", size=8.5,
     anchor="middle")
# schematic mini-heatmap
hm_x, hm_y, cell = bx2 + 46, by + 172, 13
shades = ["#27346E", "#3D55A8", "#5E7BC4", "#8FA6DB", "#C3CfEC", "#E9EDF8"]
import random
random.seed(11)
for r in range(4):
    lab = ["Phylum 1", "Phylum 2", "...", "Phylum n"][r]
    text(hm_x - 5, hm_y + r * cell + 10, lab, size=7.5, anchor="end")
    ramp = sorted(random.sample(range(6), 6))
    for c in range(9):
        col = shades[min(5, abs(c - r * 2) if r < 3 else random.randint(0, 5))]
        rect(hm_x + c * cell, hm_y + r * cell + 1, cell - 1.5, cell - 1.5,
             fill=col, rx=1.5)
# colorbar
cb_y = hm_y + 4 * cell + 8
for i, col in enumerate(reversed(shades)):
    rect(hm_x + 22 + i * 12, cb_y, 12, 6, fill=col)
text(hm_x + 18, cb_y + 6, "Low", size=7, anchor="end")
text(hm_x + 22 + 6 * 12 + 4, cb_y + 6, "High (contribution)", size=7)
arrow(bx2 + bw / 2, cb_y + 12, bx2 + bw / 2, cb_y + 22)
text(bx2 + bw / 2, cb_y + 34, "Cross-method validation", size=9.5,
     weight="bold", fill=BLUE, anchor="middle")
text(bx2 + bw / 2, cb_y + 46, "SCBD | CAP | L1  &  MS2LDA motifs", size=8.5,
     anchor="middle")

# link arrow between layers
S.append(f'<line x1="{bx1+bw+2}" y1="{by+115}" x2="{bx2-2}" y2="{by+115}" '
         f'stroke="{GREY}" stroke-width="1" stroke-dasharray="4 3" '
         'marker-end="url(#arr)" marker-start="url(#arr)"/>')

# kingdom legend
ky = by + bh + 22
lx = 570
for name, col in KING.items():
    rect(lx, ky - 8, 9, 9, fill=col, rx=2)
    text(lx + 13, ky, name, size=9)
    lx += 14 + 7.2 * len(name) + 18
S.append('</g>')

# ============================================================ PANEL C =======
S.append('<g id="panel_c">')
panel_letter(16, 402, "c", "Soil community decoding (ClimGrass experiment)")

cy0 = 420
text(30, cy0 + 10, "ClimGrass 2 × 2 factorial design", size=10.5, weight="bold")
text(30, cy0 + 22, "n = 3 per treatment (12 soils)", size=9, fill=GREY)
# 2x2 grid
tg_x, tg_y, tg_w, tg_h = 62, cy0 + 62, 78, 54
cells = [("Ambient", "No drought", "#DFF0DB"), ("Ambient", "+ drought", "#FBEFC9"),
         ("Future climate", "No drought", "#D8E6F6"), ("Future climate", "+ drought", "#FADFD2")]
text(tg_x + tg_w, tg_y - 22, "Drought", size=9, weight="bold", anchor="middle")
text(tg_x + tg_w / 2, tg_y - 8, "No", size=8.5, anchor="middle")
text(tg_x + tg_w * 1.5, tg_y - 8, "Yes", size=8.5, anchor="middle")
S.append(f'<text x="{tg_x-32}" y="{tg_y+tg_h}" font-family="{FONT}" font-size="9" '
         f'font-weight="bold" fill="{INK}" text-anchor="middle" '
         f'transform="rotate(-90 {tg_x-32} {tg_y+tg_h})">Climate</text>')
text(tg_x - 16, tg_y + tg_h / 2 + 3, "Ambient", size=8, anchor="end")
text(tg_x - 16, tg_y + tg_h * 1.5 + 3, "Future", size=8, anchor="end")
for i, (r1, r2, col) in enumerate(cells):
    x = tg_x + (i % 2) * tg_w
    y = tg_y + (i // 2) * tg_h
    rect(x, y, tg_w - 3, tg_h - 3, fill=col, stroke=LINE, sw=0.8, rx=4)
    text(x + tg_w / 2 - 1, y + 20, r1, size=8, anchor="middle")
    text(x + tg_w / 2 - 1, y + 31, r2, size=8, anchor="middle")
    text(x + tg_w / 2 - 1, y + 43, "(n = 3)", size=7.5, fill=GREY, anchor="middle")

# middle: MS/MS matching
mm_x = 250
text(mm_x + 60, cy0 + 10, "Soil lipidomes (positive mode)", size=10.5,
     weight="bold", anchor="middle")
text(mm_x + 60, cy0 + 40, "Atlas matching by MS/MS", size=9.5, weight="bold",
     anchor="middle")
text(mm_x + 60, cy0 + 52, "(spectral cosine ≥ 0.70 and", size=8.5, anchor="middle")
text(mm_x + 60, cy0 + 64, "≥ 4 matched peaks)", size=8.5, anchor="middle")
# mirrored spectra sketch
sp_x, sp_y = mm_x + 14, cy0 + 116
line(sp_x, sp_y, sp_x + 92, sp_y, stroke=INK, sw=0.8)
random.seed(7)
for i in range(9):
    fx = sp_x + 6 + i * 10
    up = 8 + random.random() * 26
    dn = 8 + random.random() * 26
    matched = i in (1, 3, 4, 6, 8)
    cu = PURPLE if matched else "#B9B9B9"
    line(fx, sp_y, fx, sp_y - up, stroke=cu, sw=1.4)
    line(fx, sp_y, fx, sp_y + dn, stroke=cu if matched else "#D3D3D3", sw=1.4)
    if matched:
        circle(fx, sp_y - up - 3, 1.6, PURPLE)
text(sp_x + 46, sp_y + 44, "Soil spectrum (top) vs reference (bottom)",
     size=7.5, fill=GREY, anchor="middle")

# right: verified matches + mini heatmap
vm_x = 408
text(vm_x + 56, cy0 + 10, f"{N_MATCHES:,} verified spectral", size=10.5,
     weight="bold", anchor="middle")
text(vm_x + 56, cy0 + 22, "matches → " + f"{N_SUBSTRATE} diagnostic",
     size=10.5, weight="bold", anchor="middle")
text(vm_x + 56, cy0 + 34, f"fingerprint features ({simper['n_phyla']} phyla)",
     size=10.5, weight="bold", anchor="middle")
fm_x, fm_y, fc_cell = vm_x + 16, cy0 + 48, 12
random.seed(3)
for r in range(6):
    for c in range(7):
        col = shades[random.randint(0, 5)]
        rect(fm_x + c * fc_cell, fm_y + r * fc_cell, fc_cell - 1.5,
             fc_cell - 1.5, fill=col, rx=1.5)
S.append(f'<text x="{fm_x-6}" y="{fm_y+36}" font-family="{FONT}" font-size="7.5" '
         f'fill="{GREY}" text-anchor="middle" '
         f'transform="rotate(-90 {fm_x-6} {fm_y+36})">Soils</text>')
text(fm_x + 42, fm_y + 6 * fc_cell + 10, "Fingerprint features", size=7.5,
     fill=GREY, anchor="middle")

# bottom: decomposition box
dc_y = cy0 + 208
arrow(150, tg_y + 2 * tg_h + 8, 150, dc_y - 4)
arrow(vm_x + 56, fm_y + 6 * fc_cell + 16, vm_x + 56, dc_y - 4)
rect(96, dc_y, 356, 34, fill=BOX, stroke=INK, sw=1, rx=5)
text(274, dc_y + 15, "Fold-change-weighted Bray–Curtis decomposition",
     size=9.5, weight="bold", anchor="middle")
text(274, dc_y + 27, "(against phylum reference fingerprints)", size=8.5,
     fill=GREY, anchor="middle")
S.append('</g>')

# ============================================================ PANEL D =======
S.append('<g id="panel_d">')
panel_letter(562, 402, "d", "Lipid framework: correction and lipid-derived output")

dy0 = 420
# correction chain box
rect(562, dy0, 236, 108, fill="white", stroke=LINE, sw=1.2, rx=8)
text(680, dy0 + 17, "Quantification correction", size=10.5, weight="bold",
     anchor="middle")
corr = [("Rule A", "out-of-range RIE → uncalibrated (1.0)"),
        ("Rule C", "enriched weight split across k phyla"),
        ("ArchLips", "restricted archaeal reference")]
for i, (t1, t2) in enumerate(corr):
    yy = dy0 + 38 + i * 22
    rect(574, yy - 11, 52, 16, fill="#F1F1F1", stroke=LINE, sw=0.8, rx=3)
    text(600, yy, t1, size=8.5, weight="bold", anchor="middle")
    text(632, yy, t2, size=8.5)
    if i < 2:
        arrow(600, yy + 6, 600, yy + 10, sw=0.9)

# framework validation box (placeholder — producer not rerun)
fv_y = dy0 + 122
rect(562, fv_y, 236, 92, fill="#FCFCFC", stroke=LINE, sw=1.2, rx=8, dash="5 4")
text(680, fv_y + 17, "Framework validation", size=10.5, weight="bold",
     anchor="middle")
text(680, fv_y + 40, "[TBD — negative-control rerun pending", size=8.5,
     fill=GREY, anchor="middle", style="italic")
text(680, fv_y + 52, "on the current substrate; old n=163 /", size=8.5,
     fill=GREY, anchor="middle", style="italic")
text(680, fv_y + 64, "80% / <2% claims are pre-correction]", size=8.5,
     fill=GREY, anchor="middle", style="italic")

# forest plot (data-driven)
fp_x, fp_y, fp_w, fp_h = 866, dy0 + 26, 176, 200
axis_max = 50.0
def fx(v):
    return fp_x + fp_w * min(v, axis_max) / axis_max

text(fp_x + fp_w / 2, dy0 + 4, "Lipid-derived proportion of", size=9.5,
     weight="bold", anchor="middle")
text(fp_x + fp_w / 2, dy0 + 16, "matched & corrected signal", size=9.5,
     weight="bold", anchor="middle")
rows = ["Bacteria", "Archaea", "Fungi", "Plantae", "Animalia", "Protozoa"]
STATUS_COL = {"within": "#3F8E3F", "above": "#C46210", "below": "#666666"}
row_h = fp_h / len(rows)
for i, k in enumerate(rows):
    yy = fp_y + row_h * (i + 0.5)
    lab = "Archaea†" if k == "Archaea" else k
    text(fp_x - 8, yy + 3, lab, size=9.5, fill=KING[k], weight="bold", anchor="end")
    lo, hi = EXPECTED[k]
    rect(fx(lo), yy - 6, fx(hi) - fx(lo), 12, fill="#DCDCDC")
    m = fc[k]
    mv, ml, mh = 100 * m["mean"], 100 * m["ci_lo"], 100 * m["ci_hi"]
    line(fx(ml), yy, fx(mh), yy, stroke=KING[k], sw=1.4)
    circle(fx(mv), yy, 3.4, KING[k])
    status = "within" if lo <= mv <= hi else ("above" if mv > hi else "below")
    sc = STATUS_COL[status]
    sx, sy_ = fx(min(mh, axis_max)) + 8, yy
    if status == "within":
        circle(sx, sy_, 2.6, sc)
    elif status == "above":
        S.append(f'<path d="M {sx:.1f} {sy_-3.4:.1f} l 3.4 6 l -6.8 0 z" '
                 f'fill="{sc}"/>')
    else:
        S.append(f'<path d="M {sx:.1f} {sy_+3.4:.1f} l 3.4 -6 l -6.8 0 z" '
                 f'fill="{sc}"/>')
# axis
axy = fp_y + fp_h + 6
line(fp_x, axy, fp_x + fp_w, axy, stroke=INK, sw=1)
for v in (0, 10, 20, 30, 40, 50):
    line(fx(v), axy, fx(v), axy + 3, stroke=INK, sw=1)
    text(fx(v), axy + 13, str(v), size=8, anchor="middle")
text(fp_x + fp_w / 2, axy + 26, "% of matched & corrected lipid signal",
     size=8.5, anchor="middle")
# legend (left column, under the validation box — clear of the plot)
lgx, lgy = 574, fv_y + 112
rect(lgx, lgy - 8, 12, 10, fill="#DCDCDC")
text(lgx + 16, lgy, "Literature-expected range", size=7.5)
circle(lgx + 5, lgy + 12, 3.2, INK)
text(lgx + 16, lgy + 15, "fc-weighted estimate (mean, 95% CI)", size=7.5)
circle(lgx + 3, lgy + 26, 2.4, STATUS_COL["within"])
text(lgx + 9, lgy + 29, "within", size=7.5, fill=GREY)
S.append(f'<path d="M {lgx+46:.0f} {lgy+23:.0f} l 3 5.4 l -6 0 z" '
         f'fill="{STATUS_COL["above"]}"/>')
text(lgx + 53, lgy + 29, "above", size=7.5, fill=GREY)
S.append(f'<path d="M {lgx+90:.0f} {lgy+28.4:.0f} l 3 -5.4 l -6 0 z" '
         f'fill="{STATUS_COL["below"]}"/>')
text(lgx + 97, lgy + 29, "below expected", size=7.5, fill=GREY)
# note
text(816, 700, "Note: values are lipid-derived proportions of matched & "
     "corrected signal, not absolute biomass fractions.", size=8, fill=GREY,
     anchor="middle", style="italic")
text(816, 712, "Archaea scale uncertain without ether-lipid RIE standards "
     "(†).  Animalia reflects the known holobiont-reference bias.",
     size=8, fill=GREY, anchor="middle", style="italic")
S.append('</g>')

S.append('</svg>')

out = OUT_DIR / "Figure1_concept_skeleton.svg"
out.write_text("\n".join(S), encoding="utf-8")
print("wrote", out)
print("tier pct:", {k: round(v, 1) for k, v in tier_pct.items()})
print("fcw mean [ci]:", {k: (round(100 * v["mean"], 1), round(100 * v["ci_lo"], 1),
                             round(100 * v["ci_hi"], 1)) for k, v in fc.items()})