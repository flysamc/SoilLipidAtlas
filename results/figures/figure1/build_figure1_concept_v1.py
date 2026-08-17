"""Figure 1 concept rebuild — editable vector skeleton (v2, concept-only).

Rebuilds Figure 1 as a true SVG: every label is live text, every shape a
vector object, grouped per panel (a-d) so the file can be edited in
Illustrator/Inkscape/Affinity.

v2 design decision (2026-08-17): Figure 1 is a WORKFLOW figure — it carries
only design-scale facts (phyla, samples, batches, feature counts) and no
statistical results. Tier percentages, validation percentages, composite
weights, alignment tolerances, and the real composition forest plot were
removed; those results live in Figs 2 & 5 and Supp Fig 5. Panel d's output
plot is an explicitly schematic sketch. The organism count was dropped
entirely (the strict S1 producer records the species count as not reliably
derivable; submitted main text and S1 disagree at 110 vs 112).

Data-driven values that remain: 11,371 biomarkers / 16 phyla, 168/195 core
samples, 736 decoded fingerprint features, 273,248 consensus features.

Output: outputs/analysis/ncbi-phylum-2026-08-04-v1/figure1_redesign_2026-08-11_v1/
        Figure1_concept_skeleton.svg

Canvas: 180 x 120 mm (Nature Comms double column), viewBox 1080 x 720.
"""

import json
import math
import random
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REL = ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
OUT_DIR = REL / "figure1_redesign_2026-08-11_v1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- data inputs
bio = json.load(open(REL / "biomarker_discovery" / "summary.json"))
N_BIOMARKERS = bio["validation"]["atlas_rows"]          # 11371
N_BIO_PHYLA = bio["validation"]["observed_phyla"]       # 16

arch = json.load(open(REL / "climgrass" / "strict16_archlips_extended_2026-08-08_v1"
                      / "RUN_SUMMARY.json"))
N_SUBSTRATE = arch["simper_features"] + arch["archlips_added"]  # 736

sub = json.load(open(REL / "figure3" / "substrate" / "substrate_summary.json"))
N_POS_CORE = sub["population"]["collection_core_samples"]["POS"]   # 168
N_NEG_CORE = sub["population"]["collection_core_samples"]["NEG"]   # 195

# ------------------------------------------------------------------- palette
KING = {
    "Bacteria": "#7B52AB", "Archaea": "#1B9E8F", "Fungi": "#D9A420",
    "Plantae": "#3E9C35", "Animalia": "#D64541", "Protozoa": "#3B6FD4",
}
# locked ecological_group policy — display labels only. Data keys stay legacy.
DISPLAY = {"Plantae": "Viridiplantae", "Protozoa": "Protists"}
def disp(n): return DISPLAY.get(n, n)
TIER_COL = {"Gold": "#E8B931", "Silver": "#B8B8B8", "Bronze": "#C98A4B",
            "Unidentified": "#E3E3E3"}
TIER_ORDER = ["Gold", "Silver", "Bronze", "Unidentified"]
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

# panel divider lines
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


def icon_paths(name, col):
    """Line-art for one organism group, drawn in LOCAL coordinates centred on
    (0,0) inside a ~32 x 32 box. The caller wraps these in a <g> that carries the
    translate/scale, stroke colour and joins, so each icon is an independently
    editable object in Illustrator/Inkscape. Bodies get a light tint of their own
    colour so the silhouette still reads at 15 mm; details stay unfilled.
    """
    A = f' fill="{col}" fill-opacity="0.12"'      # tinted body
    p = []
    if name == "Bacteria":
        # three bacilli (rod-shaped cells) + a trailing flagellum
        for tx, ty, rot, L in ((-8, -7, -28, 19), (9, -4, 18, 18), (-2, 7, -7, 21)):
            p.append(f'<rect x="{-L/2:.1f}" y="-3.9" width="{L}" height="7.8" '
                     f'rx="3.9"{A} transform="translate({tx},{ty}) rotate({rot})"/>')
        p.append('<path d="M 9.5 9 q 4.5 -3.2 8 -0.4 q 3.4 2.8 6.2 0.2"/>')
    elif name == "Archaea":
        # coccoid packet (Methanosarcina-like), hexagonal close packing
        p.append(f'<circle cx="0" cy="0" r="4.6"{A}/>')
        for k in range(6):
            a = math.radians(k * 60 - 90)
            p.append(f'<circle cx="{9.9*math.cos(a):.1f}" cy="{9.9*math.sin(a):.1f}" '
                     f'r="4.6"{A}/>')
    elif name == "Fungi":
        # capped sporocarp + stipe + basal hyphae
        p.append(f'<path d="M -12 0 C -12 -10.5 -6.5 -14.5 0 -14.5 '
                 f'C 6.5 -14.5 12 -10.5 12 0 Z"{A}/>')
        p.append(f'<path d="M -3.2 0.5 C -3.8 5 -3.8 8.5 -5 12.5 L 5 12.5 '
                 f'C 3.8 8.5 3.8 5 3.2 0.5"{A}/>')
        p.append('<path d="M -5 12.5 c -3.2 0.6 -5 2 -6.4 3.8 '
                 'M 5 12.5 c 3.2 0.6 5 2 6.4 3.8"/>')
    elif name == "Plantae":
        # three tapered grass blades from a common base
        p.append(f'<path d="M 0 15 C -5 6 -9 -2 -12 -11 C -5.5 -2.5 -2 6 0 15 Z"{A}/>')
        p.append(f'<path d="M 0 15 C -2.5 4 -1.5 -6 1 -15.5 C 3.5 -6 3 4 0 15 Z"{A}/>')
        p.append(f'<path d="M 0 15 C 5 7 8.5 0.5 12.5 -8 C 6.5 0.5 2.5 7 0 15 Z"{A}/>')
    elif name == "Animalia":
        # oribatid soil mite, dorsal view: idiosoma + gnathosoma + EIGHT legs
        # (Acari are arachnids — the previous six-legged icon was incorrect)
        p.append(f'<ellipse cx="2.5" cy="0" rx="10.5" ry="7.2"{A}/>')
        p.append('<path d="M -7.6 -3.4 C -10 -3.6 -11.6 -2 -11.6 0 '
                 'C -11.6 2 -10 3.6 -7.6 3.4"/>')          # gnathosoma (mouthparts)
        # four leg pairs; front pairs sweep forward, rear pairs back, so the
        # outline reads as an arthropod rather than a radial starburst
        for x1, y1, x2, y2, x3, y3 in (
                (-4.5, -6, -9.5, -10.5, -15, -12.5),
                (-0.5, -7, -3.5, -12, -7.5, -15.5),
                (4, -6.9, 6.5, -12, 10.5, -14.5),
                (8, -5, 12.5, -8.5, 15.5, -10.5)):
            for s in (-1, 1):
                p.append(f'<path d="M {x1} {s*y1} L {x2} {s*y2} L {x3} {s*y3}"/>')
    else:
        # amoeboid protist: lobopodia + nucleus + contractile vacuole
        p.append(f'<path d="M -13.5 -3 C -15 -9 -9.5 -13.5 -4.5 -11 C -2.5 -16 5.5 -16 7 -10 '
                 f'C 13 -12.5 16.5 -5 12 -1 C 16 3.5 12.5 10 6.5 8 '
                 f'C 4.5 13.5 -3 14.5 -6 8.5 C -12 10.5 -16 3 -13.5 -3 Z"{A}/>')
        p.append(f'<circle cx="2.5" cy="-2" r="3.2" fill="{col}" fill-opacity="0.32"/>')
        p.append('<circle cx="-6.5" cy="3.5" r="1.9"/>')
    return p


for i, (name, col) in enumerate(groups):
    cx = gx0 + 12 + i * slot_w
    _glab = disp(name)
    text(cx + 18, gy + 20, _glab, size=10 if len(_glab) > 9 else 11,
         fill=col, weight="bold", anchor="middle")
    icy = gy + 42
    S.append(f'<g id="icon_{disp(name).lower()}" '
             f'transform="translate({cx+18},{icy}) scale(0.9)" '
             f'stroke="{col}" fill="none" stroke-width="1.8" '
             'stroke-linecap="round" stroke-linejoin="round">')
    for frag in icon_paths(name, col):
        S.append("  " + frag)
    S.append('</g>')

# stats row — design-scale facts only. Organism count intentionally omitted
# (strict S1: species count not reliably derivable); per-mode sample counts
# omitted too (168 POS vs 195 NEG reads as confusing in a workflow figure —
# they live in Table S1 and Methods).
sy = gy + 78
line(30, sy - 16, 526, sy - 16, stroke="#DDDDDD")
stats = ("19 collection phyla (16 in analysis)   |   "
         "6 measurement batches, two ionisation modes")
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

# alignment summary (concept only — tolerances live in Supplementary Method 1)
ay = py + 40
rect(30, ay, 496, 68, fill="white", stroke=GREY, sw=1, rx=6, dash="5 4")
text(160, ay + 28, "Cross-batch alignment", size=10.5, weight="bold",
     anchor="middle")
text(160, ay + 42, "(retention-time correction +", size=8.5, anchor="middle")
text(160, ay + 53, "m/z consensus, both modes)", size=8.5, anchor="middle")
arrow(288, ay + 34, 314, ay + 34)
rect(320, ay + 10, 180, 50, fill="#F1F1F1", stroke=INK, sw=1.2, rx=4)
text(410, ay + 25, "Consensus atlas", size=9.5, weight="bold", anchor="middle")
text(410, ay + 38, "273,248 features (positive)", size=8.5, anchor="middle")
text(410, ay + 50, "122,571 features (negative)", size=8.5, anchor="middle")
S.append('</g>')

# ============================================================ PANEL B =======
S.append('<g id="panel_b">')
panel_letter(562, 34, "b", "Two complementary analytical layers")

bx1, bx2, by, bw, bh = 562, 822, 50, 246, 280
# --- layer 1: biomarker atlas
rect(bx1, by, bw, bh, fill="white", stroke=LINE, sw=1.2, rx=8)
text(bx1 + bw / 2, by + 20, "1. Phylum-enriched biomarker atlas", size=11,
     weight="bold", fill=PURPLE, anchor="middle")
text(bx1 + bw / 2, by + 44, "Which lipids mark each phylum?", size=9.5,
     weight="bold", anchor="middle")
text(bx1 + bw / 2, by + 57, "(composite score + indicator statistics)", size=9,
     anchor="middle")
arrow(bx1 + bw / 2, by + 66, bx1 + bw / 2, by + 80)
text(bx1 + bw / 2, by + 96, f"{N_BIOMARKERS:,} biomarkers", size=12,
     weight="bold", fill=PURPLE, anchor="middle")
text(bx1 + bw / 2, by + 109, f"across {N_BIO_PHYLA} phyla (positive mode)",
     size=8.5, anchor="middle")
arrow(bx1 + bw / 2, by + 118, bx1 + bw / 2, by + 132)
text(bx1 + bw / 2, by + 148, "Annotation confidence tiers", size=9.5,
     weight="bold", fill=PURPLE, anchor="middle")
# schematic tier chips (equal widths — real proportions are in Fig. 2a)
chip_w = (bw - 44) / 4
for i, t in enumerate(TIER_ORDER):
    cxx = bx1 + 22 + i * chip_w
    rect(cxx, by + 156, chip_w - 3, 12, fill=TIER_COL[t], stroke=LINE, sw=0.5,
         rx=2)
lg_y = by + 184
for i, t in enumerate(TIER_ORDER):
    lx = bx1 + 22 + (i % 2) * ((bw - 44) / 2)
    ly = lg_y + (i // 2) * 14
    rect(lx, ly - 8, 9, 9, fill=TIER_COL[t], stroke=LINE, sw=0.5)
    lab = {"Gold": "Gold (molecular species)", "Silver": "Silver (partial)",
           "Bronze": "Bronze (lipid class)", "Unidentified": "Unidentified"}[t]
    text(lx + 13, ly, lab, size=7.5)
arrow(bx1 + bw / 2, lg_y + 22, bx1 + bw / 2, lg_y + 34)
text(bx1 + bw / 2, lg_y + 48, "External validation in", size=9.5,
     weight="bold", fill=GREY, anchor="middle")
text(bx1 + bw / 2, lg_y + 60, "public soil data (fastMASST)", size=9.5,
     weight="bold", fill=GREY, anchor="middle")

# --- layer 2: distributed fingerprints
rect(bx2, by, bw, bh, fill="white", stroke=LINE, sw=1.2, rx=8)
text(bx2 + bw / 2, by + 20, "2. Distributed fingerprint analysis", size=11,
     weight="bold", fill=BLUE, anchor="middle")
text(bx2 + bw / 2, by + 44, "How do whole lipidomes differ", size=9.5,
     weight="bold", anchor="middle")
text(bx2 + bw / 2, by + 57, "between phyla? (all quality features)", size=9.5,
     weight="bold", anchor="middle")
arrow(bx2 + bw / 2, by + 66, bx2 + bw / 2, by + 80)
text(bx2 + bw / 2, by + 96, "Lipidome similarity (Bray–Curtis)", size=9.5,
     weight="bold", anchor="middle")
text(bx2 + bw / 2, by + 109, "clustering & ordination recover", size=8.5,
     anchor="middle")
text(bx2 + bw / 2, by + 120, "taxonomic structure", size=8.5,
     anchor="middle", fill=GREY)
arrow(bx2 + bw / 2, by + 128, bx2 + bw / 2, by + 140)
text(bx2 + bw / 2, by + 154, "Fingerprint per phylum: its top", size=9.5,
     weight="bold", fill=BLUE, anchor="middle")
text(bx2 + bw / 2, by + 166, "distinguishing features (SIMPER)", size=9.5,
     weight="bold", fill=BLUE, anchor="middle")
# schematic mini-heatmap
hm_x, hm_y, cell = bx2 + 46, by + 176, 11
shades = ["#27346E", "#3D55A8", "#5E7BC4", "#8FA6DB", "#C3CfEC", "#E9EDF8"]
random.seed(11)
for r in range(4):
    lab = ["Phylum 1", "Phylum 2", "...", "Phylum n"][r]
    text(hm_x - 5, hm_y + r * cell + 9, lab, size=7.5, anchor="end")
    for c in range(9):
        col = shades[min(5, abs(c - r * 2) if r < 3 else random.randint(0, 5))]
        rect(hm_x + c * cell, hm_y + r * cell + 1, cell - 1.5, cell - 1.5,
             fill=col, rx=1.5)
# colorbar
cb_y = hm_y + 4 * cell + 6
for i, col in enumerate(reversed(shades)):
    rect(hm_x + 22 + i * 12, cb_y, 12, 6, fill=col)
text(hm_x + 18, cb_y + 6, "Low", size=7, anchor="end")
text(hm_x + 22 + 6 * 12 + 4, cb_y + 6, "High (contribution)", size=7)
arrow(bx2 + bw / 2, cb_y + 10, bx2 + bw / 2, cb_y + 18)
text(bx2 + bw / 2, cb_y + 30, "Confirmed by independent methods", size=9.5,
     weight="bold", fill=BLUE, anchor="middle")
text(bx2 + bw / 2, cb_y + 42, "(SCBD, CAP, L1; MS2LDA motifs)", size=8.5,
     anchor="middle")

# link arrow between layers
S.append(f'<line x1="{bx1+bw+2}" y1="{by+115}" x2="{bx2-2}" y2="{by+115}" '
         f'stroke="{GREY}" stroke-width="1" stroke-dasharray="4 3" '
         'marker-end="url(#arr)" marker-start="url(#arr)"/>')

# kingdom legend
ky = by + bh + 22
lx = 562
for name, col in KING.items():
    _llab = disp(name)
    rect(lx, ky - 8, 9, 9, fill=col, rx=2)
    text(lx + 13, ky, _llab, size=9)
    lx += 14 + 7.2 * len(_llab) + 12
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

# middle: MS/MS matching (thresholds live in Methods, not here)
mm_x = 250
text(mm_x + 60, cy0 + 10, "Soil lipidomes (positive mode)", size=10.5,
     weight="bold", anchor="middle")
text(mm_x + 60, cy0 + 44, "MS/MS spectral matching", size=9.5, weight="bold",
     anchor="middle")
text(mm_x + 60, cy0 + 56, "against the reference atlas", size=8.5,
     anchor="middle")
# mirrored spectra sketch
sp_x, sp_y = mm_x + 14, cy0 + 116
line(sp_x, sp_y, sp_x + 92, sp_y, stroke=INK, sw=0.8)
random.seed(7)
for i in range(9):
    fx_ = sp_x + 6 + i * 10
    up = 8 + random.random() * 26
    dn = 8 + random.random() * 26
    matched = i in (1, 3, 4, 6, 8)
    cu = PURPLE if matched else "#B9B9B9"
    line(fx_, sp_y, fx_, sp_y - up, stroke=cu, sw=1.4)
    line(fx_, sp_y, fx_, sp_y + dn, stroke=cu if matched else "#D3D3D3", sw=1.4)
    if matched:
        circle(fx_, sp_y - up - 3, 1.6, PURPLE)
text(sp_x + 46, sp_y + 44, "Soil spectrum (top) vs reference (bottom)",
     size=7.5, fill=GREY, anchor="middle")

# right: verified fingerprint features + mini heatmap
vm_x = 408
text(vm_x + 56, cy0 + 10, "Spectrally verified soil", size=10.5,
     weight="bold", anchor="middle")
text(vm_x + 56, cy0 + 22, f"substrate: {N_SUBSTRATE} diagnostic",
     size=10.5, weight="bold", anchor="middle")
text(vm_x + 56, cy0 + 34, f"fingerprint features ({N_BIO_PHYLA} phyla)",
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
text(274, dc_y + 15, "Each soil lipidome decomposed into phylum contributions",
     size=9.5, weight="bold", anchor="middle")
text(274, dc_y + 27, "(compared against every phylum's reference fingerprint)",
     size=8.5, fill=GREY, anchor="middle")
S.append('</g>')

# ============================================================ PANEL D =======
S.append('<g id="panel_d">')
panel_letter(562, 402, "d", "From lipid signal to community composition")

dy0 = 420
# correction chain box — plain-language phrasing of the three corrections
# (internally: Rule A response calibration, Rule C shared-biomarker split,
# ArchLips-restricted archaeal reference)
rect(562, dy0, 236, 108, fill="white", stroke=LINE, sw=1.2, rx=8)
text(680, dy0 + 17, "Making lipid signals comparable", size=10.5,
     weight="bold", anchor="middle")
corr = [("Calibrate", "signal weighted by lipid response"),
        ("Share", "shared biomarkers split across phyla"),
        ("Archaea", "only validated ether lipids counted")]
for i, (t1, t2) in enumerate(corr):
    yy = dy0 + 38 + i * 22
    rect(574, yy - 11, 52, 16, fill="#F1F1F1", stroke=LINE, sw=0.8, rx=3)
    text(600, yy, t1, size=8.5, weight="bold", anchor="middle")
    text(632, yy, t2, size=8.5)
    if i < 2:
        arrow(600, yy + 6, 600, yy + 10, sw=0.9)

# framework validation box — pointer only; the results live in Supp. Fig. 5
fv_y = dy0 + 122
rect(562, fv_y, 236, 64, fill="white", stroke=LINE, sw=1.2, rx=8)
text(680, fv_y + 18, "Framework validation", size=10.5, weight="bold",
     anchor="middle")
text(680, fv_y + 34, "Leave-one-out negative control on", size=8.5,
     anchor="middle")
text(680, fv_y + 46, "pure isolates (Supp. Fig. 5)", size=8.5,
     anchor="middle")

arrow(680, fv_y + 70, 680, fv_y + 84)
rect(562, fv_y + 90, 236, 40, fill=BOX, stroke=INK, sw=1, rx=6)
text(680, fv_y + 106, "Phylum-resolved community", size=9.5, weight="bold",
     anchor="middle")
text(680, fv_y + 120, "composition estimate (Fig. 5)", size=9.5, weight="bold",
     anchor="middle")

# schematic output sketch — ILLUSTRATIVE positions only (real values: Fig. 5)
fp_x, fp_y, fp_w, fp_h = 866, dy0 + 26, 176, 200
text(fp_x + fp_w / 2, dy0 + 4, "Lipid-derived composition", size=9.5,
     weight="bold", anchor="middle")
text(fp_x + fp_w / 2, dy0 + 16, "(schematic — see Fig. 5)", size=8.5,
     fill=GREY, anchor="middle", style="italic")
rows = ["Bacteria", "Archaea", "Fungi", "Plantae", "Animalia", "Protozoa"]
# fixed illustrative geometry: (range_lo, range_hi, dot) in axis units 0-50
SKETCH = {"Bacteria": (30, 44, 37), "Archaea": (3, 10, 6),
          "Fungi": (16, 30, 24), "Plantae": (14, 28, 20),
          "Animalia": (2, 8, 5), "Protozoa": (1, 6, 3)}
row_h = fp_h / len(rows)
def fx(v):
    return fp_x + fp_w * v / 50.0
for i, k in enumerate(rows):
    yy = fp_y + row_h * (i + 0.5)
    lab = disp(k)
    text(fp_x - 8, yy + 3, lab, size=8 if len(lab) > 9 else 9.5,
         fill=KING[k], weight="bold", anchor="end")
    lo, hi, dot = SKETCH[k]
    rect(fx(lo), yy - 6, fx(hi) - fx(lo), 12, fill="#DCDCDC")
    line(fx(max(lo - 4, 0.5)), yy, fx(min(hi + 4, 49)), yy, stroke=KING[k],
         sw=1.4)
    circle(fx(dot), yy, 3.4, KING[k])
# unlabelled axis line (schematic — no scale)
axy = fp_y + fp_h + 6
line(fp_x, axy, fp_x + fp_w, axy, stroke=INK, sw=1)
text(fp_x + fp_w / 2, axy + 14, "% of matched & corrected lipid signal",
     size=8.5, anchor="middle")
# legend
lgx, lgy = fp_x, axy + 32
rect(lgx, lgy - 8, 12, 10, fill="#DCDCDC")
text(lgx + 16, lgy, "Literature-expected range", size=7.5)
circle(lgx + 5, lgy + 12, 3.2, INK)
text(lgx + 16, lgy + 15, "Lipid-derived estimate (with CI)", size=7.5)
S.append('</g>')

S.append('</svg>')

out = OUT_DIR / "Figure1_concept_skeleton.svg"
out.write_text("\n".join(S), encoding="utf-8")
print("wrote", out)
print("design-scale numbers:", dict(
    biomarkers=N_BIOMARKERS, phyla=N_BIO_PHYLA, pos=N_POS_CORE,
    neg=N_NEG_CORE, substrate=N_SUBSTRATE))
