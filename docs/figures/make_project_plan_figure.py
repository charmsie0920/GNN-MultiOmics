"""Render the full project-plan architecture as a standalone SVG.

Laid out like MoGraphDRP Fig. 1 so the two can be read side by side:
(A) feature encoders, (B) our biological-graph stage, (C) interaction and
prediction, (D) ensemble refinement. Every component carries a tag saying
whether it is ours (changed from or added to MoGraphDRP) or kept from the
benchmark. Numbers come from docs/architecture.md and src/models/.

    python docs/figures/make_project_plan_figure.py
"""

import random
from pathlib import Path

OUT = Path(__file__).with_name("project_plan_architecture.svg")
W, H = 1900, 965

INK, MUTED, PANEL, PANEL_STROKE = "#1E2530", "#5B6472", "#F4F6F9", "#D5DAE1"
CELL, CELL_L = "#2E6DB4", "#DCE8F6"
DRUG, DRUG_L = "#D9782A", "#FBE7D5"
PROT, PROT_L = "#3E9A6A", "#DDF0E5"
PPI = "#A9C9B5"
GE_C, MC_C, PR_C = "#2E6DB4", "#7B52AB", "#2A8C8C"
F_C = "#8E6BB8"
OURS = "#C0392B"
KEPT = "#6B7480"

out: list[str] = []


def add(s: str) -> None:
    out.append(s)


def rect(x, y, w, h, fill, stroke="none", rx=0, sw=1.0, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="{sw}"{d}/>')


def circle(cx, cy, r, fill, stroke="none", sw=1.0):
    add(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def line(x1, y1, x2, y2, stroke=INK, sw=1.3, arrow=None, dash=None, both=False):
    m = f' marker-end="url(#ar-{arrow})"' if arrow else ""
    if both and arrow:
        m += f' marker-start="url(#as-{arrow})"'
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{sw}"{m}{d}/>')


def poly(points, stroke=INK, sw=1.3, arrow="ink", dash=None):
    p = " ".join(f"{x},{y}" for x, y in points)
    m = f' marker-end="url(#ar-{arrow})"' if arrow else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<polyline points="{p}" fill="none" stroke="{stroke}" stroke-width="{sw}"{m}{d}/>')


def text(x, y, s, size=12, anchor="middle", weight=400, fill=INK, italic=False):
    it = ' font-style="italic"' if italic else ""
    add(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" '
        f'font-weight="{weight}" fill="{fill}"{it}>{s}</text>')


def lines(x, y, rows, size=11, anchor="middle", step=15, fill=INK, weight=400):
    for i, r in enumerate(rows):
        text(x, y + i * step, r, size, anchor, weight, fill)


def pill(x, y, kind):
    """Tag at (x, y) = left edge, vertical centre. Returns its width."""
    label = "OURS" if kind == "ours" else "MoGraphDRP"
    w = 8 + len(label) * 6.6
    if kind == "ours":
        rect(x, y - 8, w, 16, OURS, rx=8)
        text(x + w / 2, y + 3.5, label, 9.5, weight=700, fill="#FFFFFF")
    else:
        rect(x, y - 8, w, 16, "#FFFFFF", KEPT, rx=8, sw=1.1)
        text(x + w / 2, y + 3.5, label, 9.5, weight=700, fill=KEPT)
    return w


def panel(x, y, w, h, letter, title, subtitle, kind):
    rect(x, y, w, h, PANEL, PANEL_STROKE, rx=10)
    text(x + 16, y + 30, letter, 18, "start", 700)
    text(x + 38, y + 30, title, 15, "start", 700)
    pill(x + 44 + len(title.replace('&amp;', '&')) * 8.4, y + 25, kind)
    text(x + 16, y + 49, subtitle, 11, "start", 400, MUTED)


def vbar(x, y, w, h, color, n=6):
    """A vector drawn as a column of n cells."""
    step = h / n
    for i in range(n):
        rect(x, y + i * step, w, step, color if i % 2 == 0 else color, "#FFFFFF", sw=0.8)
    rect(x, y, w, h, "none", color, sw=1)


def hbar(x, y, w, h, color, n=8):
    step = w / n
    for i in range(n):
        rect(x + i * step, y, step, h, color, "#FFFFFF", sw=0.8)
    rect(x, y, w, h, "none", color, sw=1)


def heatmap(x, y, w, h, base, seed):
    rnd = random.Random(seed)
    cols, rows = 6, 4
    cw, ch = w / cols, h / rows
    for r in range(rows):
        for c in range(cols):
            op = 0.25 + 0.75 * rnd.random()
            add(f'<rect x="{x + c * cw:.1f}" y="{y + r * ch:.1f}" width="{cw:.1f}" height="{ch:.1f}" '
                f'fill="{base}" fill-opacity="{op:.2f}" stroke="#FFFFFF" stroke-width="0.6"/>')
    rect(x, y, w, h, "none", base, sw=1)


def box(x, y, w, h, rows, stroke, fill="#FFFFFF", size=11, bold_first=True, step=14):
    rect(x, y, w, h, fill, stroke, rx=6, sw=1.4)
    top = y + h / 2 - (len(rows) - 1) * step / 2 + 4
    for i, r in enumerate(rows):
        text(x + w / 2, top + i * step, r, size, weight=700 if (bold_first and i == 0) else 400)


def sup(base, s):
    return f'{base}<tspan dy="-5" font-size="0.72em">{s}</tspan><tspan dy="5"></tspan>'


def sub(base, s):
    return f'{base}<tspan dy="3" font-size="0.75em">{s}</tspan><tspan dy="-3"></tspan>'


# ---------------------------------------------------------------- header
add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
    f'font-family="Helvetica, Arial, sans-serif">')
add("<defs>")
for name, col in [("ink", INK), ("cell", CELL), ("drug", DRUG), ("prot", PROT), ("ours", OURS)]:
    add(f'<marker id="ar-{name}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="{col}"/></marker>')
    add(f'<marker id="as-{name}" viewBox="0 0 10 10" refX="1" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto"><path d="M10,0 L0,5 L10,10 z" fill="{col}"/></marker>')
add("</defs>")
rect(0, 0, W, H, "#FFFFFF")

# ================================================================ PANEL A
panel(20, 20, 760, 880, "A", "Feature encoders", "cell line and drug are encoded separately, as in MoGraphDRP", "ours")

# ---- cell line
text(40, 100, "Cell-line multi-omics · 532 cell lines", 13, "start", 700)
pill(300, 96, "ours")
text(40, 117, "MoGraphDRP: expression, mutation, methylation, pathway · COSMIC genes · per-omics MLP → concat",
     10, "start", fill=MUTED, italic=True)

omics = [
    ("Transcriptomics", "RNA-seq TPM", "log1p · scale", GE_C, "GE"),
    ("Genomics", "Mutation + CNV", "scale", MC_C, "MC"),
    ("Proteomics", "mass spec", "impute · scale", PR_C, "PR"),
]
fx, fy, fw, fh = 478, 140, 222, 236      # fusion block
node = {"GE": (512, 205), "MC": (512, 305), "PR": (578, 255)}
for i, (name, src, prep, col, key) in enumerate(omics):
    y0 = 145 + i * 80
    heatmap(40, y0, 56, 44, col, i + 1)
    text(106, y0 + 18, name, 12.5, "start", 700)
    text(106, y0 + 34, src, 10.5, "start", fill=MUTED)
    line(222, y0 + 22, 246, y0 + 22, arrow="ink")
    rect(248, y0, 118, 44, "#FFFFFF", col, rx=6, sw=1.4)
    text(307, y0 + 18, prep, 10.5)
    text(307, y0 + 34, "PCA → 128", 11.5, weight=700)
    line(366, y0 + 22, 386, y0 + 22, arrow="ink")
    vbar(389, y0 + 4, 14, 36, col, 4)
    text(396, y0 + 54, "128", 9.5, fill=MUTED)
    line(403, y0 + 22, fx - 2, y0 + 22 if i != 1 else y0 + 22, arrow="ink")

# fusion block
rect(fx, fy - 20, fw, fh + 20, CELL_L, CELL, rx=8, sw=1.4)
text(fx + 12, fy, "Cross-attention fusion", 12.5, "start", 700)
pill(fx + 158, fy - 4, "ours")
for a, b in [("GE", "MC"), ("GE", "PR"), ("MC", "PR")]:
    (x1, y1), (x2, y2) = node[a], node[b]
    dx, dy = x2 - x1, y2 - y1
    L = (dx * dx + dy * dy) ** 0.5
    ux, uy = dx / L, dy / L
    line(x1 + ux * 19, y1 + uy * 19, x2 - ux * 19, y2 - uy * 19, CELL, 1.4, arrow="cell", both=True)
for key, (cx, cy) in node.items():
    col = {"GE": GE_C, "MC": MC_C, "PR": PR_C}[key]
    circle(cx, cy, 16, col)
    text(cx, cy + 4, key, 10, weight=700, fill="#FFFFFF")
lines(606, 190, ["6 directed pairs", "q ← omics i", "k, v ← omics j", "MHA · 4 heads", "+ FFN · LayerNorm"],
      10.5, "start", 16)
text(fx + fw / 2, fy + 196, "concat 6 × 128 → Linear → 256", 10.5)
text(fx + fw / 2, fy + 211, "GELU · LayerNorm", 10.5, fill=MUTED)
line(fx + fw, 255, 716, 255, arrow="ink")
vbar(719, 205, 16, 100, CELL, 8)
text(727, 322, sub("x", "c") + " ∈ ℝ²⁵⁶", 10.5)
text(727, 336, "cell vector", 9.5, fill=MUTED)

# divider
line(40, 432, 760, 432, PANEL_STROKE, 1)

# ---- drug
text(40, 460, "Drug structure · 498 drugs", 13, "start", 700)
text(40, 477, "MoGraphDRP: Morgan + ESPF + PubChem fingerprints, attention-based drug fusion",
     10, "start", fill=MUTED, italic=True)

# SMILES hexagon
hx, hy, hr = 78, 660, 24
pts = " ".join(f"{hx + hr * c:.1f},{hy + hr * s:.1f}" for c, s in
               [(0, -1), (0.866, -0.5), (0.866, 0.5), (0, 1), (-0.866, 0.5), (-0.866, -0.5)])
add(f'<polygon points="{pts}" fill="{DRUG_L}" stroke="{DRUG}" stroke-width="1.6"/>')
circle(hx, hy, 10, "none", DRUG, 1.2)
text(hx, hy + 42, "SMILES", 11, weight=700)
text(hx, hy + 56, "PubChem", 10, fill=MUTED)
line(106, hy, 128, hy, arrow="ink")
box(130, hy - 18, 64, 36, ["RDKit"], DRUG)
line(194, hy, 214, hy)
line(214, 560, 214, 780)
line(214, 560, 238, 560, arrow="ink")
line(214, 780, 238, 780, arrow="ink")

# graph branch (kept)
text(240, 512, "Graph branch", 11.5, "start", 700)
pill(326, 508, "kept")
atoms = [(252, 548), (272, 536), (292, 548), (292, 572), (272, 584), (252, 572), (312, 536), (232, 536)]
for a, b in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (2, 6), (0, 7)]:
    line(*atoms[a], *atoms[b], MUTED, 1.2)
for i, (ax, ay) in enumerate(atoms):
    circle(ax, ay, 5, "#FFFFFF" if i % 3 else DRUG_L, DRUG, 1.2)
text(272, 608, "atom graph · 39-d atoms", 9.5, fill=MUTED)
line(322, 560, 346, 560, arrow="ink")
for k in range(3):
    rect(350 + k * 7, 530 + k * 5, 50, 52, CELL_L, CELL, rx=2, sw=1)
text(385, 608, "GCN × 3", 11, weight=700)
text(385, 622, "39 → 128", 9.5, fill=MUTED)
line(416, 560, 436, 560, arrow="ink")
box(438, 542, 70, 36, ["mean pool"], DRUG, size=10.5)
line(508, 560, 528, 560, arrow="ink")
vbar(531, 538, 14, 44, DRUG, 4)
text(538, 596, "128", 9.5, fill=MUTED)

# fingerprint branch (ours)
text(240, 728, "Fingerprint branch", 11.5, "start", 700)
pill(358, 724, "ours")
bits = [1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0]
for i, b in enumerate(bits):
    rect(240 + i * 6, 768, 6, 24, DRUG if b else "#FFFFFF", DRUG, sw=0.6)
text(282, 808, "Morgan · r = 2", 9.5, fill=MUTED)
text(282, 821, "2048 bits", 9.5, fill=MUTED)
line(326, 780, 346, 780, arrow="ink")
box(348, 758, 160, 44, ["FP encoder", "2048 → 128 → 128"], DRUG, size=10.5)
line(508, 780, 528, 780, arrow="ink")
vbar(531, 758, 14, 44, DRUG, 4)
text(538, 816, "128", 9.5, fill=MUTED)

# drug fusion
circle(610, 670, 15, "#FFFFFF", DRUG, 1.4)
text(610, 675, "‖", 13, weight=700)
poly([(545, 560), (610, 560), (610, 653)], arrow="ink")
poly([(545, 780), (610, 780), (610, 687)], arrow="ink")
text(632, 648, "concat", 10, "start", fill=MUTED)
line(625, 670, 716, 670, arrow="ink")
vbar(719, 620, 16, 100, DRUG, 8)
text(727, 737, sub("x", "d") + " ∈ ℝ²⁵⁶", 10.5)
text(727, 751, "drug vector", 9.5, fill=MUTED)

# A -> B arrows (node features into the graph)
line(735, 255, 846, 255, CELL, 1.6, arrow="cell")
line(735, 670, 846, 670, DRUG, 1.6, arrow="drug")
text(788, 247, "node", 9.5, fill=CELL)
text(788, 662, "node", 9.5, fill=DRUG)

# ================================================================ PANEL B
panel(800, 20, 520, 880, "B", "Biological graph learning", "heterogeneous GNN over the STRING PPI network · not in MoGraphDRP", "ours")

text(862, 110, "cell_line", 12, weight=700, fill=CELL)
text(990, 110, "protein", 12, weight=700, fill=PROT)
text(862, 540, "drug", 12, weight=700, fill=DRUG)

cells = [(862, 150), (862, 205), (862, 260), (862, 315), (862, 370)]
drugs = [(862, 580), (862, 635), (862, 690), (862, 745)]
prots = [(975, 150), (1040, 190), (990, 240), (1050, 290), (970, 340), (1035, 390),
         (985, 440), (1050, 490), (970, 545), (1035, 595), (985, 650), (1045, 705), (975, 760)]
ppi = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (2, 4), (3, 5), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7),
       (6, 8), (7, 9), (8, 9), (8, 10), (9, 10), (9, 11), (10, 11), (10, 12), (11, 12), (3, 7)]
hi_ppi = {(4, 6), (6, 8)}
for a, b in ppi:
    red = (a, b) in hi_ppi
    line(*prots[a], *prots[b], OURS if red else PPI, 2.4 if red else 1.4)
mut = [(0, 0), (1, 2), (2, 2), (3, 3), (4, 4)]
for c, p in mut:
    red = (c, p) == (4, 4)
    line(cells[c][0] + 12, cells[c][1], prots[p][0] - 8, prots[p][1], OURS if red else CELL, 2.4 if red else 1.1)
tgt = [(0, 8), (1, 8), (2, 10), (3, 12)]
for d, p in tgt:
    red = (d, p) == (0, 8)
    line(drugs[d][0] + 12, drugs[d][1], prots[p][0] - 8, prots[p][1], OURS if red else DRUG, 2.4 if red else 1.1)
for cx, cy in cells:
    circle(cx, cy, 12, CELL_L, CELL, 1.5)
for cx, cy in drugs:
    rect(cx - 11, cy - 11, 22, 22, DRUG_L, DRUG, rx=4, sw=1.5)
for cx, cy in prots:
    circle(cx, cy, 8, PROT_L, PROT, 1.4)

# relation legend under the graph
ly = 800
rows = [("has_mutation", "4,341 driver-mutation edges", CELL),
        ("targets", "683 GDSC drug–target edges", DRUG),
        ("interacts_with", "473,860 STRING edges, score ≥ 700", PROT)]
for i, (rel, desc, col) in enumerate(rows):
    line(830, ly + i * 20 - 4, 852, ly + i * 20 - 4, col, 2)
    text(860, ly + i * 20, rel, 10.5, "start", 700, col)
    text(952, ly + i * 20, desc, 10.5, "start", fill=MUTED)
text(830, ly + 70, "red path: a mutated protein and a drug target meet in 2 hops", 10, "start", fill=OURS)

# HeteroConv block
bx, by, bw, bh = 1090, 96, 214, 520
rect(bx, by, bw, bh, "#FFFFFF", INK, rx=8, sw=1.4)
text(bx + bw / 2, by + 22, "HeteroConv × 2 layers", 12.5, weight=700)
text(bx + bw / 2, by + 38, "one GraphSAGE conv per relation", 10, fill=MUTED)
text(bx + 14, by + 64, "node encoders → D = 128", 10.5, "start", 700)
lines(bx + 14, by + 80, ["cell_line: Linear 256 → 128",
                         "drug: Linear 256 → 128",
                         "protein: Embedding 16,214 × 128"], 10, "start", 14, MUTED)
rels = [("interacts_with", "protein → protein", PROT_L, PROT),
        ("has_mutation", "cell_line → protein", CELL_L, CELL),
        ("targets", "drug → protein", DRUG_L, DRUG),
        ("rev_has_mutation", "protein → cell_line", PROT_L, PROT),
        ("rev_targets", "protein → drug", PROT_L, PROT)]
for i, (r, d, f, s) in enumerate(rels):
    ry = by + 136 + i * 50
    rect(bx + 14, ry, bw - 28, 40, f, s, rx=5, sw=1.2)
    text(bx + bw / 2, ry + 17, r, 11, weight=700, italic=True)
    text(bx + bw / 2, ry + 31, d, 9.5, fill=MUTED)
text(bx + bw / 2, by + 402, "Σ over relations → ReLU", 11, weight=700)
lines(bx + bw / 2, by + 422, ["full-batch message passing;", "only label pairs are mini-batched"], 10, step=14, fill=MUTED)
rect(bx + 14, by + 460, bw - 28, 44, "#FFF4F2", OURS, rx=5, sw=1, dash="4 3")
lines(bx + bw / 2, by + 478, ["GAT in the proposal → GraphSAGE", "(GAT diverged in 3 runs)"], 10, step=14, fill=OURS)
line(1062, 400, bx - 2, 400, arrow="ink")

# outputs
line(bx + bw / 2, by + bh, bx + bw / 2, 640, arrow="ink")
hbar(1120, 648, 150, 16, CELL, 8)
text(1197, 682, sub("h", "c") + " ∈ ℝ¹²⁸  refined cell", 10.5)
hbar(1120, 700, 150, 16, DRUG, 8)
text(1197, 734, sub("h", "d") + " ∈ ℝ¹²⁸  refined drug", 10.5)
# B -> C
line(1270, 656, 1292, 656)
line(1270, 708, 1292, 708)
line(1292, 656, 1292, 708)
poly([(1292, 682), (1333, 682), (1333, 205), (1374, 205)], arrow="ink", sw=1.5)

# ================================================================ PANEL C
panel(1345, 20, 535, 560, "C", "Interaction &amp; prediction", "multi-head bilinear attention, then an MLP regressor", "kept")

line(1374, 150, 1374, 265)
vbar(1394, 110, 16, 80, CELL, 6)
vbar(1394, 225, 16, 80, DRUG, 6)
line(1374, 150, 1392, 150, arrow="ink")
line(1374, 265, 1392, 265, arrow="ink")
text(1402, 104, sub("h", "c"), 10.5)
text(1402, 320, sub("h", "d"), 10.5)
line(1410, 150, 1436, 150, arrow="ink")
line(1410, 265, 1436, 265, arrow="ink")
box(1438, 124, 112, 52, ["FCNet", "128 → 4 × 128", "ReLU · dropout"], INK, size=10, step=13)
box(1438, 239, 112, 52, ["FCNet", "128 → 4 × 128", "ReLU · dropout"], INK, size=10, step=13)
line(1550, 150, 1588, 190, arrow="ink")
line(1550, 265, 1588, 232, arrow="ink")
for k in range(4):
    rect(1590 + k * 8, 160 + k * 6, 70, 70, PROT_L, PROT, rx=2, sw=1)
text(1650, 270, "Bilinear attention", 11.5, weight=700)
text(1650, 285, "4 heads · c ⊙ d", 10, fill=MUTED)
line(1690, 205, 1740, 205, arrow="ink")
vbar(1743, 160, 18, 90, F_C, 8)
text(1752, 268, "f ∈ ℝ²⁵⁶", 10.5, weight=700)
text(1752, 282, "interaction", 9.5, fill=MUTED)

# MLP
line(1752, 290, 1752, 320)
poly([(1752, 320), (1470, 320), (1470, 350)], arrow="ink")
layers = [(1470, 6), (1580, 4), (1690, 1)]
ys = {}
for lx, n in layers:
    top = 440 - (n - 1) * 18
    ys[lx] = [top + j * 36 for j in range(n)]
for (x1, _), (x2, _) in zip(layers, layers[1:]):
    for y1 in ys[x1]:
        for y2 in ys[x2]:
            line(x1, y1, x2, y2, "#C4CAD3", 0.7)
for lx, n in layers:
    for yy in ys[lx]:
        if n == 1:
            circle(lx, yy, 11, INK)
        else:
            circle(lx, yy, 10, "#FFFFFF", INK, 1.3)
text(1470, 556, "256", 10.5, weight=700)
text(1580, 556, "128", 10.5, weight=700)
text(1690, 556, "1", 10.5, weight=700)
text(1712, 430, "ŷ", 15, "start", 700)
text(1712, 450, "ln(IC50)", 11, "start")
text(1810, 405, "MLP", 11.5, weight=700)
lines(1810, 421, ["BatchNorm", "ReLU", "dropout 0.4"], 10, step=14, fill=MUTED)

# ================================================================ PANEL D
panel(1345, 600, 535, 300, "D", "Ensemble refinement", "gradient-boosted trees correct the network's residual error", "kept")
# C -> D: f and y-hat join one bus that drops into Z
poly([(1761, 205), (1860, 205), (1860, 672), (1458, 672), (1458, 682)], arrow="ink", sw=1.4)
poly([(1690, 452), (1690, 520), (1860, 520)], arrow=None, sw=1.4)
circle(1860, 520, 3.5, INK)
text(1868, 196, "f", 10.5, "start", 700)
text(1698, 536, "ŷ", 11, "start", 700)

vbar(1450, 684, 16, 100, F_C, 8)
rect(1450, 784, 16, 12, INK)
text(1458, 814, "Z = f ⊕ ŷ", 10.5, weight=700)
text(1458, 828, "257-d", 9.5, fill=MUTED)
line(1470, 740, 1520, 740, arrow="ink")


def tree(tx, ty):
    pts = [(tx, ty), (tx - 26, ty + 32), (tx + 26, ty + 32), (tx - 38, ty + 64), (tx - 14, ty + 64),
           (tx + 14, ty + 64), (tx + 38, ty + 64)]
    for a, b in [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]:
        line(*pts[a], *pts[b], MUTED, 1.1)
    for i, (px, py) in enumerate(pts):
        circle(px, py, 6, INK if i in (0, 4, 5) else "#B7BEC8")


tree(1570, 700)
text(1622, 740, "···", 16, weight=700)
tree(1676, 700)
text(1622, 800, "XGBoost", 12, weight=700)
text(1622, 816, "100 trees · depth 6 · lr 0.05", 10, fill=MUTED)
line(1725, 740, 1778, 740, arrow="ink")
add(f'<polygon points="1810,712 1838,740 1810,768 1782,740" fill="{OURS}"/>')
text(1810, 744, "IC50", 9.5, weight=700, fill="#FFFFFF")
text(1810, 790, "final", 10.5, weight=700)
text(1810, 804, "ln(IC50)", 10.5)

# ---------------------------------------------------------------- legend
lx = 20
lx += pill(lx, 930, "ours") + 8
text(lx, 934, "changed from, or added to, MoGraphDRP", 11, "start")
lx += 262
lx += pill(lx, 930, "kept") + 8
text(lx, 934, "kept from MoGraphDRP (hyperparameters from their Table 1)", 11, "start")
lx += 370
for label, col in [("cell line", CELL), ("drug", DRUG), ("protein / PPI", PROT)]:
    circle(lx + 6, 930, 6, col)
    text(lx + 18, 934, label, 11, "start")
    lx += 30 + len(label) * 6.5
text(W - 20, 934, "trained on 111,799 (cell line, drug) pairs · GDSC2 ln(IC50) · split grouped by cell line 70 / 15 / 15",
     10.5, "end", fill=MUTED)

add("</svg>")
OUT.write_text("\n".join(out), encoding="utf-8")
print(f"wrote {OUT}")
