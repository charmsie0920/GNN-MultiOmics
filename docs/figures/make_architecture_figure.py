"""Render the HeteroIC50GNN (E11) architecture figure as a standalone SVG.

Four panels in the style of the MoGraphDRP / AMOGEL overview figures:
(A) inputs, (B) heterogeneous graph, (C) encoders + HeteroConv, (D) pair head.
Every number drawn here comes from docs/architecture.md.

    python docs/figures/make_architecture_figure.py
"""

import random
from pathlib import Path

OUT = Path(__file__).with_name("hetero_ic50_gnn_architecture.svg")

INK, MUTED, PANEL, PANEL_STROKE = "#1E2530", "#5B6472", "#F4F6F9", "#D5DAE1"
CELL, CELL_L = "#2E6DB4", "#DCE8F6"
DRUG, DRUG_L = "#D9782A", "#FBE7D5"
PROT, PROT_L = "#3E9A6A", "#DDF0E5"
PPI = "#A9C9B5"
HI = "#C0392B"  # the 2-hop "shared evidence" path

out: list[str] = []


def add(s: str) -> None:
    out.append(s)


def rect(x, y, w, h, fill, stroke="none", rx=0, sw=1.0):
    add(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def circle(cx, cy, r, fill, stroke="none", sw=1.0):
    add(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')


def line(x1, y1, x2, y2, stroke=INK, sw=1.2, arrow=None, dash=None):
    m = f' marker-end="url(#ar-{arrow})"' if arrow else ""
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{stroke}" stroke-width="{sw}"{m}{d}/>')


def text(x, y, s, size=12, anchor="middle", weight=400, fill=INK, italic=False):
    it = ' font-style="italic"' if italic else ""
    add(f'<text x="{x}" y="{y}" font-size="{size}" text-anchor="{anchor}" font-weight="{weight}" fill="{fill}"{it}>{s}</text>')


def sub(base, s):
    return f'{base}<tspan dy="3" font-size="0.75em">{s}</tspan><tspan dy="-3"></tspan>'


def panel(x, w, letter, title, subtitle):
    rect(x, 45, w, 570, PANEL, PANEL_STROKE, rx=10)
    text(x + 16, 72, letter, 17, "start", 700)
    text(x + 36, 72, title, 14, "start", 700)
    text(x + 16, 90, subtitle, 11, "start", 400, MUTED)


def badge(x, y, s):
    circle(x, y, 9, "#FFFFFF", HI, 1.4)
    text(x, y + 3.5, s, 9, "middle", 700, HI)


W, H = 1600, 630
add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" '
    'aria-label="HeteroIC50GNN architecture: multi-omics and drug inputs become node features on a '
    'cell line, protein and drug graph; two HeteroConv layers refine them; the cell line and drug '
    'embeddings are concatenated and an MLP predicts ln IC50." '
    'font-family="Helvetica, Arial, sans-serif">')
add("<defs>")
for name, col in [("ink", INK), ("muted", MUTED), ("cell", CELL), ("drug", DRUG), ("prot", PROT)]:
    add(f'<marker id="ar-{name}" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="7" markerHeight="7" '
        f'orient="auto-start-reverse"><path d="M0,0 L8,4 L0,8 z" fill="{col}"/></marker>')
add("</defs>")
rect(0, 0, W, H, "#FFFFFF")

# --------------------------------------------------------------- A: inputs
panel(15, 355, "A", "Inputs", "per-modality preprocessing, fitted independently")
text(30, 115, "Cell-line multi-omics · 532 cell lines", 11.5, "start", 700, MUTED)
rng = random.Random(7)
modalities = [
    (75, "#2E6DB4", "Transcriptomics", "RNA-seq TPM", "log1p · scale"),
    (190, "#7A55B0", "Mut + CNV", "DepMap", "scale"),
    (305, "#2A8C9A", "Proteomics", "mass spec", "impute · scale"),
]
for cx, col, name, src, prep in modalities:
    x0, y0 = cx - 35, 124
    for r in range(6):
        for c in range(7):
            op = round(rng.uniform(0.12, 0.95), 2)
            add(f'<rect x="{x0 + c * 10}" y="{y0 + r * 10}" width="10" height="10" fill="{col}" fill-opacity="{op}"/>')
    rect(x0, y0, 70, 60, "none", col, sw=1)
    text(cx, 198, name, 11, weight=700)
    text(cx, 211, src, 10, fill=MUTED)
    line(cx, 216, cx, 229, MUTED, 1.2, "muted")
    rect(cx - 45, 231, 90, 34, "#FFFFFF", col, rx=5)
    text(cx, 245, prep, 10)
    text(cx, 259, "PCA → 128", 10.5, weight=700)
    line(cx, 265, cx, 278, MUTED, 1.2, "muted")
    rect(cx - 45, 280, 90, 14, col)
text(190, 314, sub("x", "c") + " ∈ ℝ<tspan dy=\"-5\" font-size=\"0.7em\">384</tspan><tspan dy=\"5\">  (concatenated)</tspan>", 11.5)
line(30, 330, 355, 330, PANEL_STROKE, 1)

text(30, 352, "Drug structure · 498 drugs", 11.5, "start", 700, MUTED)
hx, hy, hr = 62, 398, 19
pts = " ".join(f"{hx + hr * dx:.1f},{hy + hr * dy:.1f}" for dx, dy in
               [(0, -1), (0.866, -0.5), (0.866, 0.5), (0, 1), (-0.866, 0.5), (-0.866, -0.5)])
add(f'<polygon points="{pts}" fill="{DRUG_L}" stroke="{DRUG}" stroke-width="1.6"/>')
circle(hx, hy, 10, "none", DRUG, 1.2)
line(hx, hy - hr, hx, hy - hr - 12, DRUG, 1.6)
text(hx, 434, "SMILES", 10, fill=MUTED)
text(hx, 446, "PubChem", 10, fill=MUTED)
line(88, 398, 110, 398, MUTED, 1.2, "muted")
rect(112, 381, 104, 34, "#FFFFFF", DRUG, rx=5)
text(164, 395, "RDKit Morgan", 10.5, weight=700)
text(164, 409, "radius 2", 10)
line(216, 398, 232, 398, MUTED, 1.2, "muted")
for i, b in enumerate("101100100101"):
    rect(234 + i * 9, 390, 9, 16, DRUG if b == "1" else "#FFFFFF", DRUG, sw=0.8)
text(288, 426, sub("x", "d") + " ∈ {0,1}<tspan dy=\"-5\" font-size=\"0.7em\">2048</tspan>", 11.5)
line(30, 460, 355, 460, PANEL_STROKE, 1)

text(30, 482, "Protein interaction network", 11.5, "start", 700, MUTED)
net = [(45, 510), (82, 500), (108, 528), (66, 540), (96, 566), (40, 568), (70, 590)]
for a, b in [(0, 1), (0, 3), (1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6), (5, 6), (3, 6)]:
    line(*net[a], *net[b], PPI, 1.4)
for x, y in net:
    circle(x, y, 6, PROT_L, PROT, 1.3)
for i, s in enumerate(["STRING v12.0, human", "combined_score ≥ 700", "16,214 proteins", "473,860 edges"]):
    text(135, 512 + i * 19, s, 11, "start", 700 if i == 0 else 400)

# ------------------------------------------------------------ B: the graph
panel(385, 360, "B", "Heterogeneous graph 𝒢", "built once offline (~20 MB); red = the 2-hop path in C")
text(440, 118, "cell_line", 12, weight=700, fill=CELL)
text(578, 118, "protein", 12, weight=700, fill=PROT)
text(705, 118, "drug", 12, weight=700, fill=DRUG)
cells = [(440, y) for y in (160, 240, 320, 400)]
prots = [(555, 145), (612, 185), (548, 228), (605, 272), (545, 318), (610, 360), (560, 402)]
drugs = [(705, y) for y in (180, 272, 370)]
ppi = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5), (4, 6), (5, 6)]
mut = [(0, 0), (0, 2), (1, 2), (2, 4), (3, 4), (3, 6)]
tgt = [(0, 1), (1, 3), (1, 5), (2, 5)]
for a, b in ppi:
    line(*prots[a], *prots[b], PPI, 1.8)
for c, p in mut:
    line(*cells[c], *prots[p], CELL, 1.3)
for d, p in tgt:
    line(*drugs[d], *prots[p], DRUG, 1.3)
line(*cells[1], *prots[2], HI, 3)
line(*prots[2], *prots[3], HI, 3)
line(*drugs[1], *prots[3], HI, 3)
for x, y in cells:
    circle(x, y, 12, CELL_L, CELL, 1.8)
for x, y in prots:
    circle(x, y, 9, PROT_L, PROT, 1.6)
for x, y in drugs:
    rect(x - 11, y - 11, 22, 22, DRUG_L, DRUG, rx=4, sw=1.8)
badge(494, 222, "ℓ1")
badge(590, 238, "ℓ2")
badge(655, 285, "ℓ1")
text(490, 440, "has_mutation", 10.5, fill=CELL, italic=True)
text(578, 440, "interacts_with", 10.5, fill=PROT, italic=True)
text(665, 440, "targets", 10.5, fill=DRUG, italic=True)

text(400, 474, "Nodes", 11, "start", 700, MUTED)
rows = [
    ("cell_line", "532", "384-d omics features", CELL, "c"),
    ("protein", "16,214", "no input features", PROT, "c"),
    ("drug", "498", "2048-bit fingerprint", DRUG, "s"),
]
for i, (n, k, d, col, shape) in enumerate(rows):
    y = 492 + i * 17
    if shape == "c":
        circle(405, y - 4, 5, col)
    else:
        rect(400, y - 9, 10, 10, col, rx=2)
    text(416, y, n, 11, "start", 700, col)
    text(538, y, k, 11, "end")
    text(548, y, d, 10.5, "start", fill=MUTED)
text(400, 552, "Relations (each also mirrored → 5 total)", 11, "start", 700, MUTED)
rels = [("has_mutation", "4,341", "cancer-driver mutations", CELL),
        ("targets", "683", "GDSC drug targets", DRUG),
        ("interacts_with", "473,860", "STRING PPI", PROT)]
for i, (n, k, d, col) in enumerate(rels):
    y = 570 + i * 17
    text(400, y, n, 11, "start", 700, col)
    text(538, y, k, 11, "end")
    text(548, y, d, 10.5, "start", fill=MUTED)

# ------------------------------------------- C: encoders + message passing
panel(765, 440, "C", "Encoders + relational message passing", "all node types projected to D = 128")
text(865, 124, "Node encoders", 11.5, weight=700, fill=MUTED)
enc = [(140, CELL, CELL_L, "Linear 384 → 128", "+ ReLU"),
       (250, PROT, PROT_L, "Embedding", "16,214 × 128"),
       (360, DRUG, DRUG_L, "Linear 2048 → 128", "+ ReLU")]
for y, col, lt, a, b in enc:
    rect(790, y, 150, 44, lt, col, rx=6, sw=1.4)
    text(865, y + 19, a, 11, weight=700)
    text(865, y + 35, b, 10.5)
    line(940, y + 22, 966, y + 22, INK, 1.3, "ink")
for i, s in enumerate(["Proteins have no features,", "so each gets a learnable", "ID vector: 74% of the", "2.82 M parameters"]):
    text(865, 440 + i * 16, s, 10.5, fill=MUTED)

rect(983, 108, 208, 318, "#FFFFFF", PANEL_STROKE, rx=8)
rect(970, 96, 208, 318, "#FFFFFF", INK, rx=8, sw=1.3)
text(1191, 446, "× 2 layers", 12, "end", 700)
text(1074, 118, "HeteroConv layer ℓ · 5 × SAGEConv", 11.5, weight=700)
convs = [(128, "interacts_with", "protein → protein", PROT, PROT_L),
         (166, "targets", "drug → protein", DRUG, DRUG_L),
         (204, "has_mutation", "cell_line → protein", CELL, CELL_L),
         (262, "rev_targets", "protein → drug", PROT, PROT_L),
         (320, "rev_has_mutation", "protein → cell_line", PROT, PROT_L)]
for y, rel, path, col, lt in convs:
    rect(980, y, 140, 32, lt, col, rx=4)
    text(1050, y + 13, rel, 10, weight=700, italic=True)
    text(1050, y + 26, path, 9.5, fill=MUTED)
circle(1140, 182, 11, "#FFFFFF", INK, 1.3)
text(1140, 187, "Σ", 13, weight=700)
for y in (144, 182, 220):
    line(1120, y, 1128, y + (182 - y) * 0.7, INK, 1, "ink")
line(1151, 182, 1158, 182, INK, 1.2, "ink")
circle(1166, 182, 6, PROT)
line(1120, 278, 1156, 278, INK, 1.2, "ink")
rect(1160, 272, 12, 12, DRUG, rx=2)
line(1120, 336, 1157, 336, INK, 1.2, "ink")
circle(1166, 336, 6, CELL)
for i, s in enumerate(["sum over relations per", "destination type, then ReLU;", "independent weights per relation"]):
    text(1074, 374 + i * 14, s, 9.5, fill=MUTED)

tx = 985
text(tx, 478, "What two layers buy", 11, "start", 700, MUTED)
badge(tx + 8, 497, "ℓ1")
text(tx + 22, 501, "mutation + target signal", 10.5, "start")
text(tx + 22, 515, "move onto protein nodes", 10.5, "start")
badge(tx + 8, 535, "ℓ2")
text(tx + 22, 539, "one PPI hop, then back out", 10.5, "start")
text(tx + 22, 553, "to cell_line and drug", 10.5, "start")
text(tx, 580, "Full-batch: whole graph every step;", 10.5, "start", fill=MUTED)
text(tx, 594, "only label pairs are mini-batched", 10.5, "start", fill=MUTED)

# ------------------------------------------------------------ D: pair head
panel(1225, 360, "D", "Pair prediction", "for each (cell line i, drug j) pair in a batch of 1,024")
for i in range(8):
    rect(1245, 132 + i * 14, 16, 14, CELL_L if i % 2 else "#BBD1EC", CELL, sw=0.8)
    rect(1245, 300 + i * 14, 16, 14, DRUG_L if i % 2 else "#F4C9A2", DRUG, sw=0.8)
text(1253, 262, sub("h", "c") + "[i]", 11, weight=700, fill=CELL)
text(1253, 276, "128", 10, fill=MUTED)
text(1253, 430, sub("h", "d") + "[j]", 11, weight=700, fill=DRUG)
text(1253, 444, "128", 10, fill=MUTED)
for i in range(16):
    col, lt = (CELL, "#BBD1EC") if i < 8 else (DRUG, "#F4C9A2")
    rect(1302, 175 + i * 12, 16, 12, lt, col, sw=0.8)
line(1263, 188, 1298, 215, INK, 1.2, "ink")
line(1263, 356, 1298, 330, INK, 1.2, "ink")
text(1310, 166, "concat", 10.5, weight=700)
text(1310, 386, sub("z", "ij"), 11, weight=700)
text(1310, 400, "256", 10, fill=MUTED)

c1 = [(1370, 190 + i * 26.7) for i in range(7)]
c2 = [(1430, 215 + i * 27.5) for i in range(5)]
o = (1490, 270)
for i, (x, y) in enumerate(c1):
    line(1320, 180 + i * 29.5, x - 8, y, "#C4CAD3", 0.8)
for a in c1:
    for b in c2:
        line(*a, *b, "#C4CAD3", 0.7)
for b in c2:
    line(*b, *o, "#C4CAD3", 0.8)
for x, y in c1 + c2:
    circle(x, y, 8, "#FFFFFF", INK, 1.3)
circle(*o, 9, INK)
text(1506, 267, "ŷ", 15, "start", 700, HI)
text(1506, 283, "ln(IC50)", 10.5, "start", 700)
for x, n in [(1370, "256"), (1430, "128"), (1490, "1")]:
    text(x, 392, n, 10.5, weight=700)
text(1432, 420, "hidden layer = Linear → BatchNorm", 10.5, fill=MUTED)
text(1432, 434, "→ ReLU → Dropout 0.3", 10.5, fill=MUTED)

rect(1240, 460, 330, 142, "#FFFFFF", PANEL_STROKE, rx=6)
text(1254, 481, "Training", 11.5, "start", 700)
train = [("Loss", "MSE against GDSC2 ln(IC50)"),
         ("Optimiser", "AdamW, lr 1e-3, weight decay 1e-4"),
         ("Batch", "1,024 label pairs"),
         ("Split", "70 / 15 / 15, grouped by cell line"),
         ("Checkpoint", "best validation RMSE")]
for i, (k, v) in enumerate(train):
    y = 502 + i * 20
    text(1254, y, k, 10.5, "start", 700, MUTED)
    text(1334, y, v, 10.5, "start")

# chevrons between panels
for x in (372, 752, 1212):
    add(f'<polygon points="{x - 4},318 {x + 6},330 {x - 4},342" fill="{MUTED}"/>')

add("</svg>")
OUT.write_text("\n".join(out), encoding="utf-8")
print(f"wrote {OUT}")
