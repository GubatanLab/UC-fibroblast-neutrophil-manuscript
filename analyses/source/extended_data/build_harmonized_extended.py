from pathlib import Path
import hashlib
import json
import shutil
import sys

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parent
PKG = ROOT / "tmp/canonical_nature_20260903/packages"
sys.path.insert(0, str(PKG))

import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

MM = 72 / 25.4
PAGE_W_MM = 180
FONT = "C:/Windows/Fonts/arial.ttf"
BOLD = "C:/Windows/Fonts/arialbd.ttf"

ACTIVE = ROOT / "output/Archive_before_extended_data_harmonization_2026-09-03/Extended Data"
OUT = ROOT / "output/Nature_Extended_Data_Manuscript_Matched_2026-09-03"
FIG_OUT = OUT / "extended_data"
COMPONENTS = WORK / "components"
QA = WORK / "qa"
TABLES = ROOT / "output/additional_analyses_2026-09-02/tables"
OLD = ROOT / "output/Archive_before_canonical_promotion_2026-09-03/Supplementary Figures/Ordered Individual Figures"

HUM_TABLES = Path("input_data/human/UC_fibroblast_neutrophil_analysis/tables")
COCULTURE = Path("input_data/mouse/Neutrophil_Fibroblast_Coculture_Blockade_Figures/scVI_Harmony_All_Coculture")
NEUT_SOURCE = COCULTURE / "Neutrophil_Subclustering/source_data"

for directory in [FIG_OUT, COMPONENTS, QA, OUT / "reports", OUT / "source_data"]:
    directory.mkdir(parents=True, exist_ok=True)

# Manuscript-wide visual language.
INK = "#2F3337"
MUTED = "#6B7280"
GRID = "#E5E7EB"
CONTROL = "#3B75B9"
NONINFLAMED = "#C9A227"
INFLAMED = "#D44B50"
DSS = "#E76F00"
BLOCKADE = "#1FA987"
STATE = {"CXCR4": "#6F4C9B", "MX1/ISG": "#3B75B9", "OSM": "#D44B50", "PADI4": "#E6862B"}
DIVERGING = LinearSegmentedColormap.from_list("manuscript_diverging", [CONTROL, "#F7F7F7", "#B53A45"])

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 6.5,
        "axes.titlesize": 7,
        "axes.labelsize": 6.5,
        "xtick.labelsize": 5.8,
        "ytick.labelsize": 5.8,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "lines.linewidth": 0.75,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "savefig.facecolor": "white",
    }
)


def mm_rect(values):
    return fitz.Rect(*(v * MM for v in values))


def spans(page):
    return [
        span
        for block in page.get_text("dict")["blocks"]
        if "lines" in block
        for line in block["lines"]
        for span in line["spans"]
    ]


def text(page, x_mm, y_mm, value, size=6.5, bold=False, color=INK):
    name = "ArialB" if bold else "Arial"
    page.insert_font(fontname=name, fontfile=BOLD if bold else FONT)
    rgb = tuple(int(color[i : i + 2], 16) / 255 for i in (1, 3, 5)) if isinstance(color, str) else color
    page.insert_text((x_mm * MM, y_mm * MM), str(value), fontname=name, fontsize=size, color=rgb)


def panel_header(page, x_mm, y_mm, label, title, subtitle=None, title_size=8):
    """y_mm is the top of a consistent 9-mm header band."""
    text(page, x_mm, y_mm + 3.2, label, 8, True)
    text(page, x_mm + 6.2, y_mm + 3.2, title, title_size, True)
    if subtitle:
        text(page, x_mm + 6.2, y_mm + 7.0, subtitle, 6, False, MUTED)


def new_page(height_mm):
    doc = fitz.open()
    page = doc.new_page(width=PAGE_W_MM * MM, height=height_mm * MM)
    return doc, page


def tight(doc, clip=None, pad=1.5, threshold=244):
    page = doc[0]
    region = fitz.Rect(clip) if clip else page.rect
    region &= page.rect
    pix = page.get_pixmap(matrix=fitz.Matrix(1.1, 1.1), clip=region, colorspace=fitz.csRGB, alpha=False)
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
    yy, xx = np.where(arr.min(axis=2) < threshold)
    if not len(xx):
        return region
    scale = 1.1
    found = fitz.Rect(
        region.x0 + xx.min() / scale - pad,
        region.y0 + yy.min() / scale - pad,
        region.x0 + (xx.max() + 1) / scale + pad,
        region.y0 + (yy.max() + 1) / scale + pad,
    )
    return found & region


def place(page, source_doc, clip, box_mm, tight_crop=True, align="center"):
    source = tight(source_doc, clip) if tight_crop else (fitz.Rect(clip) if clip else source_doc[0].rect)
    target = mm_rect(box_mm)
    scale = min(target.width / source.width, target.height / source.height)
    width, height = source.width * scale, source.height * scale
    dx = (target.width - width) / 2
    dy = (target.height - height) / 2
    if align == "top":
        dy = 0
    dest = fitz.Rect(target.x0 + dx, target.y0 + dy, target.x0 + dx + width, target.y0 + dy + height)
    page.show_pdf_page(dest, source_doc, 0, clip=source, keep_proportion=True)
    return dest


def current(n):
    return fitz.open(ACTIVE / f"Extended_Data_Figure_{n:02d}.pdf")


def old(n):
    return fitz.open(OLD / f"Supplementary_Figure_S{n:02d}.pdf")


def component(fig, name):
    path = COMPONENTS / f"{name}.pdf"
    fig.savefig(path, transparent=False)
    plt.close(fig)
    return fitz.open(path)


inventory = []


def save_figure(doc, number, sources, height_mm, changes):
    stem = FIG_OUT / f"Extended_Data_Figure_{number:02d}"
    doc.set_metadata(
        {
            "title": f"Extended Data Fig. {number}",
            "author": "UC fibroblast-neutrophil study",
            "subject": "Nature-style harmonized figure; observations and statistical annotations preserved",
            "keywords": "Extended Data; Nature; harmonized; no FACS reanalysis",
        }
    )
    doc.save(stem.with_suffix(".pdf"), garbage=4, deflate=True)

    # High-resolution review PNG and production TIFF. Fall back to 300 dpi only if LZW exceeds 10 MB.
    pix = doc[0].get_pixmap(matrix=fitz.Matrix(450 / 72, 450 / 72), alpha=False)
    pix.save(stem.with_suffix(".png"))
    with Image.open(stem.with_suffix(".png")) as image:
        image.save(stem.with_suffix(".tiff"), compression="tiff_lzw", dpi=(450, 450))
    if stem.with_suffix(".tiff").stat().st_size > 10_000_000:
        pix = doc[0].get_pixmap(matrix=fitz.Matrix(300 / 72, 300 / 72), alpha=False)
        image = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
        image.save(stem.with_suffix(".tiff"), compression="tiff_lzw", dpi=(300, 300))
    doc[0].get_pixmap(matrix=fitz.Matrix(2.6, 2.6), alpha=False).save(QA / f"ED{number:02d}.png")

    inventory.append(
        {
            "figure": number,
            "pages": len(doc),
            "width_mm": PAGE_W_MM,
            "height_mm": height_mm,
            "pdf_bytes": stem.with_suffix(".pdf").stat().st_size,
            "tiff_bytes": stem.with_suffix(".tiff").stat().st_size,
            "sources": sources,
            "changes": changes,
            "scientific_content": "unchanged; no new experiments or FACS reanalysis",
        }
    )


# -----------------------------------------------------------------------------
# Extended Data Figure 1: fully native vector redraw from exported tables.
# -----------------------------------------------------------------------------
fib = pd.read_csv(HUM_TABLES / "09_fibroblast_biopsy_programs.tsv", sep="\t")
programs = [
    ("Alpha5Beta1_adhesion", "α5β1 adhesion"),
    ("ECM_remodeling", "ECM remodeling"),
    ("FAP_inflammatory", "FAP inflammatory"),
    ("NFkB_AP1", "NF-κB/AP-1"),
    ("Neutrophil_recruitment", "Neutrophil recruitment"),
    ("OSM_response", "OSM response"),
    ("TGFb_response", "TGF-β response"),
    ("YAP_mechanotransduction", "YAP mechanotransduction"),
]
paired_ids = [pid for pid, frame in fib.groupby("PatientID") if {"Inflamed UC", "Uninflamed UC"}.issubset(set(frame.condition))]

fig, axes = plt.subplots(2, 4, figsize=(174 / 25.4, 72 / 25.4))
fig.subplots_adjust(left=0.055, right=0.995, top=0.95, bottom=0.16, wspace=0.32, hspace=0.62)
for index, (column, title) in enumerate(programs):
    ax = axes.flat[index]
    for pid in paired_ids:
        frame = fib[fib.PatientID.eq(pid)].set_index("condition")
        values = [frame.loc["Uninflamed UC", column], frame.loc["Inflamed UC", column]]
        ax.plot([0, 1], values, color="#C7CBD1", lw=0.65, zorder=1)
    left = [fib[(fib.PatientID.eq(pid)) & (fib.condition.eq("Uninflamed UC"))][column].iloc[0] for pid in paired_ids]
    right = [fib[(fib.PatientID.eq(pid)) & (fib.condition.eq("Inflamed UC"))][column].iloc[0] for pid in paired_ids]
    ax.scatter(np.zeros(len(left)), left, s=11, color=NONINFLAMED, edgecolor="white", linewidth=0.25, zorder=3)
    ax.scatter(np.ones(len(right)), right, s=11, color=INFLAMED, edgecolor="white", linewidth=0.25, zorder=3)
    ax.set_xlim(-0.22, 1.22)
    ax.set_xticks([0, 1], ["Noninflamed", "Inflamed"])
    ax.tick_params(axis="x", length=0, pad=1)
    ax.set_title(title, loc="left", fontweight="bold", pad=3)
    ax.grid(axis="y", color=GRID, lw=0.45)
    ax.set_axisbelow(True)
    if index in (0, 4):
        ax.set_ylabel("Mean module expression")
    else:
        ax.set_ylabel("")
native_ed1a = component(fig, "ED1_paired_programs")

neut = pd.read_csv(HUM_TABLES / "22_neutrophil_state_program_summary.tsv", sep="\t")
prog_cols = [c for c in neut.columns if c not in ("neut_state", "n")]
matrix = neut[prog_cols].astype(float).to_numpy()
matrix = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0)
row_labels = [
    value.replace("Neutrophil ", "").replace("Monocyte-derived macrophage", "Monocyte-derived\nmacrophage")
    for value in neut.neut_state
]
column_labels = [
    "Recruitment /\nmigration",
    "Retention /\naging",
    "Degranulation",
    "NET-associated",
    "Oxidative\nburst",
    "Inflammatory",
    "Interferon",
    "Survival /\nimmaturity",
    "Tissue injury",
]
fig, ax = plt.subplots(figsize=(174 / 25.4, 45 / 25.4))
fig.subplots_adjust(left=0.22, right=0.89, top=0.97, bottom=0.25)
im = ax.imshow(matrix, cmap=DIVERGING, vmin=-1.8, vmax=1.8, aspect="auto")
ax.set_xticks(range(len(column_labels)), column_labels, rotation=42, ha="right")
ax.set_yticks(range(len(row_labels)), row_labels)
ax.tick_params(length=0)
ax.set_xticks(np.arange(-0.5, len(column_labels), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=0.7)
for spine in ax.spines.values():
    spine.set_visible(False)
cax = fig.add_axes([0.92, 0.29, 0.015, 0.58])
cb = fig.colorbar(im, cax=cax)
cb.set_label("Program z-score", fontsize=6)
cb.ax.tick_params(labelsize=5.5, width=0.4)
native_ed1b = component(fig, "ED1_state_heatmap")

doc, page = new_page(145)
panel_header(page, 3, 2, "a", "Paired fibroblast program changes", f"Matched noninflamed and inflamed UC biopsies; lines connect {len(paired_ids)} patients")
place(page, native_ed1a, None, (3, 12, 177, 84), False)
panel_header(page, 3, 88, "b", "Neutrophil functional programs by state", "Column-wise z-scores; non-neutrophil and low-RNA labels retained as annotation controls")
place(page, native_ed1b, None, (3, 98, 177, 143), False)
save_figure(doc, 1, [str(HUM_TABLES / "09_fibroblast_biopsy_programs.tsv"), str(HUM_TABLES / "22_neutrophil_state_program_summary.tsv")], 145, "Native paired plots and heatmap; manuscript color semantics; readable biological labels; compact two-tier layout.")


# -----------------------------------------------------------------------------
# Extended Data Figure 2: CODEX atlas plus a native validation heatmap.
# -----------------------------------------------------------------------------
phenotype = np.array(
    [
        [-1.2, 0.7, 1.1, 0.4, -0.8, -1.3, -1.6, -1.2],
        [1.8, -0.9, -1.0, -0.7, -1.4, 0.9, 1.3, 1.7],
        [-0.2, 0.3, 1.3, 0.3, 1.3, -0.2, -0.3, -0.3],
        [0.2, -1.4, -0.9, -1.5, 0.7, -0.7, 0.8, 0.3],
        [-0.5, 1.3, -0.5, 1.5, 0.2, 1.3, -0.2, -0.6],
    ]
)
fig, ax = plt.subplots(figsize=(72 / 25.4, 56 / 25.4))
fig.subplots_adjust(left=0.26, right=0.88, top=0.96, bottom=0.22)
im = ax.imshow(phenotype, cmap=DIVERGING, vmin=-1.8, vmax=1.8, aspect="auto")
ax.set_xticks(range(8), ["PADI4", "CXCR4", "MX1", "OSM", "CD16", "CD11b", "CD15", "CD66b"], rotation=47, ha="right")
ax.set_yticks(range(5), ["Unpolarized", "PADI4+", "MX1+", "CXCR4+", "OSM+"])
ax.tick_params(length=0)
for i in range(phenotype.shape[0]):
    for j in range(phenotype.shape[1]):
        ax.text(j, i, f"{phenotype[i, j]:.1f}", ha="center", va="center", fontsize=5.7, color="white" if abs(phenotype[i, j]) > 1.0 else INK)
ax.set_xticks(np.arange(-0.5, 8, 1), minor=True)
ax.set_yticks(np.arange(-0.5, 5, 1), minor=True)
ax.grid(which="minor", color="white", lw=0.6)
for spine in ax.spines.values():
    spine.set_visible(False)
cax = fig.add_axes([0.91, 0.29, 0.022, 0.55])
cb = fig.colorbar(im, cax=cax)
cb.set_label("Marker z-score", fontsize=6)
cb.ax.tick_params(labelsize=5.5, width=0.4)
native_ed2b = component(fig, "ED2_phenotype_heatmap")

src2 = current(2)
doc, page = new_page(137)
panel_header(page, 3, 2, "a", "Broad CODEX cell atlas", "245,020 cells; 12 broad annotations")
place(page, src2, mm_rect((3, 13, 99, 88)), (3, 12, 99, 78), True)
panel_header(page, 102, 2, "b", "Neutrophil phenotype validation", "Median marker expression, z-scored across states")
place(page, native_ed2b, None, (104, 12, 177, 78), False)
panel_header(page, 3, 82, "c", "Region-resolved cell composition", "Six control, nine noninflamed UC and nine inflamed UC region identifiers")
text(page, 30, 94.5, "Control", 5.8, True, CONTROL)
text(page, 78, 94.5, "UC noninflamed", 5.8, True, NONINFLAMED)
text(page, 139, 94.5, "UC inflamed", 5.8, True, INFLAMED)
place(page, src2, mm_rect((3, 105, 177, 147)), (3, 97, 177, 135), True)
save_figure(doc, 2, ["Prior canonical CODEX UMAP and region-composition source vectors", "Original displayed phenotype heatmap values"], 137, "Rebuilt asymmetric atlas layout; native phenotype heatmap; consistent cohort colors and hierarchy; no cell-level reanalysis.")


# -----------------------------------------------------------------------------
# Extended Data Figure 3: spatial specificity and current adjusted inference.
# -----------------------------------------------------------------------------
spatial = pd.read_csv(TABLES / "spatial_UC_pooled_effects.csv")
restricted = pd.read_csv(TABLES / "spatial_restricted_label_tests.csv")
forest = spatial[spatial.model.eq("Density + anatomy adjusted")].merge(restricted[["marker", "q_restricted"]], on="marker")
order = ["CXCR4", "OSM", "CD16", "CD11b"]
forest["marker"] = pd.Categorical(forest.marker, order, ordered=True)
forest = forest.sort_values("marker")
fig, ax = plt.subplots(figsize=(174 / 25.4, 31 / 25.4))
fig.subplots_adjust(left=0.11, right=0.61, top=0.95, bottom=0.28)
for y, row in enumerate(forest.itertuples()):
    significant = row.q_restricted < 0.05
    ax.errorbar(
        row.effect,
        y,
        xerr=[[row.effect - row.low], [row.high - row.effect]],
        fmt="o",
        ms=3.8,
        color=CONTROL,
        mfc=CONTROL if significant else "white",
        mew=0.9,
        capsize=2,
    )
    yy = 0.90 - (y + 0.5) * (0.72 / len(forest))
    fig.text(0.73, yy, f"{row.q:.3f}", ha="center", va="center", fontsize=6.2)
    fig.text(0.91, yy, f"{row.q_restricted:.3f}", ha="center", va="center", fontsize=6.2, fontweight="bold" if significant else "normal")
ax.set_yticks(range(len(forest)), forest.marker.astype(str))
ax.set_ylim(len(forest) - 0.5, -0.5)
ax.axvline(0, color="#9CA3AF", ls=(0, (2, 2)), lw=0.7)
ax.set_xlabel("Adjusted proximity coefficient (95% CI)")
ax.grid(axis="x", color=GRID, lw=0.45)
fig.text(0.73, 0.94, "Model q", ha="center", fontsize=6.5, fontweight="bold")
fig.text(0.91, 0.94, "Spatial-null q", ha="center", fontsize=6.5, fontweight="bold")
native_ed3c = component(fig, "ED3_spatial_forest")

src3 = current(3)
doc, page = new_page(170)
panel_header(page, 3, 2, "a", "Immune-lineage spatial specificity", "Median log2 observed/expected proximity across FAP/α5β1 strata")
place(page, src3, mm_rect((3, 13, 177, 69)), (3, 12, 177, 59), True)
panel_header(page, 3, 63, "b", "Threshold sensitivity", "FAP and α5β1 percentile cutoffs across tissue conditions")
place(page, src3, mm_rect((3, 80, 177, 135)), (3, 73, 177, 121), True)
panel_header(page, 3, 125, "c", "Adjusted neutrophil proximity effects", "Filled point: restricted spatial-null q < 0.05; bars show 95% CI")
place(page, native_ed3c, None, (3, 135, 177, 169), False)
save_figure(doc, 3, ["Prior canonical spatial heatmaps", str(TABLES / "spatial_UC_pooled_effects.csv"), str(TABLES / "spatial_restricted_label_tests.csv")], 170, "Common diverging palette and headers; current adjusted forest plot redrawn natively; compact three-tier reading order.")


# -----------------------------------------------------------------------------
# Extended Data Figure 4: fibroblast annotation audit.
# -----------------------------------------------------------------------------
src4 = current(4)
doc, page = new_page(153)
panel_header(page, 3, 2, "a", "Fibroblast subclusters and annotations", "Donor-harmonized fibroblast embedding; annotation resolution 0.4")
place(page, src4, mm_rect((3, 19, 108, 77)), (3, 12, 109, 75), True)
panel_header(page, 113, 2, "b", "Target transcripts by condition", "FAP, ITGA5 and ITGB1 expression")
place(page, src4, mm_rect((115, 32, 177, 77)), (113, 12, 177, 75), True)
panel_header(page, 3, 80, "c", "Marker programs supporting annotations", "Dot size: percent expressed; color: scaled mean expression")
place(page, src4, mm_rect((3, 98, 177, 167)), (3, 90, 177, 151), True)
save_figure(doc, 4, ["Prior canonical fibroblast UMAP, target-expression and marker-program source vectors"], 153, "Recomposed source vectors under one typographic hierarchy; enlarged top panels; shortened page and removed figure-level dead space.")


# -----------------------------------------------------------------------------
# Extended Data Figure 5: neutrophil subcluster annotation and cleanup audit.
# -----------------------------------------------------------------------------
marker = pd.read_csv(NEUT_SOURCE / "neutrophil_annotation_marker_dotplot_source.csv")
genes = marker.gene.drop_duplicates().tolist()
annotations = marker.annotation.drop_duplicates().tolist()
short = {
    "Low-signal neutrophil": "Low signal",
    "Activated OSM/CXCR4 neutrophil": "Activated OSM/CXCR4",
    "IL1R2+ inflammatory neutrophil": "IL1R2+ inflammatory",
    "DHFR+ proliferative/stress neutrophil": "DHFR+ stress/proliferative",
    "PADI4/translation-high neutrophil": "PADI4 / translation high",
    "LTF/BPI immature neutrophil": "LTF/BPI immature",
    "CCL3/CCL4 inflammatory neutrophil": "CCL3/CCL4 inflammatory",
    "RORA+ atypical neutrophil": "RORA+ atypical",
}
fig, ax = plt.subplots(figsize=(174 / 25.4, 51 / 25.4))
fig.subplots_adjust(left=0.25, right=0.82, top=0.97, bottom=0.27)
xmap, ymap = {v: i for i, v in enumerate(genes)}, {v: i for i, v in enumerate(annotations)}
scatter = ax.scatter(
    marker.gene.map(xmap),
    marker.annotation.map(ymap),
    s=marker.pct_expressed * 0.52,
    c=marker.scaled_expression,
    cmap=DIVERGING,
    vmin=-2,
    vmax=2,
    edgecolors="none",
)
ax.set_xticks(range(len(genes)), genes, rotation=58, ha="right")
ax.set_yticks(range(len(annotations)), [short[v] for v in annotations])
ax.set_ylim(len(annotations) - 0.5, -0.5)
ax.set_xlim(-0.5, len(genes) - 0.5)
ax.tick_params(length=0)
ax.grid(color=GRID, lw=0.38)
ax.set_axisbelow(True)
cax = fig.add_axes([0.85, 0.57, 0.016, 0.31])
cb = fig.colorbar(scatter, cax=cax)
cb.set_label("Scaled mean", fontsize=6)
cb.ax.tick_params(labelsize=5.3, width=0.4)
for value in [25, 50, 100]:
    ax.scatter([], [], s=value * 0.52, color="#4B5563", label=str(value))
ax.legend(title="Detected (%)", loc="lower left", bbox_to_anchor=(1.03, -0.02), frameon=False, fontsize=5.5, title_fontsize=6, labelspacing=0.7)
native_ed5b = component(fig, "ED5_marker_dotplot")

cleanup = pd.read_csv(NEUT_SOURCE / "neutrophil_cleanup_audit.csv")
cleanup = cleanup.groupby(["cleanup_action", "cleanup_reason"], as_index=False)["cells"].sum()
cleanup["label"] = cleanup.cleanup_reason.replace(
    {
        "Retained after cell-state cleanup": "Retained neutrophils",
        "Low-signal/unresolved neutrophil state": "Low-signal / unresolved state",
        "Doublet, contaminant-like, or singleton/outlier community": "Doublet / contaminant / outlier",
    }
)
cleanup = cleanup.sort_values("cells", ascending=True)
fig, ax = plt.subplots(figsize=(174 / 25.4, 22 / 25.4))
fig.subplots_adjust(left=0.31, right=0.92, top=0.96, bottom=0.28)
colors = [BLOCKADE if action == "Retained" else "#9CA3AF" for action in cleanup.cleanup_action]
ax.barh(range(len(cleanup)), cleanup.cells, color=colors, height=0.62)
ax.set_yticks(range(len(cleanup)), cleanup.label)
ax.set_xlabel("")
ax.set_xlim(0, cleanup.cells.max() * 1.13)
for y, value in enumerate(cleanup.cells):
    ax.text(value + cleanup.cells.max() * 0.018, y, f"{value:,}", va="center", fontsize=6)
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.grid(axis="x", color=GRID, lw=0.4)
ax.set_axisbelow(True)
native_ed5c = component(fig, "ED5_cleanup")

src5 = current(5)
doc, page = new_page(168)
panel_header(page, 3, 2, "a", "Neutrophil subcluster atlas", "Retained cells colored by transcriptomic subcluster annotation")
place(page, src5, mm_rect((3, 12, 177, 71)), (3, 12, 177, 67), True)
# Small orientation axes, matching the main UMAPs.
page.draw_line((29 * MM, 61 * MM), (38 * MM, 61 * MM), color=(0.25, 0.27, 0.29), width=0.55)
page.draw_line((29 * MM, 61 * MM), (29 * MM, 52 * MM), color=(0.25, 0.27, 0.29), width=0.55)
text(page, 38.5, 61.5, "UMAP 1", 5.2, False, MUTED)
text(page, 27.0, 51.5, "UMAP 2", 5.2, False, MUTED)
panel_header(page, 3, 70, "b", "Marker support for subcluster annotations", "Dot size: detected cells; color: scaled mean expression")
place(page, native_ed5b, None, (3, 80, 177, 126), False)
panel_header(page, 3, 129, "c", "Cell-retention and cleanup inventory", "Original counts summarized by disposition")
place(page, native_ed5c, None, (3, 137, 177, 166), False)
text(page, 88, 166.5, "Cells", 6, False, INK)
save_figure(doc, 5, ["Prior canonical neutrophil UMAP", str(NEUT_SOURCE / "neutrophil_annotation_marker_dotplot_source.csv"), str(NEUT_SOURCE / "neutrophil_cleanup_audit.csv")], 168, "UMAP enlarged and oriented; marker and cleanup panels redrawn natively; consistent state-color and typography system.")


# -----------------------------------------------------------------------------
# Extended Data Figure 6: donor-aware effects and OSM sensitivity.
# -----------------------------------------------------------------------------
composition = pd.read_csv(TABLES / "composition_donor_contrasts.csv")
composition = composition[composition.min_cells.eq(0)]
contrast_order = [
    "control_fibroblasts_vs_alone",
    "UC_vs_alone",
    "UC_vs_control_fibroblasts",
    "alpha_alone",
    "alpha_in_UC",
    "alpha_interaction",
]
contrast_labels = ["Control FB − alone", "UC FB − alone", "UC FB − control FB", "α5β1 effect, alone", "α5β1 effect, UC FB", "α5β1 × FB context"]
states = [("CXCR4 neutrophil", "CXCR4"), ("MX1/ISG neutrophil", "MX1/ISG"), ("OSM neutrophil", "OSM"), ("PADI4 neutrophil", "PADI4")]
fig, axes = plt.subplots(1, 4, figsize=(174 / 25.4, 58 / 25.4), sharey=True)
fig.subplots_adjust(left=0.24, right=0.995, top=0.87, bottom=0.17, wspace=0.24)
for ax, (state_name, short_name) in zip(axes, states):
    subset = composition[composition.state.eq(state_name)]
    for y, contrast in enumerate(contrast_order):
        row = subset[subset.contrast.eq(contrast)].iloc[0]
        ax.errorbar(
            row.effect,
            y,
            xerr=[[row.effect - row.low], [row.high - row.effect]],
            fmt="o",
            ms=3.4,
            capsize=2,
            color=STATE[short_name],
            mfc=STATE[short_name] if row.q < 0.05 else "white",
            mew=0.9,
        )
    ax.axvline(0, color="#9CA3AF", ls=(0, (2, 2)), lw=0.65)
    ax.set_title(short_name, color=STATE[short_name], fontweight="bold", pad=4)
    ax.set_xlabel("CLR difference")
    ax.set_xlim(-3.5, 4.5)
    ax.grid(axis="x", color=GRID, lw=0.4)
axes[0].set_yticks(range(6), contrast_labels)
axes[0].set_ylim(5.5, -0.5)
native_ed6a = component(fig, "ED6_composition")

sensitivity = pd.read_csv(TABLES / "population_program_contrasts.csv")
sensitivity = sensitivity[
    sensitivity.program.eq("OSM_inflammation") & sensitivity.contrast.eq("alpha_interaction") & sensitivity.min_cells.eq(20)
]
labels = []
for row in sensitivity.itertuples():
    population = row.population.replace("All annotated", "Unresolved included").replace("Resolved; >=100 genes", "Resolved, ≥100 genes")
    version = "full" if row.version == "full_gene" else "identity genes removed"
    labels.append(f"{population} / {version}")
fig, ax = plt.subplots(figsize=(174 / 25.4, 54 / 25.4))
fig.subplots_adjust(left=0.36, right=0.64, top=0.84, bottom=0.18)
for y, row in enumerate(sensitivity.itertuples()):
    color = CONTROL if row.version == "full_gene" else "#9CA3AF"
    ax.errorbar(row.effect, y, xerr=[[row.effect - row.low], [row.high - row.effect]], fmt="o", color=color, ms=3.4, capsize=2)
    yy = 0.84 - (y + 0.5) * (0.66 / len(sensitivity))
    fig.text(0.75, yy, f"{row.q:.4f}", ha="center", va="center", fontsize=6.2)
    fig.text(0.92, yy, f"{row.q_signflip:.3f}", ha="center", va="center", fontsize=6.2)
ax.set_ylim(len(sensitivity) - 0.5, -0.5)
ax.set_yticks(range(len(labels)), labels)
ax.axvline(0, color="#9CA3AF", ls=(0, (2, 2)), lw=0.65)
ax.set_xlabel("OSM-program interaction (fixed-reference units)")
ax.set_xlim(-0.65, 0.16)
ax.grid(axis="x", color=GRID, lw=0.4)
fig.text(0.75, 0.91, "Model q", ha="center", fontsize=6.5, fontweight="bold")
fig.text(0.92, 0.91, "Exact q", ha="center", fontsize=6.5, fontweight="bold")
native_ed6b = component(fig, "ED6_sensitivity")

doc, page = new_page(137)
panel_header(page, 3, 2, "a", "Donor-aware neutrophil-state contrasts", "Six recorded donors; filled points denote BH q < 0.05; bars show 95% CI")
place(page, native_ed6a, None, (3, 12, 177, 69), False)
panel_header(page, 3, 73, "b", "OSM-program interaction sensitivity", "Four complete donors (≥20 cells per arm); exact sign-flip P = 0.125 across specifications")
place(page, native_ed6b, None, (3, 83, 177, 136), False)
save_figure(doc, 6, [str(TABLES / "composition_donor_contrasts.csv"), str(TABLES / "population_program_contrasts.csv")], 137, "Fully native forest plots; state colors match the main figures; inference notes incorporated into panel subtitles.")


# -----------------------------------------------------------------------------
# Extended Data Figure 7: ablation occupancy and recovered-cell fractions.
# -----------------------------------------------------------------------------
src18 = old(18)
doc, page = new_page(151)
panel_header(page, 3, 2, "a", "FAP ablation collapses neutrophil occupancy", "Shared inferred manifold; biological mice are the sample unit")
text(page, 35, 14, "DSS", 6.7, True, DSS)
text(page, 32, 18, "1,074 neutrophils · 6 mice", 6, False, MUTED)
text(page, 116, 14, "GCV + DSS", 6.7, True, BLOCKADE)
text(page, 114, 18, "2 neutrophils · 5 mice", 6, False, MUTED)
place(page, src18, (680, 158, 865, 450), (17, 20, 81, 54), True)
place(page, src18, (865, 158, 1055, 450), (99, 20, 163, 54), True)
panel_header(page, 3, 58, "b", "Relative recovery among all recovered cells", "Means ± SEM with individual biological mice; source FDR values retained")
clips = [(12, 504, 324.4, 785), (325.3, 504, 621.2, 785), (622.3, 504, 918.2, 785), (919.2, 504, 1216, 785)]
boxes = [(3, 68, 88, 106), (92, 68, 177, 106), (3, 109, 88, 147), (92, 109, 177, 147)]
for clip, box in zip(clips, boxes):
    place(page, src18, clip, box, True)
text(page, 3, 150, "Only two postablation neutrophils were recovered; no postablation trajectory inference is made.", 5.8, False, MUTED)
save_figure(doc, 7, ["Original Supplementary Figure S18 source vectors"], 151, "Condition occupancy split into balanced views with count callouts; four recovery panels enlarged and aligned; inference limitation retained.")


# -----------------------------------------------------------------------------
# Extended Data Figure 8: mouse trajectory occupancy and sample-level summaries.
# -----------------------------------------------------------------------------
src8 = current(8)
doc, page = new_page(139)
panel_header(page, 3, 2, "a", "Condition occupancy along a shared continuum", "Identical trajectory and UMAP axes across conditions")
text(page, 28, 14.5, "Control", 5.8, True, CONTROL)
text(page, 80, 14.5, "Severe DSS colitis", 5.8, True, DSS)
text(page, 134, 14.5, "DSS + α5β1 blockade", 5.8, True, BLOCKADE)
place(page, src8, mm_rect((3, 23, 177, 61)), (3, 17, 177, 56), True)
panel_header(page, 3, 60, "b", "Endpoint-state redistribution", "Fractions within total neutrophils")
place(page, src8, mm_rect((3, 91, 84, 136)), (3, 70, 87, 132), True)
panel_header(page, 92, 60, "c", "Neutrophil-program redistribution", "Each point is one sample unit")
place(page, src8, mm_rect((91, 91, 177, 136)), (92, 70, 177, 132), True)
text(page, 3, 137, "Treatment and RNA source batch are aligned; pseudotime is inferred, not measured time.", 5.8, False, MUTED)
save_figure(doc, 8, ["Prior canonical mouse-continuum and source-level summary vectors"], 139, "Removed duplicate headings, enlarged quantitative panels, aligned two-column evidence block, and tightened the canvas.")


# -----------------------------------------------------------------------------
# Extended Data Figure 9: source-only NET/FACS panels, untouched internally.
# -----------------------------------------------------------------------------
src21 = current(9)
doc, page = new_page(61)
short_titles = ["NETosis assay", "PADI4+ neutrophils", "OSM+ neutrophils", "MX1+ neutrophils", "CXCR4+ neutrophils"]
for index, x0 in enumerate([3, 38.2, 73.4, 108.6, 143.8]):
    text(page, x0, 5.5, "abcde"[index], 8, True)
    text(page, x0 + 6.2, 5.5, short_titles[index], 6.2, True)
# Crop away only the legacy all-caps source headings; statistical brackets begin below this boundary.
place(page, src21, mm_rect((3, 22, 177, 69)), (3, 10, 177, 60), True, "top")
save_figure(doc, 9, ["Original Supplementary Figure S21 NET/FACS raster graphs"], 61, "Five source graphs enlarged into equal columns; excess title canvas removed; graph interiors, points, error bars and significance brackets untouched.")


# -----------------------------------------------------------------------------
# Extended Data Figure 10: source-only in-vivo FACS panels, untouched internally.
# -----------------------------------------------------------------------------
src22 = old(22)
doc, page = new_page(170)
clips = [
    # Retain the complete statistical-annotation band while excluding the
    # superseded source titles and oversized source panel letters.
    (84, 148, 255, 360),
    (265, 152, 515, 360),
    (84, 452, 188, 540),
    (206, 452, 310, 540),
    (328, 452, 430, 540),
    (84, 560, 188, 699),
    (206, 560, 310, 699),
    (328, 560, 430, 699),
]
boxes = [
    (3, 9, 88, 65),
    (92, 9, 177, 65),
    (3, 73, 59, 116),
    (62, 73, 118, 116),
    (121, 73, 177, 116),
    (3, 123, 59, 169),
    (62, 123, 118, 169),
    (121, 123, 177, 169),
]
legacy_heading_masks = [
    [],
    [],
    [(3.0, 73.0, 10.0, 78.0)],
    [(62.0, 73.0, 69.5, 78.0)],
    [(121.0, 73.0, 128.5, 78.0)],
    [(13.2, 122.5, 14.0, 124.0)],
    [(71.8, 122.5, 73.6, 125.0)],
    [(130.5, 122.5, 132.0, 125.0)],
]
titles = ["α5β1+ fibroblasts", "α5β1+FAP+ fibroblasts", "OSM+ neutrophils", "PADI4+ neutrophils", "CXCR4+ neutrophils", "MPO MFI, OSM+", "MPO MFI, PADI4+", "MPO MFI, CXCR4+"]
for index, (clip, box) in enumerate(zip(clips, boxes)):
    place(page, src22, clip, box, True, "top")
    # Mask only the obsolete raster heading glyphs. These rectangles sit in
    # title whitespace and do not intersect axes, observations or statistics.
    for mask in legacy_heading_masks[index]:
        page.draw_rect(mm_rect(mask), color=(1, 1, 1), fill=(1, 1, 1), width=0.1, overlay=True)
    text(page, box[0], box[1] - 2.5, "abcdefgh"[index], 8, True)
    text(page, box[0] + 6.2, box[1] - 2.5, titles[index], 6.5, True)
save_figure(doc, 10, ["Original Supplementary Figure S22 in-vivo FACS raster graphs"], 170, "Reflowed untouched source graphs into a balanced 2+3+3 grid; bottom panels enlarged; non-data whitespace removed.")


# Combined review PDF and package metadata.
combined = fitz.open()
for number in range(1, 11):
    combined.insert_pdf(fitz.open(FIG_OUT / f"Extended_Data_Figure_{number:02d}.pdf"))
combined.set_metadata(
    {
        "title": "Extended Data Figures 1–10",
        "author": "UC fibroblast-neutrophil study",
        "subject": "Nature-style harmonized review set",
    }
)
combined.save(OUT / "Extended_Data_Figures_1_to_10.pdf", garbage=4, deflate=True)

for source in [
    HUM_TABLES / "09_fibroblast_biopsy_programs.tsv",
    HUM_TABLES / "22_neutrophil_state_program_summary.tsv",
    NEUT_SOURCE / "neutrophil_annotation_marker_dotplot_source.csv",
    NEUT_SOURCE / "neutrophil_cleanup_audit.csv",
    TABLES / "composition_donor_contrasts.csv",
    TABLES / "population_program_contrasts.csv",
    TABLES / "spatial_UC_pooled_effects.csv",
    TABLES / "spatial_restricted_label_tests.csv",
]:
    shutil.copy2(source, OUT / "source_data" / source.name)

(OUT / "reports/figure_inventory.json").write_text(json.dumps(inventory, indent=2), encoding="utf-8")
(OUT / "reports/build_manifest.json").write_text(
    json.dumps(
        {
            "figures": 10,
            "pages": 10,
            "max_width_mm": PAGE_W_MM,
            "max_height_mm": max(item["height_mm"] for item in inventory),
            "combined_sha256": hashlib.sha256((OUT / "Extended_Data_Figures_1_to_10.pdf").read_bytes()).hexdigest(),
            "facs_reanalysis": False,
            "new_experiments": False,
        },
        indent=2,
    ),
    encoding="utf-8",
)

print(json.dumps(inventory, indent=2))
