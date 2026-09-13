from pathlib import Path
import json
import re
import shutil
import tempfile
import zipfile

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image


Image.MAX_IMAGE_PIXELS = None

ROOT = Path(r"input_data/codex\giotto_codex_results")
CELL_FILE = ROOT / "additional_analyses" / "codex_extended_cells.csv"
RAW_QPTIFF = Path(r"input_data/codex_images\20240420_Pos1_IBD_Scan1-001.qptiff")
CONTROL_ARCHIVE = Path(r"input_data/codex_images\Control-20240731T174939Z-001.zip")
CONTROL_MEMBER = "Control/20240420_Pos2_control_Scan1.qptiff"
OUT = Path(r"input_data/mouse\CODEX_Figure_2")
OUT.mkdir(parents=True, exist_ok=True)

GROUPS = ["Control", "UC_Noninflamed", "UC_Inflamed"]
GROUP_LABEL = {"Control": "Control", "UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}
GROUP_COLOR = {"Control": "#4C78A8", "UC_Noninflamed": "#C9A227", "UC_Inflamed": "#D95F5F"}
REPRESENTATIVE = {"Control": "P06", "UC_Noninflamed": "P17", "UC_Inflamed": "P12"}

BROAD_ORDER = [
    "Epithelial Cell", "Enteroendocrine Cell", "Endothelial Cell", "Fibroblast",
    "Macrophage", "Dendritic Cell", "Neutrophil", "CD4 T", "CD8 T", "TReg",
    "B Cell", "Plasma B Cell"
]
BROAD_LABEL = {
    "Epithelial Cell": "Epithelial", "Enteroendocrine Cell": "Enteroendocrine",
    "Endothelial Cell": "Endothelial", "Fibroblast": "Fibroblast",
    "Macrophage": "Macrophage", "Dendritic Cell": "Dendritic",
    "Neutrophil": "Neutrophil", "CD4 T": "CD4 T", "CD8 T": "CD8 T",
    "TReg": "Treg", "B Cell": "B cell", "Plasma B Cell": "Plasma B"
}
BROAD_COLOR = dict(zip(BROAD_ORDER, [
    "#4E79A7", "#A0CBE8", "#59A14F", "#8CD17D", "#E15759", "#FF9D9A",
    "#F28E2B", "#B6992D", "#EDC948", "#B07AA1", "#76B7B2", "#9C755F"
]))
NEUT_ORDER = ["Neutrophil", "Neutrophil_PADI4", "Neutrophil_MX1", "Neutrophil_CXCR4", "Neutrophil_OSM"]
NEUT_LABEL = {
    "Neutrophil": "Unpolarized", "Neutrophil_PADI4": "PADI4+",
    "Neutrophil_MX1": "MX1+", "Neutrophil_CXCR4": "CXCR4+", "Neutrophil_OSM": "OSM+"
}
NEUT_COLOR = {
    "Neutrophil": "#D9D9D9", "Neutrophil_PADI4": "#A56CC1", "Neutrophil_MX1": "#28B6A8",
    "Neutrophil_CXCR4": "#56A0D3", "Neutrophil_OSM": "#FF5A5F"
}
CHANNEL_PAGES_QUARTER = {"DAPI": 111, "CD66b": 118, "FAP": 142, "a5B1": 151}
CHANNEL_COLOR = {
    "DAPI": np.array([0.13, 0.25, 1.00]), "CD66b": np.array([0.00, 0.95, 1.00]),
    "FAP": np.array([1.00, 0.00, 0.52]), "a5B1": np.array([0.15, 1.00, 0.20])
}

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 10,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"
})


def parse_scale_factor(description):
    match = re.search(r"<ScaleFactor>(.*?)</ScaleFactor>", description or "")
    if not match:
        raise ValueError("QPTIFF ScaleFactor was not found")
    return float(match.group(1))


def normalize_plane(array, low, high, gamma=.72):
    x = np.clip((array.astype(np.float32) - low) / max(high - low, 1e-6), 0, 1)
    return np.power(x, gamma)


def composite_planes(planes, scaling):
    rgb = np.zeros((*next(iter(planes.values())).shape, 3), dtype=np.float32)
    for channel, array in planes.items():
        signal = normalize_plane(array, *scaling[channel])
        rgb += signal[..., None] * CHANNEL_COLOR[channel]
    return 1 - np.exp(-1.10 * rgb)


def style_image_axis(ax, group, show_border=True):
    ax.set_xticks([]); ax.set_yticks([])
    if show_border:
        for spine in ax.spines.values():
            spine.set_color(GROUP_COLOR[group]); spine.set_linewidth(1.15)


def scale_bar(ax, xmin, xmax, ymin, ymax, color="white"):
    bar = 500
    bx = xmin + .055 * (xmax - xmin)
    by = ymax - .055 * (ymax - ymin)
    ax.plot([bx, bx + bar], [by, by], color=color, lw=3.2, solid_capstyle="butt")
    ax.text(bx + bar / 2, by - .023 * (ymax - ymin), "500 µm", color=color,
            ha="center", va="top", fontsize=7.2, fontweight="bold")


cells = pd.read_csv(CELL_FILE, usecols=["cell_ID", "PatientID", "Diagnosis2", "cell_type", "x", "y"], low_memory=False)
cells["broad_type"] = np.where(cells.cell_type.str.startswith("Neutrophil", na=False), "Neutrophil", cells.cell_type)

temp_root = Path(tempfile.mkdtemp(prefix="codex_control_qptiff_"))
try:
    with zipfile.ZipFile(CONTROL_ARCHIVE) as archive:
        archive.extract(CONTROL_MEMBER, temp_root)
    control_qptiff = temp_root / CONTROL_MEMBER
    uc_raw = Image.open(RAW_QPTIFF)
    control_raw = Image.open(control_qptiff)
    uc_raw.seek(0); full_scale = parse_scale_factor(str(uc_raw.tag_v2.get(270, "")))
    control_raw.seek(0); control_scale = parse_scale_factor(str(control_raw.tag_v2.get(270, "")))
    quarter_scale = full_scale / 4.0

    crop_boxes = {}; raw_crops = {}
    for group in GROUPS:
        z = cells[cells.PatientID.eq(REPRESENTATIVE[group])]
        pad = 90.0
        xmin = max(0, float(z.x.min() - pad)); xmax = float(z.x.max() + pad)
        ymin = max(0, float(z.y.min() - pad)); ymax = float(z.y.max() + pad)
        box = (int(xmin * quarter_scale), int(ymin * quarter_scale),
               int(xmax * quarter_scale), int(ymax * quarter_scale))
        crop_boxes[group] = (xmin, xmax, ymin, ymax, box)
        raw = control_raw if group == "Control" else uc_raw
        for channel, page in CHANNEL_PAGES_QUARTER.items():
            raw.seek(page)
            raw_crops[(group, channel)] = np.asarray(raw.crop(box), dtype=np.uint8)

    scaling = {}
    for channel in CHANNEL_PAGES_QUARTER:
        sampled = np.concatenate([raw_crops[(g, channel)][::8, ::8].ravel() for g in GROUPS])
        scaling[channel] = (float(np.percentile(sampled, 1.0)), float(np.percentile(sampled, 99.7)))
    uc_raw.close(); control_raw.close()
finally:
    shutil.rmtree(temp_root, ignore_errors=True)


fig = plt.figure(figsize=(18, 18.2), facecolor="white")
gs = GridSpec(3, 3, figure=fig, left=.05, right=.985, bottom=.115, top=.91,
              hspace=.14, wspace=.035)

for row, group in enumerate(GROUPS):
    patient = REPRESENTATIVE[group]
    z = cells[cells.PatientID.eq(patient)]
    xmin, xmax, ymin, ymax, _ = crop_boxes[group]
    planes = {ch: raw_crops[(group, ch)] for ch in CHANNEL_PAGES_QUARTER}
    rgb = composite_planes(planes, scaling)

    # Raw CODEX fluorescence.
    ax = fig.add_subplot(gs[row, 0])
    ax.imshow(rgb, extent=[xmin, xmax, ymax, ymin], interpolation="nearest", aspect="equal", rasterized=True)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymax, ymin); style_image_axis(ax, group)
    scale_bar(ax, xmin, xmax, ymin, ymax)
    ax.text(.025, .975, f"{GROUP_LABEL[group]}  |  {patient}", transform=ax.transAxes,
            ha="left", va="top", color="white", fontsize=9.2, fontweight="bold",
            bbox=dict(boxstyle="round,pad=.25", facecolor="#161616", edgecolor="none", alpha=.72))

    # All segmented cell annotations over a dimmed fluorescence reference.
    ax = fig.add_subplot(gs[row, 1])
    ax.imshow(rgb * .30, extent=[xmin, xmax, ymax, ymin], interpolation="nearest", aspect="equal", rasterized=True)
    for celltype in BROAD_ORDER:
        q = z[z.broad_type.eq(celltype)]
        ax.scatter(q.x, q.y, s=5.0, c=BROAD_COLOR[celltype], alpha=.77,
                   linewidths=0, rasterized=True)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymax, ymin); style_image_axis(ax, group)

    # Neutrophil subtype map at the identical registered coordinates.
    ax = fig.add_subplot(gs[row, 2])
    gray = np.clip(rgb.mean(axis=2), 0, 1)
    ax.imshow(gray, cmap="gray", vmin=0, vmax=1, extent=[xmin, xmax, ymax, ymin],
              interpolation="nearest", aspect="equal", alpha=.38, rasterized=True)
    other = z[~z.cell_type.isin(NEUT_ORDER)]
    if len(other) > 0:
        ax.scatter(other.x, other.y, s=1.2, c="#B6B6B6", alpha=.16, linewidths=0, rasterized=True)
    for subtype in NEUT_ORDER:
        q = z[z.cell_type.eq(subtype)]
        ax.scatter(q.x, q.y, s=13.0, c=NEUT_COLOR[subtype], alpha=.92,
                   edgecolors="#202020", linewidths=.20, rasterized=True)
    ax.set_xlim(xmin, xmax); ax.set_ylim(ymax, ymin); style_image_axis(ax, group)
    ax.text(.975, .965, f"n = {z.cell_type.isin(NEUT_ORDER).sum():,} neutrophils", transform=ax.transAxes,
            ha="right", va="top", fontsize=7.2, color="#202020",
            bbox=dict(boxstyle="round,pad=.22", facecolor="white", edgecolor="none", alpha=.82))

for col, title in enumerate([
    "CODEX fluorescence", "Complete cell annotations", "Neutrophil subtype annotations"
]):
    ax = fig.axes[col]
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax.text(-.075, 1.055, chr(ord('A') + col), transform=ax.transAxes, fontsize=17,
            fontweight="bold", va="top", ha="left", clip_on=False)

channel_handles = [mpl.lines.Line2D([], [], marker='o', ls='', color=CHANNEL_COLOR[ch], label=ch, markersize=5)
                   for ch in CHANNEL_PAGES_QUARTER]
broad_handles = [mpl.lines.Line2D([], [], marker='o', ls='', color=BROAD_COLOR[x], label=BROAD_LABEL[x], markersize=5)
                 for x in BROAD_ORDER]
neut_handles = [mpl.lines.Line2D([], [], marker='o', ls='', markerfacecolor=NEUT_COLOR[x],
                                markeredgecolor="#202020", markeredgewidth=.3,
                                label=NEUT_LABEL[x], markersize=5.5) for x in NEUT_ORDER]

fig.legend(handles=channel_handles, loc="lower left", bbox_to_anchor=(.055, .071),
           ncol=4, frameon=False, fontsize=7.5, title="Fluorescence channels", title_fontsize=8.2)
fig.legend(handles=broad_handles, loc="lower center", bbox_to_anchor=(.52, .035),
           ncol=6, frameon=False, fontsize=7.2, title="Complete cell annotations", title_fontsize=8.2,
           handletextpad=.25, columnspacing=.8)
fig.legend(handles=neut_handles, loc="lower right", bbox_to_anchor=(.98, .071),
           ncol=3, frameon=False, fontsize=7.5, title="Neutrophil subtypes", title_fontsize=8.2,
           handletextpad=.25, columnspacing=.8)

fig.suptitle("Figure 2—figure supplement 3. Registered CODEX microscopy and spatial cell annotations",
             x=.05, y=.968, ha="left", fontsize=16, fontweight="bold")
fig.text(.05, .944,
         "Matched fields show raw DAPI/CD66b/FAP/α5β1 fluorescence, all segmented cell types, and neutrophil phenotypic states.",
         ha="left", fontsize=9)

stem = "Figure_2_supplement_3_CODEX_microscopy_full_annotations"
for ext in ["png", "pdf", "svg"]:
    kwargs = {"dpi": 400} if ext == "png" else {}
    fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight", facecolor="white", **kwargs)
plt.close(fig)

manifest = {
    "source_cells": str(CELL_FILE), "uc_qptiff": str(RAW_QPTIFF),
    "control_qptiff_archive": str(CONTROL_ARCHIVE), "representative_patients": REPRESENTATIVE,
    "columns": ["CODEX fluorescence", "complete cell annotations", "neutrophil subtype annotations"],
    "broad_annotations": BROAD_ORDER, "neutrophil_subtypes": NEUT_ORDER,
    "outputs": [f"{stem}.png", f"{stem}.pdf", f"{stem}.svg"]
}
(OUT / "Figure_2_supplement_3_manifest.json").write_text(json.dumps(manifest, indent=2))
