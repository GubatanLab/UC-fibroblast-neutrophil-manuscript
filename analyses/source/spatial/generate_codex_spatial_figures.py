from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.style.use("dark_background")
BG = "#07090D"

OUT = Path("giotto_codex_results/nolan_neighborhoods")
FIG = Path("giotto_codex_results/figures")
FIG.mkdir(parents=True, exist_ok=True)
GROUPS = ["Control", "UC_Noninflamed", "UC_Inflamed"]
GROUP_LABELS = {"Control": "Control", "UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}

cells = pd.read_csv(OUT / "cells_with_nolan_neighborhoods.csv")
freq = pd.read_csv(OUT / "neighborhood_frequency_by_patient.csv")

# Select the patient nearest the multivariate group median neighborhood profile.
wide = freq.pivot(index=["PatientID", "Diagnosis2"], columns="neighborhood", values="frequency").fillna(0)
representatives = {}
for group in GROUPS:
    sub = wide.xs(group, level="Diagnosis2")
    center = sub.median(axis=0)
    representatives[group] = ((sub - center) ** 2).sum(axis=1).idxmin()
# Use the repaired P07 specimen as the designated inflamed representative.
representatives["UC_Inflamed"] = "P07"

neighborhoods = sorted(cells.neighborhood.unique(), key=lambda x: int(x.split(":")[0][2:]))
neighborhood_palette = [
    "#00C8FF", "#FF7A00", "#35E65C", "#FF3154", "#A66BFF",
    "#D98C5F", "#FF5EC4", "#B9C0CA", "#F1E400", "#00E0C6"
]
neigh_colors = dict(zip(neighborhoods, neighborhood_palette))
cell_types = sorted(cells.cell_type.unique())
cell_palette = [
    "#00C8FF", "#5B8CFF", "#FF7A00", "#35E65C", "#8EF07C", "#FF3154",
    "#A66BFF", "#C89BFF", "#FF5EC4", "#FF9ACD", "#B9C0CA", "#F1E400",
    "#E8E27A", "#00E0C6", "#73E6FF", "#D98C5F", "#E24AFF", "#79FFB2"
]
cell_colors = {v: cell_palette[i % len(cell_palette)] for i, v in enumerate(cell_types)}

def spatial_panel(ax, dat, color_col, color_map, size=0.32):
    ax.set_facecolor(BG)
    ax.scatter(dat.x, dat.y, c=dat[color_col].map(color_map), s=size, alpha=.96,
               linewidths=0, rasterized=True)
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

# Representative tissue maps colored by cellular neighborhood.
fig, axes = plt.subplots(1, 3, figsize=(18, 7))
fig.patch.set_facecolor(BG)
for ax, group in zip(axes, GROUPS):
    patient = representatives[group]
    dat = cells[cells.PatientID == patient]
    spatial_panel(ax, dat, "neighborhood", neigh_colors, size=1.05)
    ax.set_title(f"{GROUP_LABELS[group]} ({patient})", fontweight="bold")
handles = [Line2D([0], [0], marker="o", linestyle="", markersize=7,
                  markerfacecolor=neigh_colors[n], markeredgecolor="none", label=n.split(":")[0])
           for n in neighborhoods]
fig.legend(handles=handles, loc="lower center", ncol=10, frameon=False, title="Cellular neighborhood")
fig.suptitle("Representative CODEX tissue maps: cellular neighborhoods", fontweight="bold", fontsize=18)
fig.tight_layout(rect=[0, .09, 1, .94])
fig.savefig(FIG / "Figure_10_representative_CODEX_neighborhood_maps.png", dpi=300, bbox_inches="tight", facecolor=BG)
fig.savefig(FIG / "Figure_10_representative_CODEX_neighborhood_maps.pdf", bbox_inches="tight", facecolor=BG)
plt.close(fig)

# Same representative tissues colored by annotated cell type.
fig, axes = plt.subplots(1, 3, figsize=(18, 8))
fig.patch.set_facecolor(BG)
for ax, group in zip(axes, GROUPS):
    patient = representatives[group]
    dat = cells[cells.PatientID == patient]
    spatial_panel(ax, dat, "cell_type", cell_colors, size=1.00)
    ax.set_title(f"{GROUP_LABELS[group]} ({patient})", fontweight="bold")
handles = [Line2D([0], [0], marker="o", linestyle="", markersize=6,
                  markerfacecolor=cell_colors[n], markeredgecolor="none", label=n)
           for n in cell_types]
fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, title="Cell annotation", fontsize=8)
fig.suptitle("Representative CODEX tissue maps: cell annotations", fontweight="bold", fontsize=18)
fig.tight_layout(rect=[0, .19, 1, .94])
fig.savefig(FIG / "Figure_11_representative_CODEX_celltype_maps.png", dpi=300, bbox_inches="tight", facecolor=BG)
fig.savefig(FIG / "Figure_11_representative_CODEX_celltype_maps.pdf", bbox_inches="tight", facecolor=BG)
plt.close(fig)

# Cohort-wide neighborhood map, ordered by diagnostic group and patient.
patient_info = cells[["PatientID", "Diagnosis2"]].drop_duplicates()
patient_info["group_order"] = patient_info.Diagnosis2.map({g: i for i, g in enumerate(GROUPS)})
patient_info = patient_info.sort_values(["group_order", "PatientID"])
fig, axes = plt.subplots(6, 4, figsize=(16, 22))
fig.patch.set_facecolor(BG)
for ax, row in zip(axes.flat, patient_info.itertuples(index=False)):
    dat = cells[cells.PatientID == row.PatientID]
    spatial_panel(ax, dat, "neighborhood", neigh_colors, size=0.68)
    ax.set_title(f"{row.PatientID} | {GROUP_LABELS[row.Diagnosis2]}", fontsize=10, fontweight="bold")
for ax in axes.flat[len(patient_info):]:
    ax.axis("off")
fig.legend(handles=handles[:0], loc="lower center", frameon=False)
nh = [Line2D([0], [0], marker="o", linestyle="", markersize=7,
             markerfacecolor=neigh_colors[n], markeredgecolor="none", label=n.split(":")[0])
      for n in neighborhoods]
fig.legend(handles=nh, loc="lower center", ncol=10, frameon=False, title="Cellular neighborhood")
fig.suptitle("CODEX cellular-neighborhood maps across the cohort", fontweight="bold", fontsize=18)
fig.tight_layout(rect=[0, .035, 1, .975])
fig.savefig(FIG / "Figure_12_all_patient_CODEX_neighborhood_maps.png", dpi=300, bbox_inches="tight", facecolor=BG)
fig.savefig(FIG / "Figure_12_all_patient_CODEX_neighborhood_maps.pdf", bbox_inches="tight", facecolor=BG)
plt.close(fig)

pd.DataFrame([{"Diagnosis2": g, "representative_patient": p} for g, p in representatives.items()]).to_csv(
    OUT / "representative_CODEX_patients.csv", index=False
)
print(representatives)
