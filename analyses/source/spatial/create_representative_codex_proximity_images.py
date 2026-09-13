from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial import cKDTree

ROOT = Path("giotto_codex_results")
OUT = ROOT / "fap_high_a5b1_high_rerun"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GROUPS = ("UC_Noninflamed", "UC_Inflamed")
GROUP_LABEL = {"UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}
TARGETS = ("FAP-high", "a5B1-high", "FAP-high/a5B1-high")
SUBTYPES = (
    "Neutrophil", "Neutrophil_CXCR4", "Neutrophil_MX1",
    "Neutrophil_OSM", "Neutrophil_PADI4",
)
SUBTYPE_LABEL = {
    "Neutrophil": "Unassigned neutrophil",
    "Neutrophil_CXCR4": "CXCR4+ neutrophil",
    "Neutrophil_MX1": "MX1+ neutrophil",
    "Neutrophil_OSM": "OSM+ neutrophil",
    "Neutrophil_PADI4": "PADI4+ neutrophil",
}
SUBTYPE_COLOR = {
    "Neutrophil": "#F2F2F2",
    "Neutrophil_CXCR4": "#00D5FF",
    "Neutrophil_MX1": "#3478F6",
    "Neutrophil_OSM": "#FF9F1C",
    "Neutrophil_PADI4": "#7BE141",
}
TARGET_COLOR = "#FF3EA5"
CROP_HALF_WIDTH = 400.0
GRID_STEP = 100.0
SCALEBAR_UM = 100.0


def choose_representative_patients(subtype_results):
    wide = subtype_results.pivot_table(
        index=["PatientID", "Diagnosis2"],
        columns=["target", "neutrophil_subtype"],
        values="median_log2_observed_random_distance",
    )
    group_median = wide.groupby(level="Diagnosis2").transform("median")
    filled = wide.fillna(group_median)
    scores = ((filled - group_median) ** 2).mean(axis=1)
    return {
        group: scores.xs(group, level="Diagnosis2").idxmin()
        for group in GROUPS
    }, scores


def crop_mask(frame, cx, cy, half):
    return (
        frame.x.between(cx - half, cx + half)
        & frame.y.between(cy - half, cy + half)
    )


def select_representative_crop(frame, target, half=CROP_HALF_WIDTH):
    target_xy = frame.loc[frame[target], ["x", "y"]].to_numpy()
    neut = frame[frame.cell_type.isin(SUBTYPES)]
    neut_xy = neut[["x", "y"]].to_numpy()
    if not len(target_xy) or not len(neut_xy):
        raise ValueError(f"Missing target or neutrophil cells for {target}")

    global_median = float(np.median(cKDTree(target_xy).query(neut_xy)[0]))
    xmin, xmax = frame.x.min(), frame.x.max()
    ymin, ymax = frame.y.min(), frame.y.max()
    xs = np.arange(xmin + half, xmax - half + GRID_STEP, GRID_STEP)
    ys = np.arange(ymin + half, ymax - half + GRID_STEP, GRID_STEP)
    if not len(xs):
        xs = np.array([(xmin + xmax) / 2])
    if not len(ys):
        ys = np.array([(ymin + ymax) / 2])

    target_tree = cKDTree(target_xy)
    candidates = []
    for cx in xs:
        for cy in ys:
            mask = crop_mask(frame, cx, cy, half)
            local = frame[mask]
            local_neut = local[local.cell_type.isin(SUBTYPES)]
            n_target = int(local[target].sum())
            n_neutrophils = len(local_neut)
            n_cells = len(local)
            if n_cells < 250 or n_target < 8 or n_neutrophils < 12:
                continue
            local_median = float(np.median(target_tree.query(local_neut[["x", "y"]].to_numpy())[0]))
            candidates.append({
                "cx": float(cx), "cy": float(cy), "half_width_um": float(half),
                "n_cells": n_cells, "n_target_fibroblasts": n_target,
                "n_neutrophils": n_neutrophils,
                "highlighted_cells": n_target + n_neutrophils,
                "local_median_neutrophil_target_distance_um": local_median,
                "patient_median_neutrophil_target_distance_um": global_median,
                "representativeness": abs(np.log1p(local_median) - np.log1p(global_median)),
            })
    if not candidates:
        center = np.vstack([target_xy, neut_xy]).mean(axis=0)
        mask = crop_mask(frame, center[0], center[1], half)
        local = frame[mask]
        local_neut = local[local.cell_type.isin(SUBTYPES)]
        local_median = float(np.median(target_tree.query(local_neut[["x", "y"]].to_numpy())[0]))
        return {
            "cx": float(center[0]), "cy": float(center[1]), "half_width_um": float(half),
            "n_cells": len(local), "n_target_fibroblasts": int(local[target].sum()),
            "n_neutrophils": len(local_neut),
            "highlighted_cells": int(local[target].sum()) + len(local_neut),
            "local_median_neutrophil_target_distance_um": local_median,
            "patient_median_neutrophil_target_distance_um": global_median,
            "representativeness": abs(np.log1p(local_median) - np.log1p(global_median)),
        }

    c = pd.DataFrame(candidates)
    density_cut = c.highlighted_cells.quantile(.65)
    dense = c[c.highlighted_cells.ge(density_cut)]
    return dense.sort_values(["representativeness", "highlighted_cells"], ascending=[True, False]).iloc[0].to_dict()


usecols = [
    "cell_ID", "PatientID", "Diagnosis2", "cell_type", "x", "y",
    "FAPa_cell_clr", "a5B1_cell_clr",
]
cells = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv", usecols=usecols)
cells = cells[cells.Diagnosis2.isin(GROUPS)].copy()
thresholds = pd.read_csv(TAB / "00_patient_high_thresholds.csv")
cells = cells.merge(
    thresholds[["PatientID", "FAP_q75_clr", "a5B1_q75_clr"]],
    on="PatientID", how="left", validate="many_to_one",
)
cells["is_fibroblast"] = cells.cell_type.eq("Fibroblast")
cells["FAP-high"] = cells.is_fibroblast & cells.FAPa_cell_clr.ge(cells.FAP_q75_clr)
cells["a5B1-high"] = cells.is_fibroblast & cells.a5B1_cell_clr.ge(cells.a5B1_q75_clr)
cells["FAP-high/a5B1-high"] = cells["FAP-high"] & cells["a5B1-high"]

subtype_results = pd.read_csv(TAB / "17_neutrophil_subtype_proximity_by_patient.csv")
representatives, medoid_scores = choose_representative_patients(subtype_results)

crop_rows = []
for group in GROUPS:
    patient = representatives[group]
    frame = cells[cells.PatientID.eq(patient)].copy()
    for target in TARGETS:
        selected = select_representative_crop(frame, target)
        selected.update({
            "Diagnosis2": group,
            "group_label": GROUP_LABEL[group],
            "PatientID": patient,
            "target": target,
            "patient_medoid_score": float(medoid_scores.loc[(patient, group)]),
        })
        crop_rows.append(selected)

crops = pd.DataFrame(crop_rows)
crops.to_csv(TAB / "19_representative_codex_proximity_fields.csv", index=False)

plt.rcParams.update({"font.family": "Arial"})
fig, axes = plt.subplots(2, 3, figsize=(17, 11), facecolor="black")
for row_i, group in enumerate(GROUPS):
    patient = representatives[group]
    full = cells[cells.PatientID.eq(patient)]
    for col_i, target in enumerate(TARGETS):
        ax = axes[row_i, col_i]
        ax.set_facecolor("black")
        crop = crops[crops.Diagnosis2.eq(group) & crops.target.eq(target)].iloc[0]
        mask = crop_mask(full, crop.cx, crop.cy, crop.half_width_um)
        local = full[mask]

        ax.scatter(local.x, local.y, s=2.0, c="#6F7782", alpha=.20, linewidths=0, rasterized=True)
        target_cells = local[local[target]]
        ax.scatter(
            target_cells.x, target_cells.y, s=15, c=TARGET_COLOR,
            marker="D", alpha=.95, linewidths=.28, edgecolors="black", rasterized=True,
        )
        for subtype in SUBTYPES:
            q = local[local.cell_type.eq(subtype)]
            ax.scatter(
                q.x, q.y, s=11, c=SUBTYPE_COLOR[subtype], marker="o",
                alpha=.92, linewidths=.25, edgecolors="black", rasterized=True,
            )

        ax.set_xlim(crop.cx - crop.half_width_um, crop.cx + crop.half_width_um)
        ax.set_ylim(crop.cy + crop.half_width_um, crop.cy - crop.half_width_um)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#8A929D"); spine.set_linewidth(1.0)

        x0 = crop.cx - crop.half_width_um + 45
        y0 = crop.cy + crop.half_width_um - 45
        ax.plot([x0, x0 + SCALEBAR_UM], [y0, y0], color="white", linewidth=4, solid_capstyle="butt")
        ax.text(x0 + SCALEBAR_UM / 2, y0 - 22, f"{int(SCALEBAR_UM)} µm", color="white", ha="center", va="top", fontsize=9)
        ax.text(
            .98, .025,
            f"{int(crop.n_target_fibroblasts)} target fibroblasts · {int(crop.n_neutrophils)} neutrophils",
            transform=ax.transAxes, ha="right", va="bottom", color="white", fontsize=8,
            bbox={"facecolor": "black", "edgecolor": "none", "alpha": .65, "pad": 2},
        )

        if row_i == 0:
            ax.set_title(target, color="white", fontweight="bold", fontsize=15, pad=10)
        if col_i == 0:
            ax.text(
                -.055, .5, f"{GROUP_LABEL[group]}\n{patient}", transform=ax.transAxes,
                ha="right", va="center", rotation=90, color="white",
                fontweight="bold", fontsize=13,
            )

handles = [
    Line2D([0], [0], marker="D", color="none", label="High-state fibroblast",
           markerfacecolor=TARGET_COLOR, markeredgecolor="black", markersize=8),
]
handles.extend(
    Line2D([0], [0], marker="o", color="none", label=SUBTYPE_LABEL[s],
           markerfacecolor=SUBTYPE_COLOR[s], markeredgecolor="black", markersize=8)
    for s in SUBTYPES
)
legend = fig.legend(
    handles=handles, loc="upper center", bbox_to_anchor=(.5, .927), ncol=3,
    frameon=False, fontsize=9, labelcolor="white", handletextpad=.4, columnspacing=1.25,
)
fig.suptitle(
    "Representative CODEX fields: neutrophil subtype proximity to high-state fibroblasts",
    color="white", fontweight="bold", fontsize=20, y=.987,
)
fig.text(
    .5, .018,
    "Each dot is one segmented CODEX cell. Fields were selected for high target/neutrophil visibility while matching the patient-wide median proximity.",
    ha="center", color="#D0D5DC", fontsize=10,
)
fig.subplots_adjust(left=.07, right=.99, bottom=.06, top=.85, wspace=.055, hspace=.075)
png = FIG / "Figure_Representative_CODEX_Neutrophil_Subtype_Proximity.png"
pdf = FIG / "Figure_Representative_CODEX_Neutrophil_Subtype_Proximity.pdf"
tif = FIG / "Figure_Representative_CODEX_Neutrophil_Subtype_Proximity.tiff"
fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="black")
fig.savefig(pdf, bbox_inches="tight", facecolor="black")
fig.savefig(tif, dpi=300, bbox_inches="tight", facecolor="black", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

manifest = {
    "included_groups": list(GROUPS),
    "representative_patients": representatives,
    "representative_patient_method": "minimum mean squared deviation from group-median subtype proximity profile",
    "field_selection_method": "dense field with local median neutrophil-to-target distance closest to patient-wide median",
    "crop_width_um": 2 * CROP_HALF_WIDTH,
    "scale_bar_um": SCALEBAR_UM,
    "fibroblast_targets": list(TARGETS),
    "neutrophil_subtypes": list(SUBTYPES),
    "high_definition": "within-patient fibroblast CLR 75th percentile",
}
(OUT / "representative_codex_proximity_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print(f"Representative patients: {representatives}")
print(f"Saved {png}")
