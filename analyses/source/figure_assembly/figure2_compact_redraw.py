"""Native Figure 2 redraw: original CODEX fields and unchanged patient-level data."""
from __future__ import annotations
import json
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
WORK = HERE.parent
DATA = HERE / "source_data" / "figure2"
WIDTH, HEIGHT = 183., 168.
INK, MUTED, GRID = "#25282C", "#5F666D", "#E4E7EB"
GROUPS = ["Control", "UC_Noninflamed", "UC_Inflamed"]
COLORS = dict(zip(GROUPS, ["#3B6FB6", "#C6A43A", "#C84C4C"]))
MARKERS = dict(zip(GROUPS, ["o", "^", "s"]))
GROUP_LABELS = ["Control", "UC\nnoninfl.", "UC\ninflamed"]
STATES = ["FAP+/a5B1-", "FAP-/a5B1+", "FAP+/a5B1+", "FAP-/a5B1-"]
STATE_LABELS = ["FAP+ / α5β1-", "FAP- / α5β1+", "FAP+ / α5β1+", "FAP- / α5β1-"]
CELL_COLORS = {
    "Epithelial": "#4E79A7", "Enteroendocrine": "#A0CBE8",
    "Endothelial": "#59A14F", "Fibroblast": "#8CD17D",
    "Macrophage": "#E15759", "Dendritic": "#FF9D9A",
    "Neutrophil": "#F28E2B", "CD4 T": "#B6992D",
    "CD8 T": "#EDC948", "Treg": "#B07AA1",
    "B cell": "#76B7B2", "Plasma B": "#9C755F",
}
FIELD_SPECS = [
    ((2471, 611, 4670, 2078), "Control | P06", "Control", 2624.1),
    ((2607, 2471, 4534, 4241), "UC noninflamed | P17", "UC_Noninflamed", 3338.5),
    ((2471, 4734, 4670, 6001), "UC inflamed | P12", "UC_Inflamed", 3144.9),
]


def text(fig, x, y, value, size=5.3, weight="normal", color=INK, **kwargs):
    return fig.text(x / WIDTH, 1 - y / HEIGHT, value, fontsize=size, fontweight=weight,
                    color=color, va=kwargs.pop("va", "top"), **kwargs)


def axis(fig, x, y, width, height, name):
    return fig.add_axes([x / WIDTH, 1 - (y + height) / HEIGHT, width / WIDTH, height / HEIGHT], label=name)


def heading(fig, x, y, letter, caption):
    text(fig, x, y - .2, letter, size=7, weight="bold")
    text(fig, x + 4.5, y, caption, size=6.5, weight="bold")


def style_axis(ax, grid="y"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#60666D")
    ax.tick_params(length=2, width=.5, pad=1.4, labelsize=5.1, colors=INK)
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color=GRID, linewidth=.35)


def read_data():
    names = {"cohort": "01_patient_cell_counts.csv", "abundance": "01_fibroblast_state_abundance_by_patient.csv",
             "spatial": "03_spatial_enrichment_by_patient.csv", "markers": "08_neutrophil_marker_proximity_effects_by_patient.csv",
             "multiscale": "analysis_3_multiscale_patient_effects.csv", "distance": "analysis_4_patient_distance_bin_marker_gradients.csv"}
    return {key: pd.read_csv(DATA / name) for key, name in names.items()}


def letterbox(field):
    # Preserve the existing source-field pixels, orientation and 3:2 letterboxing.
    width, height = field.size
    if np.isclose(width / height, 1.5, rtol=0, atol=.001):
        return field
    size = (int(np.ceil(height * 1.5)), height) if width / height < 1.5 else (width, int(np.ceil(width / 1.5)))
    canvas = Image.new("RGB", size, (0, 0, 0))
    canvas.paste(field, ((size[0] - width) // 2, (size[1] - height) // 2))
    return canvas


def draw_images(fig):
    heading(fig, 3, 2.7, "a", "Representative CODEX cell-annotation maps")
    source_path = WORK / "Figure 2/Figure_2_supplement_3_CODEX_microscopy_full_annotations.png"
    source = Image.open(source_path).convert("RGB")
    records = []
    for index, (box, label, group, field_width_um) in enumerate(FIELD_SPECS):
        field = source.crop(box)
        field = field.crop((6, 6, field.width - 6, field.height - 6))
        active_width, active_height = field.size
        field = letterbox(field)
        x = 7.5 + index * 57.5
        text(fig, x + 27.75, 7.2, label, size=5.6, weight="bold", ha="center")
        ax = axis(fig, x, 11, 55.5, 37, f"a_{group}")
        ax.imshow(field, interpolation="lanczos", resample=True)
        ax.set_axis_off()
        ax.set_box_aspect(2 / 3)
        ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, fill=False,
                               edgecolor=COLORS[group], linewidth=.65, clip_on=False))
        pad_left = (field.width - active_width) / 2
        pad_bottom = (field.height - active_height) / 2
        bar_x0 = (pad_left + active_width * .055) / field.width
        bar_x1 = bar_x0 + 500 / field_width_um * active_width / field.width
        bar_y = (pad_bottom + active_height * .055) / field.height
        ax.plot([bar_x0, bar_x1], [bar_y, bar_y], transform=ax.transAxes,
                color="black", linewidth=2.5, solid_capstyle="butt", zorder=6)
        ax.plot([bar_x0, bar_x1], [bar_y, bar_y], transform=ax.transAxes,
                color="white", linewidth=1.4, solid_capstyle="butt", zorder=7)
        ax.text((bar_x0 + bar_x1) / 2, bar_y + .025, "500 µm", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=5, fontweight="bold", color="white",
                bbox=dict(facecolor="black", edgecolor="none", alpha=.65, pad=.12), zorder=7)
        records.append({"group": group, "label": label, "source_crop_px": list(box),
                        "source_field_width_um": field_width_um, "scale_bar_um": 500,
                        "image_pixels": list(field.size), "display_mm": [55.5, 37]})
    source.close()
    for index, (label, color) in enumerate(CELL_COLORS.items()):
        col, row = divmod(index, 2)
        x, y = 8.5 + col * 28.5, 50.8 + row * 3
        fig.add_artist(Line2D([x / WIDTH], [1 - y / HEIGHT], linestyle="", marker="o", markersize=2.8,
                              color=color, transform=fig.transFigure))
        text(fig, x + 1.8, y, label, size=5, va="center")
    return records


def draw_abundance(fig, data, x, letter, title, column, ylabel, name, upper):
    heading(fig, x, 58.5, letter, title)
    ax = axis(fig, x + 11, 66, 43, 29, name)
    values = [pd.to_numeric(data.loc[data.Diagnosis2.eq(g), column], errors="coerce").dropna().to_numpy() * 100 for g in GROUPS]
    boxes = ax.boxplot(values, positions=np.arange(3), widths=.52, patch_artist=True, showfliers=False,
                       medianprops={"color": INK, "linewidth": .85},
                       whiskerprops={"color": "#68707A", "linewidth": .6},
                       capprops={"color": "#68707A", "linewidth": .6})
    rng = np.random.default_rng(20260901)
    for index, (group, vals, box) in enumerate(zip(GROUPS, values, boxes["boxes"])):
        box.set(facecolor=mpl.colors.to_rgba(COLORS[group], .16), edgecolor=COLORS[group], linewidth=.65)
        jitter = rng.uniform(-.12, .12, len(vals))
        ax.scatter(index + jitter, vals, s=10, color=COLORS[group], marker=MARKERS[group],
                   edgecolors="white", linewidths=.3, zorder=3)
        text(fig, x + 11 + (index + .5) * 43 / 3, 62.6, f"n = {len(vals)}", size=5.1, ha="center", color=MUTED)
    ax.set_xlim(-.5, 2.5)
    ax.set_ylim(-.8, upper)
    ax.set_xticks(range(3), GROUP_LABELS)
    ax.set_yticks(np.arange(0, upper, 5))
    style_axis(ax)
    text(fig, x + 2, 80.5, ylabel, size=5.2, rotation=90, va="center", ha="center")
    return {g: len(v) for g, v in zip(GROUPS, values)}


def matrix_plot(fig, values, labels, x, y, width, height, cb_x, name):
    ax = axis(fig, x, y, width, height, name)
    matrix = np.asarray(values, dtype=float)
    limit = np.nanmax(np.abs(matrix))
    im = ax.pcolormesh(matrix, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
                       edgecolors="white", linewidth=.4)
    ax.set_xlim(0, matrix.shape[1])
    ax.set_ylim(matrix.shape[0], 0)
    ax.set_axis_off()
    for row in range(matrix.shape[0]):
        text(fig, x - 1.6, y + (row + .5) * height / matrix.shape[0], labels[row], size=5.1, ha="right", va="center")
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            ax.text(col + .5, row + .5, f"{value:.2f}", fontsize=5.1, ha="center", va="center",
                    color="white" if abs(value) > .55 * limit else INK)
    cax = axis(fig, cb_x, y, 1.3, height, f"{name}_colorbar")
    ticks = [-1, 0, 1] if limit > 1 else [-.4, 0, .4]
    cb = fig.colorbar(im, cax=cax, ticks=ticks)
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=1.3, width=.4, labelsize=5, pad=1)
    return ax


def draw_spatial(fig, spatial):
    heading(fig, 123, 58.5, "d", "Spatial enrichment, 50 µm")
    text(fig, 127.5, 62.6, "Patient-median log2 observed / expected", size=5, color=MUTED)
    mat = spatial.loc[spatial.radius_um.eq(50)].groupby(["fibroblast_state", "Diagnosis2"]).log2_observed_expected.median().unstack().reindex(index=STATES, columns=GROUPS)
    matrix_plot(fig, mat, STATE_LABELS, 146, 69.8, 28.5, 25.2, 177, "d_spatial")
    for index, label in enumerate(GROUP_LABELS):
        text(fig, 146 + (index + .5) * 9.5, 96.1, label, size=5, ha="center")
    return mat.to_dict()


def draw_cohort_key(fig):
    text(fig, 7.5, 105.2, "Cohort colors / symbols (b, c, e, g)", size=5.1, weight="bold", va="center")
    for x, group, label in [(64, GROUPS[0], "Control"), (97, GROUPS[1], "UC noninflamed"), (145, GROUPS[2], "UC inflamed")]:
        fig.add_artist(Line2D([x / WIDTH], [1 - 105.2 / HEIGHT], marker=MARKERS[group],
                              color=COLORS[group], markersize=3, linestyle="", transform=fig.transFigure))
        text(fig, x + 2.1, 105.2, label, size=5.1, va="center")


def draw_multiscale(fig, multiscale):
    heading(fig, 3, 112, "e", "Multiscale niche topology")
    text(fig, 7.5, 116, "Patient-median log2 observed / expected", size=5, color=MUTED)
    for x, ls, caption in [(7.5, "-", "FAP+ / α5β1-"), (34, "--", "FAP+ / α5β1+")]:
        fig.add_artist(Line2D([x / WIDTH, (x + 4) / WIDTH], [1 - 120 / HEIGHT] * 2,
                              color=INK, linestyle=ls, linewidth=.8, transform=fig.transFigure))
        text(fig, x + 5, 120, caption, size=5, va="center")
    ax = axis(fig, 14, 124.5, 43, 30, "e_multiscale")
    records = []
    for state, ls in [(STATES[0], "-"), (STATES[2], "--")]:
        for group in GROUPS:
            sub = multiscale.loc[multiscale.fibroblast_state.eq(state) & multiscale.Diagnosis2.eq(group)]
            med = sub.groupby("radius_um").log2_observed_expected.median().sort_index()
            ax.plot(med.index, med.values, color=COLORS[group], linestyle=ls, linewidth=1,
                    label=f"{group}|{state}")
            records.append({"group": group, "state": state, "radii_um": med.index.tolist(), "medians": med.tolist()})
    ax.axhline(0, color="#808790", linewidth=.6)
    ax.set_xlim(12, 153)
    ax.set_ylim(-2.65, .3)
    ax.set_xticks([25, 50, 100, 150])
    ax.set_yticks([-2.5, -2, -1.5, -1, -.5, 0])
    style_axis(ax)
    text(fig, 35.5, 158.3, "Radius (µm)", size=5.3, ha="center")
    return records


def draw_distance(fig, distance):
    heading(fig, 63, 112, "f", "Distance-linked phenotype")
    text(fig, 67.5, 116, "Inflamed UC; FAP+ / α5β1- niche", size=5, color=MUTED)
    text(fig, 78, 121.2, "Patient-median robust z", size=5, color=MUTED)
    sub = distance.loc[distance.Diagnosis2.eq("UC_Inflamed") & distance.fibroblast_state.eq(STATES[0])]
    bins = ["0–25", "25–50", "50–100", ">100"]
    markers = ["CXCR4", "OSM", "CD16", "CD11b"]
    mat = sub.groupby(["marker", "distance_bin"]).median_marker_robust_z.median().unstack().reindex(index=markers, columns=bins)
    matrix_plot(fig, mat, markers, 78, 125.5, 35.2, 25.2, 116, "f_distance")
    for index, label in enumerate(["0-25", "25-50", "50-100", ">100"]):
        text(fig, 78 + (index + .5) * 8.8, 152.1, label, size=5, ha="center")
    text(fig, 95.6, 158.3, "Distance (µm)", size=5.3, ha="center")
    return mat.to_dict()


def coupling_data(data):
    spatial, markers = data["spatial"], data["markers"]
    x = spatial.loc[spatial.fibroblast_state.eq(STATES[0]) & spatial.radius_um.eq(50), ["PatientID", "Diagnosis2", "log2_observed_expected"]]
    y = markers.loc[markers.fibroblast_state.eq(STATES[0]) & markers.feature.eq("NF-niche program") & markers.method.eq("Proximal-vs-distant median difference"), ["PatientID", "effect"]]
    return x.merge(y, on="PatientID", how="inner").dropna()


def draw_coupling(fig, data):
    heading(fig, 123, 112, "g", "Spatial-phenotypic coupling")
    joined = coupling_data(data)
    rho, pval = spearmanr(joined.log2_observed_expected, joined.effect)
    text(fig, 127.5, 116, f"ρ = {rho:.2f}; P = {pval:.3g}; n = {len(joined)}", size=5.1, color=MUTED)
    ax = axis(fig, 136, 124.5, 43, 30, "g_coupling")
    for group in GROUPS:
        sub = joined.loc[joined.Diagnosis2.eq(group)]
        ax.scatter(sub.log2_observed_expected, sub.effect, s=14, color=COLORS[group], marker=MARKERS[group],
                   edgecolors="white", linewidths=.35, zorder=3)
    coefficients = np.polyfit(joined.log2_observed_expected, joined.effect, 1)
    xx = np.linspace(joined.log2_observed_expected.min(), joined.log2_observed_expected.max(), 50)
    ax.plot(xx, np.polyval(coefficients, xx), color=INK, linewidth=.8)
    ax.set_xlim(-1.55, 1.4)
    ax.set_ylim(-.42, 2.1)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([0, .5, 1, 1.5, 2])
    style_axis(ax)
    text(fig, 127.5, 139.5, "Neutrophil niche-program effect", size=5.1, rotation=90, va="center", ha="center")
    text(fig, 157.5, 158.3, "50 µm spatial enrichment (log2 O/E)", size=5.1, ha="center")
    return {"n_regions": len(joined), "spearman_rho": float(rho), "p_value": float(pval),
            "patient_ids": joined.PatientID.tolist(), "linear_visual_summary_coefficients": coefficients.tolist()}


def build_figure2(save_callback):
    with mpl.rc_context({"font.family": "Arial", "font.size": 5.3, "axes.linewidth": .55,
                         "pdf.fonttype": 42, "ps.fonttype": 42}):
        data = read_data()
        fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor="white")
        audit = {"page_mm": [WIDTH, HEIGHT], "canonical_outputs_modified": False,
                 "analysis_formulas_changed": False, "abundance_display_units": "Percent (source fractions times 100)",
                 "source_image_pixels_changed": False, "cell_annotation_colors": CELL_COLORS}
        audit["representative_fields"] = draw_images(fig)
        audit["neutrophil_region_counts"] = draw_abundance(fig, data["cohort"], 3, "b", "Neutrophil abundance",
                                                          "fraction_all_cells_neutrophils", "Neutrophils (% all cells)", "b_abundance", 21.8)
        sub = data["abundance"].loc[data["abundance"].fibroblast_state.eq(STATES[2])]
        audit["fibroblast_region_counts"] = draw_abundance(fig, sub, 63, "c", "FAP+ / α5β1+ fibroblasts",
                                                          "fraction_fibroblasts", "State (% fibroblasts)", "c_abundance", 27.8)
        audit["spatial_medians"] = draw_spatial(fig, data["spatial"])
        draw_cohort_key(fig)
        audit["multiscale_curves"] = draw_multiscale(fig, data["multiscale"])
        audit["distance_medians"] = draw_distance(fig, data["distance"])
        audit["coupling"] = draw_coupling(fig, data)
        fig.canvas.draw()
        (HERE / "Figure_2_REDRAW_NONCANONICAL_manifest.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
        save_callback(fig, 2)
