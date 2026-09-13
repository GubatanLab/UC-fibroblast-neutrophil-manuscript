from pathlib import Path
import json
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image

Image.MAX_IMAGE_PIXELS = None

RAW_QPTIFF = Path(r"input_data/codex_images\20240420_Pos1_IBD_Scan1-001.qptiff")
ROOT = Path("giotto_codex_results")
OUT = ROOT / "fap_high_a5b1_high_rerun"
TAB = OUT / "tables"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

CROP_TABLE = TAB / "19_representative_codex_proximity_fields.csv"
GROUPS = ("UC_Noninflamed", "UC_Inflamed")
GROUP_LABEL = {"UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}
TARGETS = ("FAP-high", "a5B1-high", "FAP-high/a5B1-high")
CHANNEL_PAGES = {
    "DAPI": 0,
    "CD66b": 7,
    "CD15": 8,
    "OSM": 25,
    "CD11b": 28,
    "FAPa": 31,
    "CD16": 33,
    "CXCR4": 35,
    "PADI4": 38,
    "a5B1": 40,
    "MX1": 47,
}
CHANNEL_LABEL = {
    "DAPI": "DAPI", "CD66b": "CD66b", "CD15": "CD15", "OSM": "OSM",
    "CD11b": "CD11b", "FAPa": "FAPα", "CD16": "CD16", "CXCR4": "CXCR4",
    "PADI4": "PADI4", "a5B1": "a5B1", "MX1": "MX-1",
}
CONTEXT_COLOR = {
    "DAPI": np.array([0.12, 0.25, 1.00]),
    "FAPa": np.array([1.00, 0.00, 0.55]),
    "a5B1": np.array([0.15, 1.00, 0.20]),
    "Neutrophil-associated": np.array([0.00, 0.90, 1.00]),
}
SUBTYPE_COLOR = {
    "OSM": np.array([1.00, 0.50, 0.00]),
    "PADI4": np.array([0.45, 1.00, 0.10]),
    "MX1": np.array([0.72, 0.22, 1.00]),
    "CXCR4": np.array([0.00, 0.88, 1.00]),
}
NEUTROPHIL_CHANNELS = ("CD66b", "CD15", "CD11b", "CD16")
SUBTYPE_CHANNELS = ("OSM", "PADI4", "MX1", "CXCR4")
ATLAS_CHANNELS = NEUTROPHIL_CHANNELS + SUBTYPE_CHANNELS
DISPLAY_LOW_PERCENTILE = 1.0
DISPLAY_HIGH_PERCENTILE = 99.7
DISPLAY_GAMMA = 0.72
SCALEBAR_UM = 100.0


def parse_scale_factor(description):
    match = re.search(r"<ScaleFactor>(.*?)</ScaleFactor>", description or "")
    if not match:
        raise ValueError("QPTIFF ScaleFactor was not found")
    return float(match.group(1))


def pixel_box(row, scale):
    left = int(round((row.cx - row.half_width_um) * scale))
    top = int(round((row.cy - row.half_width_um) * scale))
    width = int(round(2 * row.half_width_um * scale))
    return left, top, left + width, top + width


def normalize(array, low, high, gamma=DISPLAY_GAMMA):
    x = np.clip((array.astype(np.float32) - low) / max(high - low, 1e-6), 0, 1)
    return np.power(x, gamma)


def saturating_composite(layers):
    rgb = np.zeros((*next(iter(layers.values()))[0].shape, 3), dtype=np.float32)
    for array, color, gain in layers.values():
        rgb += array[..., None] * color * gain
    return 1.0 - np.exp(-1.12 * rgb)


def add_scalebar(ax, shape, scale_factor):
    h, _ = shape[:2]
    length = SCALEBAR_UM * scale_factor
    x0, y0 = 55, h - 55
    ax.plot([x0, x0 + length], [y0, y0], color="white", linewidth=4, solid_capstyle="butt")
    ax.text(x0 + length / 2, y0 - 24, f"{int(SCALEBAR_UM)} µm",
            color="white", ha="center", va="bottom", fontsize=8)


crops = pd.read_csv(CROP_TABLE)
crops = crops[crops.Diagnosis2.isin(GROUPS) & crops.target.isin(TARGETS)].copy()
crops["panel_order"] = crops.Diagnosis2.map({g: i for i, g in enumerate(GROUPS)}) * len(TARGETS)
crops["panel_order"] += crops.target.map({t: i for i, t in enumerate(TARGETS)})
crops = crops.sort_values("panel_order").reset_index(drop=True)

raw = Image.open(RAW_QPTIFF)
raw.seek(CHANNEL_PAGES["DAPI"])
scale_factor = parse_scale_factor(str(raw.tag_v2.get(270, "")))
planes = {}
for channel, page in CHANNEL_PAGES.items():
    raw.seek(page)
    for row in crops.itertuples():
        planes[(row.Diagnosis2, row.target, channel)] = np.asarray(
            raw.crop(pixel_box(row, scale_factor)), dtype=np.uint8
        )
raw.close()

# Shared marker-specific display scaling across disease states and fields.
scaling = {}
normalized = {}
for channel in CHANNEL_PAGES:
    arrays = [array for key, array in planes.items() if key[2] == channel]
    sampled = np.concatenate([array[::8, ::8].ravel() for array in arrays])
    low = float(np.percentile(sampled, DISPLAY_LOW_PERCENTILE))
    high = float(np.percentile(sampled, DISPLAY_HIGH_PERCENTILE))
    scaling[channel] = (low, high)
    for key, array in planes.items():
        if key[2] == channel:
            normalized[key] = normalize(array, low, high)


def context_composite(group, target):
    neutrophil = np.maximum.reduce([
        normalized[(group, target, channel)] for channel in NEUTROPHIL_CHANNELS
    ])
    return saturating_composite({
        "DAPI": (normalized[(group, target, "DAPI")], CONTEXT_COLOR["DAPI"], .72),
        "FAPa": (normalized[(group, target, "FAPa")], CONTEXT_COLOR["FAPa"], 1.00),
        "a5B1": (normalized[(group, target, "a5B1")], CONTEXT_COLOR["a5B1"], .95),
        "Neutrophil-associated": (neutrophil, CONTEXT_COLOR["Neutrophil-associated"], .92),
    })


def subtype_composite(group, target):
    layers = {
        "DAPI": (normalized[(group, target, "DAPI")], np.array([.55, .62, .80]), .32),
    }
    for channel in SUBTYPE_CHANNELS:
        layers[channel] = (normalized[(group, target, channel)], SUBTYPE_COLOR[channel], .90)
    return saturating_composite(layers)


plt.rcParams.update({"font.family": "Arial"})
fig, axes = plt.subplots(4, 3, figsize=(17, 20), facecolor="black")
for group_i, group in enumerate(GROUPS):
    patient = crops[crops.Diagnosis2.eq(group)].PatientID.iloc[0]
    for target_i, target in enumerate(TARGETS):
        for layer_i, (layer_name, maker) in enumerate((
            ("Fibroblast + neutrophil markers", context_composite),
            ("Neutrophil subtype markers", subtype_composite),
        )):
            row_i = group_i * 2 + layer_i
            ax = axes[row_i, target_i]
            rgb = maker(group, target)
            ax.imshow(rgb, interpolation="nearest")
            ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor("black")
            for spine in ax.spines.values():
                spine.set_color("#8A929D"); spine.set_linewidth(1.0)
            add_scalebar(ax, rgb.shape, scale_factor)
            if row_i == 0:
                ax.set_title(target, color="white", fontweight="bold", fontsize=15, pad=10)
            if target_i == 0:
                ax.text(
                    -.055, .5, f"{GROUP_LABEL[group]} · {patient}\n{layer_name}",
                    transform=ax.transAxes, ha="right", va="center", rotation=90,
                    color="white", fontweight="bold", fontsize=11,
                )

context_handles = [
    Line2D([0], [0], marker="o", linestyle="none", label="DAPI",
           markerfacecolor=CONTEXT_COLOR["DAPI"], markeredgecolor="none", markersize=7),
    Line2D([0], [0], marker="o", linestyle="none", label="FAPα",
           markerfacecolor=CONTEXT_COLOR["FAPa"], markeredgecolor="none", markersize=7),
    Line2D([0], [0], marker="o", linestyle="none", label="a5B1",
           markerfacecolor=CONTEXT_COLOR["a5B1"], markeredgecolor="none", markersize=7),
    Line2D([0], [0], marker="o", linestyle="none",
           label="Neutrophil-associated: CD66b/CD15/CD11b/CD16",
           markerfacecolor=CONTEXT_COLOR["Neutrophil-associated"], markeredgecolor="none", markersize=7),
]
subtype_handles = [
    Line2D([0], [0], marker="o", linestyle="none", label=CHANNEL_LABEL[channel],
           markerfacecolor=SUBTYPE_COLOR[channel], markeredgecolor="none", markersize=7)
    for channel in SUBTYPE_CHANNELS
]
fig.legend(
    handles=context_handles + subtype_handles, loc="upper center", bbox_to_anchor=(.5, .948),
    ncol=4, frameon=False, fontsize=9, labelcolor="white", handletextpad=.4, columnspacing=1.1,
)
fig.suptitle(
    "Raw CODEX fluorescence: fibroblast proximity and neutrophil-state markers",
    color="white", fontweight="bold", fontsize=20, y=.988,
)
fig.text(
    .5, .018,
    "Full-resolution QPTIFF planes · shared per-marker scaling across UC states · neutrophil-associated signal is the pixelwise maximum of four marker planes",
    ha="center", color="#D0D5DC", fontsize=9.5,
)
fig.subplots_adjust(left=.075, right=.995, bottom=.045, top=.91, wspace=.045, hspace=.055)
montage_png = FIG / "Figure_Representative_CODEX_Raw_Fluorescence_All_Neutrophil_Markers.png"
montage_pdf = FIG / "Figure_Representative_CODEX_Raw_Fluorescence_All_Neutrophil_Markers.pdf"
montage_tif = FIG / "Figure_Representative_CODEX_Raw_Fluorescence_All_Neutrophil_Markers.tiff"
fig.savefig(montage_png, dpi=300, bbox_inches="tight", facecolor="black")
fig.savefig(montage_pdf, bbox_inches="tight", facecolor="black")
fig.savefig(montage_tif, dpi=300, bbox_inches="tight", facecolor="black", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

# Individually labeled raw marker atlas using the double-high representative field.
atlas_colors = {
    "CD66b": np.array([0.00, .90, 1.00]), "CD15": np.array([0.00, .90, 1.00]),
    "CD11b": np.array([0.00, .90, 1.00]), "CD16": np.array([0.00, .90, 1.00]),
    **SUBTYPE_COLOR,
}
fig, axes = plt.subplots(2, len(ATLAS_CHANNELS), figsize=(25, 7.2), facecolor="black")
for row_i, group in enumerate(GROUPS):
    patient = crops[crops.Diagnosis2.eq(group)].PatientID.iloc[0]
    target = "FAP-high/a5B1-high"
    dapi = normalized[(group, target, "DAPI")]
    for col_i, channel in enumerate(ATLAS_CHANNELS):
        ax = axes[row_i, col_i]
        rgb = saturating_composite({
            "DAPI": (dapi, np.array([.58, .62, .72]), .28),
            channel: (normalized[(group, target, channel)], atlas_colors[channel], 1.08),
        })
        ax.imshow(rgb, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_facecolor("black")
        for spine in ax.spines.values():
            spine.set_color("#8A929D"); spine.set_linewidth(.8)
        if row_i == 0:
            ax.set_title(CHANNEL_LABEL[channel], color="white", fontweight="bold", fontsize=12, pad=8)
        if col_i == 0:
            ax.text(
                -.08, .5, f"{GROUP_LABEL[group]}\n{patient}", transform=ax.transAxes,
                ha="right", va="center", rotation=90, color="white", fontweight="bold", fontsize=11,
            )
        if col_i == 0:
            add_scalebar(ax, rgb.shape, scale_factor)

fig.suptitle(
    "Raw neutrophil and subtype-marker planes in the double-high fibroblast field",
    color="white", fontweight="bold", fontsize=19, y=.985,
)
fig.text(
    .5, .02, "Each panel shows one raw marker plane over dim DAPI; intensity scaling is shared across UC states for that marker.",
    ha="center", color="#D0D5DC", fontsize=10,
)
fig.subplots_adjust(left=.045, right=.998, bottom=.08, top=.88, wspace=.025, hspace=.04)
atlas_png = FIG / "Figure_Representative_CODEX_Raw_Neutrophil_Marker_Atlas.png"
atlas_pdf = FIG / "Figure_Representative_CODEX_Raw_Neutrophil_Marker_Atlas.pdf"
atlas_tif = FIG / "Figure_Representative_CODEX_Raw_Neutrophil_Marker_Atlas.tiff"
fig.savefig(atlas_png, dpi=300, bbox_inches="tight", facecolor="black")
fig.savefig(atlas_pdf, bbox_inches="tight", facecolor="black")
fig.savefig(atlas_tif, dpi=300, bbox_inches="tight", facecolor="black", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

manifest = {
    "source_qptiff": str(RAW_QPTIFF),
    "source_kind": "raw full-resolution QPTIFF fluorescence planes",
    "included_groups": list(GROUPS),
    "channel_pages_zero_based": CHANNEL_PAGES,
    "neutrophil_associated_composite": list(NEUTROPHIL_CHANNELS),
    "subtype_marker_composite": list(SUBTYPE_CHANNELS),
    "shared_channel_display_scaling": {
        channel: {"low": scaling[channel][0], "high": scaling[channel][1]}
        for channel in scaling
    },
    "display_percentiles": [DISPLAY_LOW_PERCENTILE, DISPLAY_HIGH_PERCENTILE],
    "display_gamma": DISPLAY_GAMMA,
    "coordinate_scale_pixels_per_um": scale_factor,
    "crop_table": str(CROP_TABLE),
    "marker_atlas_field": "FAP-high/a5B1-high",
    "scale_bar_um": SCALEBAR_UM,
}
(OUT / "representative_codex_raw_neutrophil_marker_manifest.json").write_text(
    json.dumps(manifest, indent=2), encoding="utf-8"
)

print(f"Saved {montage_png}")
print(f"Saved {atlas_png}")
