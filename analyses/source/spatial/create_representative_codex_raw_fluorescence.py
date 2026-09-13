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
CHANNEL_PAGES = {"DAPI": 0, "CD66b": 7, "FAPa": 31, "a5B1": 40}
CHANNEL_LABEL = {"DAPI": "DAPI", "CD66b": "CD66b", "FAPa": "FAPα", "a5B1": "a5B1"}
CHANNEL_COLOR = {
    "DAPI": np.array([0.12, 0.25, 1.00]),
    "CD66b": np.array([0.00, 0.90, 1.00]),
    "FAPa": np.array([1.00, 0.00, 0.55]),
    "a5B1": np.array([0.15, 1.00, 0.20]),
}
CHANNEL_GAIN = {"DAPI": .82, "CD66b": 1.05, "FAPa": 1.05, "a5B1": 1.05}
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


def normalize_plane(array, low, high, gamma=DISPLAY_GAMMA):
    x = np.clip((array.astype(np.float32) - low) / max(high - low, 1e-6), 0, 1)
    return np.power(x, gamma)


def composite(planes, scaling):
    rgb = np.zeros((*next(iter(planes.values())).shape, 3), dtype=np.float32)
    for channel, array in planes.items():
        low, high = scaling[channel]
        signal = normalize_plane(array, low, high) * CHANNEL_GAIN[channel]
        rgb += signal[..., None] * CHANNEL_COLOR[channel]
    return 1.0 - np.exp(-1.15 * rgb)


crops = pd.read_csv(CROP_TABLE)
crops = crops[
    crops.Diagnosis2.isin(GROUPS) & crops.target.isin(TARGETS)
].copy()
crops["panel_order"] = crops.Diagnosis2.map({g: i for i, g in enumerate(GROUPS)}) * len(TARGETS)
crops["panel_order"] += crops.target.map({t: i for i, t in enumerate(TARGETS)})
crops = crops.sort_values("panel_order").reset_index(drop=True)

raw = Image.open(RAW_QPTIFF)
raw.seek(CHANNEL_PAGES["DAPI"])
scale_factor = parse_scale_factor(str(raw.tag_v2.get(270, "")))

planes = {}
boxes = {}
for channel, page in CHANNEL_PAGES.items():
    raw.seek(page)
    for row in crops.itertuples():
        key = (row.Diagnosis2, row.target, channel)
        box = pixel_box(row, scale_factor)
        boxes[(row.Diagnosis2, row.target)] = box
        planes[key] = np.asarray(raw.crop(box), dtype=np.uint8)
raw.close()

# Apply identical display scaling to both UC states and all fields for each marker.
scaling = {}
for channel in CHANNEL_PAGES:
    sampled = np.concatenate([
        array[::8, ::8].ravel()
        for key, array in planes.items() if key[2] == channel
    ])
    low = float(np.percentile(sampled, DISPLAY_LOW_PERCENTILE))
    high = float(np.percentile(sampled, DISPLAY_HIGH_PERCENTILE))
    scaling[channel] = (low, high)

plt.rcParams.update({"font.family": "Arial"})
fig, axes = plt.subplots(2, 3, figsize=(17, 11), facecolor="black")
for row_i, group in enumerate(GROUPS):
    for col_i, target in enumerate(TARGETS):
        ax = axes[row_i, col_i]
        panel_planes = {
            channel: planes[(group, target, channel)]
            for channel in CHANNEL_PAGES
        }
        rgb = composite(panel_planes, scaling)
        ax.imshow(rgb, interpolation="nearest")
        ax.set_facecolor("black")
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color("#8A929D"); spine.set_linewidth(1.0)

        if row_i == 0:
            ax.set_title(target, color="white", fontweight="bold", fontsize=15, pad=10)
        if col_i == 0:
            patient = crops[
                crops.Diagnosis2.eq(group) & crops.target.eq(target)
            ].PatientID.iloc[0]
            ax.text(
                -.055, .5, f"{GROUP_LABEL[group]}\n{patient}", transform=ax.transAxes,
                ha="right", va="center", rotation=90, color="white",
                fontweight="bold", fontsize=13,
            )

        scale_pixels = SCALEBAR_UM * scale_factor
        h, w = rgb.shape[:2]
        x0, y0 = 55, h - 55
        ax.plot([x0, x0 + scale_pixels], [y0, y0], color="white", linewidth=4, solid_capstyle="butt")
        ax.text(x0 + scale_pixels / 2, y0 - 24, f"{int(SCALEBAR_UM)} µm",
                color="white", ha="center", va="bottom", fontsize=9)

handles = [
    Line2D([0], [0], marker="o", linestyle="none", label=CHANNEL_LABEL[channel],
           markerfacecolor=CHANNEL_COLOR[channel], markeredgecolor="none", markersize=8)
    for channel in CHANNEL_PAGES
]
fig.legend(
    handles=handles, loc="upper center", bbox_to_anchor=(.5, .927), ncol=4,
    frameon=False, fontsize=10, labelcolor="white", handletextpad=.45, columnspacing=1.4,
)
fig.suptitle(
    "Representative raw CODEX fluorescence: neutrophil proximity to high-state fibroblasts",
    color="white", fontweight="bold", fontsize=20, y=.987,
)
fig.text(
    .5, .018,
    "Full-resolution QPTIFF planes · identical per-marker intensity scaling across UC states · fields match the proximity analysis",
    ha="center", color="#D0D5DC", fontsize=10,
)
fig.subplots_adjust(left=.07, right=.99, bottom=.06, top=.85, wspace=.055, hspace=.075)

png = FIG / "Figure_Representative_CODEX_Raw_Fluorescence_Proximity.png"
pdf = FIG / "Figure_Representative_CODEX_Raw_Fluorescence_Proximity.pdf"
tif = FIG / "Figure_Representative_CODEX_Raw_Fluorescence_Proximity.tiff"
fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="black")
fig.savefig(pdf, bbox_inches="tight", facecolor="black")
fig.savefig(tif, dpi=300, bbox_inches="tight", facecolor="black", pil_kwargs={"compression": "tiff_lzw"})
plt.close(fig)

manifest = {
    "source_qptiff": str(RAW_QPTIFF),
    "source_kind": "raw full-resolution QPTIFF fluorescence planes",
    "included_groups": list(GROUPS),
    "representative_patients": {
        group: crops[crops.Diagnosis2.eq(group)].PatientID.iloc[0]
        for group in GROUPS
    },
    "channel_pages_zero_based": CHANNEL_PAGES,
    "channel_display_colors_rgb": {k: v.tolist() for k, v in CHANNEL_COLOR.items()},
    "shared_channel_display_scaling": {
        channel: {"low": scaling[channel][0], "high": scaling[channel][1]}
        for channel in scaling
    },
    "display_percentiles": [DISPLAY_LOW_PERCENTILE, DISPLAY_HIGH_PERCENTILE],
    "display_gamma": DISPLAY_GAMMA,
    "coordinate_scale_pixels_per_um": scale_factor,
    "crop_table": str(CROP_TABLE),
    "pixel_boxes_left_top_right_bottom": {
        f"{group}|{target}": list(boxes[(group, target)])
        for group in GROUPS for target in TARGETS
    },
    "scale_bar_um": SCALEBAR_UM,
}
(OUT / "representative_codex_raw_fluorescence_manifest.json").write_text(
    json.dumps(manifest, indent=2), encoding="utf-8"
)

print(f"Scale factor: {scale_factor:.4f} pixels/um")
print(f"Shared scaling: {scaling}")
print(f"Saved {png}")
