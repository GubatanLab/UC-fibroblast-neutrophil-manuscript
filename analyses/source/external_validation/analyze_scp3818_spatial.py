import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


ROOT = pathlib.Path(r"input_data/mouse\Figure 6 Spatial External Validation")
DATA = ROOT / "data" / "SCP3818"
OUT = ROOT / "results" / "SCP3818"
OUT.mkdir(parents=True, exist_ok=True)

SECTIONS = {
    "UC1 inflamed": "UC1_inflamed",
    "UC1 less inflamed": "UC1_less_inflamed",
}
FIBRO = ["FAP", "PDPN", "THY1", "ITGA5"]
RECRUIT = ["CXCL2", "CXCL5", "CXCL6", "CSF3", "IL6", "ICAM1"]
NEUT = ["CSF3R", "CXCR1", "CXCR2", "FCGR3B", "CEACAM8", "MPO"]
FEEDBACK = ["OSM", "IL1B", "TNF"]


def load_section(section, stem):
    cluster = json.loads((DATA / f"{stem}_cluster.json").read_text(encoding="utf-8"))
    df = pd.DataFrame({
        "cell_id": pd.Series(cluster["data"]["cells"], dtype=str),
        "x": cluster["data"]["x"],
        "y": cluster["data"]["y"],
        "cell_type": cluster["data"]["annotations"],
    }).set_index("cell_id")
    for gene in FIBRO + RECRUIT + NEUT + FEEDBACK:
        path = DATA / f"{stem}_{gene}.tsv"
        if path.exists():
            g = pd.read_csv(path, sep="\t", dtype={"cell_id": str}).set_index("cell_id")
            df[gene] = g["expression"].reindex(df.index).fillna(0.0)
    for genes, name in [(RECRUIT, "recruitment_score"), (NEUT, "neutrophil_score"), (FEEDBACK, "feedback_score")]:
        present = [g for g in genes if g in df]
        ranks = pd.concat([df[g].rank(pct=True) for g in present], axis=1)
        df[name] = ranks.mean(axis=1)
    df["activated_fap"] = (
        (df.cell_type == "Stromal") & (df.FAP > 0) & ((df.PDPN > 0) | (df.THY1 > 0))
    )
    df["alpha5_fap"] = df["activated_fap"] & (df.ITGA5 > 0)
    df["neutrophil_like"] = (
        (df.cell_type == "Myeloid") & (df[[g for g in NEUT if g in df]] > 0).sum(axis=1).ge(2)
    )
    df["section"] = section
    return df


def spatial_test(df, rng, n_perm=1000):
    stromal = df[df.cell_type.eq("Stromal")]
    active = stromal[stromal.activated_fap]
    neut = df[df.neutrophil_like]
    if len(active) < 2 or len(neut) < 2:
        return None, None
    tree = cKDTree(neut[["x", "y"]].to_numpy())
    obs_dist = tree.query(active[["x", "y"]].to_numpy(), k=1)[0]
    obs = float(np.median(obs_dist))
    null = np.empty(n_perm)
    stromal_xy = stromal[["x", "y"]].to_numpy()
    for i in range(n_perm):
        take = rng.choice(len(stromal_xy), size=len(active), replace=False)
        null[i] = np.median(tree.query(stromal_xy[take], k=1)[0])
    result = {
        "n_stromal": len(stromal), "n_active_fap": len(active), "n_neutrophil_like": len(neut),
        "observed_median_nearest_distance": obs,
        "null_median": float(np.median(null)),
        "distance_ratio_observed_to_null": float(obs / np.median(null)),
        "empirical_p_closer": float((1 + np.sum(null <= obs)) / (n_perm + 1)),
    }
    return result, null


rng = np.random.default_rng(3818)
frames, summaries, nulls = [], [], {}
for section, stem in SECTIONS.items():
    df = load_section(section, stem)
    frames.append(df)
    test, null = spatial_test(df, rng)
    if null is not None:
        nulls[section] = null
    stromal = df[df.cell_type.eq("Stromal")]
    myeloid = df[df.cell_type.eq("Myeloid")]
    row = {
        "section": section,
        "cells_subsampled": len(df),
        "stromal_cells": len(stromal),
        "myeloid_cells": len(myeloid),
        "activated_fap_n": int(stromal.activated_fap.sum()),
        "activated_fap_percent_stromal": float(100 * stromal.activated_fap.mean()),
        "alpha5_positive_percent_activated_fap": float(100 * stromal.loc[stromal.activated_fap, "alpha5_fap"].mean()),
        "neutrophil_like_n": int(myeloid.neutrophil_like.sum()),
        "neutrophil_like_percent_myeloid": float(100 * myeloid.neutrophil_like.mean()),
        "recruitment_score_activated_fap": float(stromal.loc[stromal.activated_fap, "recruitment_score"].mean()),
        "feedback_score_neutrophil_like": float(myeloid.loc[myeloid.neutrophil_like, "feedback_score"].mean()),
    }
    if test:
        row.update(test)
    summaries.append(row)

summary = pd.DataFrame(summaries)
summary.to_csv(OUT / "SCP3818_spatial_summary.csv", index=False)
(OUT / "SCP3818_spatial_summary.json").write_text(summary.to_json(orient="records", indent=2), encoding="utf-8")

fig = plt.figure(figsize=(12, 7.2), constrained_layout=True)
gs = fig.add_gridspec(2, 4, width_ratios=[1.35, 1.35, 1, 1])
colors = {"UC1 inflamed": "#D55E00", "UC1 less inflamed": "#4C78A8"}
for i, (section, df) in enumerate(zip(SECTIONS, frames)):
    ax = fig.add_subplot(gs[:, i])
    ax.scatter(df.x, df.y, s=0.12, c="#D8D8D8", rasterized=True, linewidths=0)
    n = df[df.neutrophil_like]
    a = df[df.activated_fap]
    ax.scatter(n.x, n.y, s=1.2, c="#0072B2", alpha=.75, rasterized=True, linewidths=0, label="neutrophil-like")
    ax.scatter(a.x, a.y, s=2.0, c="#D55E00", alpha=.85, rasterized=True, linewidths=0, label="activated FAP+ stromal")
    ax.set_title(section.replace("UC1 ", ""), loc="left", fontweight="bold")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
    if i == 0:
        ax.legend(frameon=False, markerscale=3.5, loc="lower left", fontsize=8)

ax = fig.add_subplot(gs[0, 2])
metrics = ["activated_fap_percent_stromal", "neutrophil_like_percent_myeloid"]
labels = ["activated FAP+\n(% stromal)", "neutrophil-like\n(% myeloid)"]
x = np.arange(2)
for j, row in summary.iterrows():
    ax.plot(x, row[metrics].astype(float), "o-", color=colors[row.section], label=row.section.replace("UC1 ", ""), lw=2)
ax.set_xticks(x, labels)
ax.set_ylabel("Cells (%)")
ax.legend(frameon=False, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)

ax = fig.add_subplot(gs[0, 3])
metric2 = ["recruitment_score_activated_fap", "feedback_score_neutrophil_like"]
labels2 = ["fibroblast\nrecruitment", "neutrophil\nfeedback"]
for _, row in summary.iterrows():
    ax.plot(x, row[metric2].astype(float), "o-", color=colors[row.section], lw=2)
ax.set_xticks(x, labels2)
ax.set_ylabel("Mean module percentile")
ax.spines[["top", "right"]].set_visible(False)

ax = fig.add_subplot(gs[1, 2:])
positions = []
for k, section in enumerate(SECTIONS):
    vals = nulls[section]
    pos = k + 1
    positions.append(pos)
    parts = ax.violinplot(vals, positions=[pos], widths=.65, showextrema=False)
    for body in parts["bodies"]:
        body.set_facecolor("#BDBDBD"); body.set_edgecolor("none"); body.set_alpha(.7)
    obs = summary.loc[summary.section.eq(section), "observed_median_nearest_distance"].iloc[0]
    p = summary.loc[summary.section.eq(section), "empirical_p_closer"].iloc[0]
    ax.scatter(pos, obs, s=55, c=colors[section], edgecolor="white", linewidth=.7, zorder=3)
    ax.text(pos, np.percentile(vals, 96), f"P={p:.3g}", ha="center", va="bottom", fontsize=8)
ax.set_xticks(positions, [s.replace("UC1 ", "") for s in SECTIONS])
ax.set_ylabel("Median nearest distance\n(activated FAP+ to neutrophil-like)")
ax.text(.01, .02, "Grey: 1,000 size-matched random stromal sets; colored: observed", transform=ax.transAxes, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle("SCP3818 Xenium: paired within-donor spatial validation", fontsize=14, fontweight="bold")
fig.savefig(OUT / "SCP3818_spatial_validation.png", dpi=300, bbox_inches="tight")
plt.close(fig)

print(summary.to_string(index=False))
