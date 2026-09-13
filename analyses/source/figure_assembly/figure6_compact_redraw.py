"""Figure 6: native, compact layout preserving the seven-panel redraw's values."""
from __future__ import annotations

import json
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "source_data" / "figure6"
WIDTH, HEIGHT = 183.0, 164.0
INK, MUTED, GRID = "#25282C", "#5F666D", "#E4E7EB"
PURPLE, TEAL, ORANGE, BLUE = "#7651A6", "#168C88", "#D97724", "#177AA9"
FEATURES = ["activated FAP fibroblast fraction", "FAP inflammatory", "alpha5beta1 adhesion",
            "neutrophil recruitment", "OSM response", "myeloid feedback proxy"]
PATHWAYS = ["TNF/NFkB", "Inflammatory response", "Degranulation", "Chemotaxis", "Phagocytosis",
            "Oxidative phosphorylation", "Adhesion/ECM"]
CONTRASTS = ["Human alpha5beta1i x UC fibroblast context", "Mouse alpha5beta1 blockade", "Mouse FAP ablation"]
ENDPOINTS = ["Cxcr4+ aged/retained neutrophil", "Osm+ inflammatory neutrophil",
             "Padi4+ antimicrobial/NET-competent neutrophil", "Mx1+ interferon-responsive neutrophil"]


def text(fig, x, y, value, size=5.5, weight="normal", color=INK, **kwargs):
    return fig.text(x / WIDTH, 1 - y / HEIGHT, value, fontsize=size, fontweight=weight,
                    color=color, va=kwargs.pop("va", "top"), **kwargs)


def axis(fig, x, y, width, height, name):
    return fig.add_axes([x / WIDTH, 1 - (y + height) / HEIGHT, width / WIDTH, height / HEIGHT], label=name)


def heading(fig, x, y, letter, caption):
    text(fig, x, y - .2, letter, size=7, weight="bold")
    text(fig, x + 4.5, y, caption, size=6.8, weight="bold")


def style_axis(ax, grid="x"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#60666D")
    ax.tick_params(length=2, width=.5, pad=1.5, labelsize=5.2, colors=INK)
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color=GRID, linewidth=.35)


def marker_key(fig, x, y, color, label, marker="o"):
    fig.add_artist(Line2D([x / WIDTH], [1 - y / HEIGHT], marker=marker,
                          markersize=3, color=color, linestyle="", transform=fig.transFigure))
    text(fig, x + 2.1, y, label, size=5.1, va="center")


def read_data():
    names = {
        "baseline": "taurus_uc_baseline_inflammation_patient_deltas.tsv",
        "stats": "taurus_uc_baseline_inflammation_statistics.tsv",
        "state": "taurus_uc_alpha5beta1_state_contrast_by_patient.tsv",
        "paired": "taurus_uc_patient_paired_deltas.tsv",
        "spatial": "SCP3818_spatial_summary.csv",
        "pathways": "cross_system_pathway_evidence.csv",
        "invivo": "cross_system_neutrophil_fibroblast_evidence.csv",
    }
    return {key: pd.read_csv(DATA / name, sep="\t" if name.endswith(".tsv") else ",")
            for key, name in names.items()}


def baseline_effects(baseline, stats):
    # Exact prior-redraw algorithm/seed/order; this is an unweighted sign effect.
    rng = np.random.default_rng(20260901)
    records = []
    for feature in FEATURES:
        values = baseline.loc[baseline.feature.eq(feature), "delta"].dropna().to_numpy(float)
        values = values[values != 0]
        effect = ((values > 0).sum() - (values < 0).sum()) / len(values)
        bootstraps = []
        for _ in range(1800):
            sample = values[rng.integers(0, len(values), len(values))]
            bootstraps.append(((sample > 0).sum() - (sample < 0).sum()) / len(sample))
        low, high = np.quantile(bootstraps, [.025, .975])
        q = stats.loc[stats.feature.eq(feature), "FDR"].iloc[0]
        records.append(dict(feature=feature, effect=float(effect), ci_low=float(low),
                            ci_high=float(high), q=float(q), n=int(len(values))))
    return records


def draw_chain(fig):
    heading(fig, 3, 2.7, "a", "Human discovery, independent validation and perturbation")
    boxes = [
        ("Discovery UC", "Disease activation", "#C73E4D"),
        ("TAURUS UC", "Independent cohort", PURPLE),
        ("Spatial UC", "Tissue proximity", BLUE),
        ("Human coculture", "Pathway blockade", "#39769B"),
        ("FAP ablation", "Niche necessity", PURPLE),
        ("α5β1 blockade", "State control", TEAL),
    ]
    for index, (title, subtitle, color) in enumerate(boxes):
        x, y, width, height = 7.5 + 28.6 * index, 8.3, 26.4, 8.1
        fig.add_artist(Rectangle((x / WIDTH, 1 - (y + height) / HEIGHT), width / WIDTH, height / HEIGHT,
                                  transform=fig.transFigure, facecolor=mpl.colors.to_rgba(color, .055),
                                  edgecolor="#D5D9DE", linewidth=.45))
        fig.add_artist(Line2D([x / WIDTH, (x + width) / WIDTH], [1 - y / HEIGHT] * 2,
                              transform=fig.transFigure, color=color, linewidth=1.1))
        text(fig, x + width / 2, y + 1.5, title, size=5.6, weight="bold", ha="center")
        text(fig, x + width / 2, y + 4.7, subtitle, size=5.0, ha="center")
        if index < 5:
            fig.add_artist(FancyArrowPatch(((x + width + .15) / WIDTH, 1 - (y + 4.1) / HEIGHT),
                                          ((x + width + 2.05) / WIDTH, 1 - (y + 4.1) / HEIGHT),
                                          transform=fig.transFigure, arrowstyle="-|>", mutation_scale=4,
                                          color=MUTED, linewidth=.55, shrinkA=0, shrinkB=0))


def draw_baseline(fig, data):
    heading(fig, 3, 21.5, "b", "TAURUS UC: paired inflammation")
    text(fig, 7.5, 25.3, "Inflamed vs noninflamed; 95% bootstrap CI", size=5.1, color=MUTED)
    records = baseline_effects(data["baseline"], data["stats"])
    labels = ["FAP+ fraction", "FAP inflammatory", "α5β1 adhesion", "Neutrophil recruitment", "OSM response", "Myeloid proxy"]
    ax = axis(fig, 39, 30.8, 50, 27.7, "b_effects")
    ax.set_xlim(-.045, 1.055)
    ax.set_ylim(58.5, 30.8)
    ax.set_xticks([0, .25, .5, .75, 1], ["0", "0.25", "0.50", "0.75", "1.00"])
    ax.set_yticks([])
    style_axis(ax)
    ax.spines["left"].set_visible(False)
    ax.axvline(0, color="#777F88", linewidth=.65)
    text(fig, 97, 29.0, "q", size=5.2, weight="bold", ha="center")
    text(fig, 106, 29.0, "n", size=5.2, weight="bold", ha="center")
    for y, label, row in zip([33.5, 38.1, 42.7, 47.3, 51.9, 56.5], labels, records):
        color = "#B1782B" if row["feature"] == "myeloid feedback proxy" else PURPLE
        text(fig, 8.5, y, label, size=5.3, va="center")
        ax.hlines(y, row["ci_low"], row["ci_high"], color=color, linewidth=.9)
        ax.scatter(row["effect"], y, s=20, color=color, edgecolors="white", linewidths=.3, zorder=3)
        text(fig, 97, y, f'{row["q"]:.3f}', size=5.1, ha="center", va="center",
             weight="bold" if row["q"] < .05 else "normal")
        text(fig, 106, y, str(row["n"]), size=5.1, ha="center", va="center")
    text(fig, 64, 62.1, "Paired directional effect", size=5.5, ha="center")
    text(fig, 7.5, 65.7, "Myeloid: monocyte proxy; TAURUS did not retain neutrophils.", size=5.0, color=MUTED)
    return records


def draw_state(fig, state):
    heading(fig, 112, 21.5, "c", "FAP+ fibroblast adhesion state")
    delta = (state.activated - state.other).median()
    text(fig, 116.5, 25.3, f"16 patients; median paired Δ = {delta:.2f}", size=5.1, color=MUTED)
    ax = axis(fig, 128, 31, 50, 27.5, "c_state")
    for _, row in state.iterrows():
        ax.plot([0, 1], [row.other, row.activated], color="#B2BBC5", linewidth=.5, alpha=.85)
    ax.scatter(np.zeros(len(state)), state.other, color="#929CA7", s=11, edgecolors="white", linewidths=.25, zorder=3)
    ax.scatter(np.ones(len(state)), state.activated, color=PURPLE, s=12, edgecolors="white", linewidths=.25, zorder=3)
    ax.plot([0, 1], [state.other.median(), state.activated.median()], color=INK, linewidth=1.3, zorder=4)
    ax.set_xlim(-.14, 1.14)
    ax.set_ylim(.53, .96)
    ax.set_xticks([0, 1], ["Other\nfibroblasts", "Activated FAP+\nfibroblasts"])
    ax.set_yticks([.6, .7, .8, .9])
    style_axis(ax, grid="y")
    text(fig, 119, 44.75, "α5β1 adhesion module", size=5.3, rotation=90, ha="center", va="center")
    return {"n_patients": len(state), "median_paired_difference": float(delta)}


def draw_spatial(fig, spatial):
    heading(fig, 3, 71, "d", "Paired spatial validation in UC")
    text(fig, 7.5, 75.1, "One donor; paired UC sections (SCP3818)", size=5.1, color=MUTED)
    inf = spatial.loc[spatial.section.eq("UC1 inflamed")].iloc[0]
    less = spatial.loc[spatial.section.eq("UC1 less inflamed")].iloc[0]
    ratios = [inf.activated_fap_percent_stromal / less.activated_fap_percent_stromal,
              inf.neutrophil_like_percent_myeloid / less.neutrophil_like_percent_myeloid,
              inf.distance_ratio_observed_to_null]
    ax = axis(fig, 39, 80, 52, 24, "d_spatial")
    ax.set_xlim(0, 3.15)
    ax.set_ylim(104, 80)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_yticks([])
    style_axis(ax)
    ax.spines["left"].set_visible(False)
    ax.axvline(1, color="#777F88", linewidth=.7, linestyle="--")
    rows = [("Activated FAP+", "Inflamed / less inflamed", ORANGE),
            ("Neutrophil-like", "Inflamed / less inflamed", BLUE),
            ("Cell proximity", "Observed / null distance", TEAL)]
    for y, value, (label, detail, color) in zip([83.5, 91.5, 99.5], ratios, rows):
        text(fig, 8.5, y - 1.4, label, size=5.3, va="center")
        text(fig, 8.5, y + 1.4, detail, size=5.0, color=MUTED, va="center")
        ax.hlines(y, 0, value, color=color, linewidth=1.1)
        ax.scatter(value, y, s=25, color=color, edgecolors="white", linewidths=.3, zorder=3)
        text(fig, 101.5, y, f"{value:.2f}" + ("×" if y != 99.5 else ""), size=5.5,
             weight="bold", ha="center", va="center")
    text(fig, 65, 107.6, "Ratio", size=5.5, ha="center")
    text(fig, 7.5, 111.0, f"Proximity: P = {inf.empirical_p_closer:.3f}; 1,000 size-matched null sets", size=5.0, color=MUTED)
    return {"ratios": [float(v) for v in ratios], "empirical_p": float(inf.empirical_p_closer)}


def draw_outcomes(fig, paired):
    heading(fig, 112, 71, "e", "Activated FAP+ cells after anti-TNF")
    text(fig, 116.5, 75.1, "Site-matched pre/post adalimumab", size=5.1, color=MUTED)
    z = paired.loc[paired.feature.eq(FEATURES[0])]
    counts = {}
    for index, (outcome, label, color) in enumerate([("Remission", "Remission", TEAL), ("Non_Remission", "Nonremission", ORANGE)]):
        sub = z.loc[z.Remission_status.eq(outcome)]
        x = 128 + index * 28.5
        text(fig, x + 10.75, 79, f"{label} (n = {len(sub)})", size=5.1, weight="bold", ha="center")
        ax = axis(fig, x, 83.4, 21.5, 20.6, f"e_{outcome}")
        for _, row in sub.iterrows():
            vals = np.array([row.pre, row.post]) * 100
            ax.plot([0, 1], vals, color=color, alpha=.48, linewidth=.55, marker="o", markersize=1.45)
        ax.plot([0, 1], np.array([sub.pre.median(), sub.post.median()]) * 100,
                color=color, linewidth=1.6, marker="o", markersize=2.8, zorder=4)
        ax.set_xlim(-.18, 1.18)
        ax.set_ylim(-.35, 13.1)
        ax.set_xticks([0, 1], ["Pre", "Post"])
        ax.set_yticks([0, 4, 8, 12])
        style_axis(ax, grid="y")
        if index:
            ax.tick_params(axis="y", labelleft=False, length=0)
            ax.spines["left"].set_visible(False)
        counts[outcome] = len(sub)
    text(fig, 119, 93.7, "Activated FAP+ (% stroma)", size=5.2, rotation=90, ha="center", va="center")
    text(fig, 116.5, 111, "Lines: patients; bold lines: group medians", size=5.0, color=MUTED)
    return counts


def draw_pathways(fig, pathways):
    heading(fig, 3, 117, "f", "Cross-system pathway responses")
    labels = ["Human coculture\nα5β1 blockade", "Mouse DSS\nα5β1 blockade", "Mouse DSS\nFAP ablation"]
    for index, label in enumerate(labels):
        text(fig, 51.5 + index * 19, 121.3, label, size=5.1, ha="center")
    mat = pathways.pivot(index="axis", columns="system_contrast", values="NES").reindex(index=PATHWAYS, columns=CONTRASTS)
    qs = pathways.pivot(index="axis", columns="system_contrast", values="FDR").reindex(index=PATHWAYS, columns=CONTRASTS)
    ax = axis(fig, 42, 128, 57, 25.2, "f_pathways")
    cmap = mpl.colormaps["RdBu_r"].copy()
    cmap.set_bad("#F2F3F5")
    im = ax.pcolormesh(np.ma.masked_invalid(mat.to_numpy()), cmap=cmap,
                       norm=TwoSlopeNorm(vmin=-2.8, vcenter=0, vmax=2.8), edgecolors="white", linewidth=.45)
    ax.set_ylim(7, 0)
    ax.set_xlim(0, 3)
    ax.set_axis_off()
    for row, label in enumerate(PATHWAYS):
        text(fig, 40.2, 128 + (row + .5) * 3.6,
             label.replace("TNF/NFkB", "TNF / NFκB"), size=5.1, ha="right", va="center")
        for col in range(3):
            value = mat.iloc[row, col]
            q = qs.iloc[row, col]
            caption = "NA" if pd.isna(value) else f"{value:.1f}" + ("*" if q < .05 else "")
            color = MUTED if pd.isna(value) else ("white" if abs(value) > 1.55 else INK)
            ax.text(col + .5, row + .5, caption, fontsize=5.4, ha="center", va="center", color=color)
    cbax = axis(fig, 17, 155.4, 25, 1.25, "f_colorbar")
    cb = fig.colorbar(im, cax=cbax, orientation="horizontal", ticks=[-2, 0, 2])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=1.3, width=.4, labelsize=5.0, pad=1)
    text(fig, 7.5, 155.2, "NES", size=5.1)
    text(fig, 50, 155.4, "*q < 0.05; NA: not reported", size=5.0, color=MUTED)
    return {"reported_NES_values": int(mat.notna().sum().sum()), "not_reported_cells": int(mat.isna().sum().sum())}


def draw_invivo(fig, invivo):
    heading(fig, 112, 117, "g", "In-vivo neutrophil-state control")
    marker_key(fig, 117, 122.1, PURPLE, "FAP ablation")
    marker_key(fig, 146.5, 122.1, TEAL, "α5β1 blockade", marker="D")
    d = invivo.loc[invivo.system.eq("Mouse DSS scRNA-seq")]
    ax = axis(fig, 128, 128, 50, 25.2, "g_invivo")
    ax.set_xlim(-110, 7)
    ax.set_ylim(-.45, 3.45)
    ax.set_xticks([-100, -50, 0])
    ax.set_yticks([3, 2, 1, 0], ["CXCR4", "OSM", "PADI4", "MX1"])
    style_axis(ax)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=2.5)
    ax.axvline(0, color="#777F88", linewidth=.65)
    records = []
    for experiment, color, offset, marker in [("FAP-TK GCV ablation", PURPLE, .11, "o"),
                                              ("alpha5beta1 blockade", TEAL, -.11, "D")]:
        sub = d.loc[d.experiment.eq(experiment)].set_index("endpoint").reindex(ENDPOINTS)
        y = np.arange(4)[::-1] + offset
        ax.hlines(y, 0, sub.effect, color=mpl.colors.to_rgba(color, .55), linewidth=1.05)
        ax.scatter(sub.effect, y, s=20, color=color, marker=marker, edgecolors="white", linewidths=.25, zorder=3)
        records += [{"experiment": experiment, "endpoint": endpoint, "effect": float(value)}
                    for endpoint, value in zip(ENDPOINTS, sub.effect)]
    text(fig, 153, 157.1, "Treatment - DSS (cells per 10,000)", size=5.3, ha="center")
    return records


def build_figure6(save_callback):
    with mpl.rc_context({"font.family": "Arial", "font.size": 5.5, "axes.linewidth": .55,
                         "pdf.fonttype": 42, "ps.fonttype": 42}):
        data = read_data()
        fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor="white")
        draw_chain(fig)
        audit = {"page_mm": [WIDTH, HEIGHT], "canonical_outputs_modified": False,
                 "baseline_algorithm_changed": False, "statistical_tests_recomputed": False,
                 "baseline_axis_label_correction": "Unweighted paired directional/sign effect, not rank-biserial"}
        audit["baseline_effects"] = draw_baseline(fig, data)
        audit["state_contrast"] = draw_state(fig, data["state"])
        audit["spatial"] = draw_spatial(fig, data["spatial"])
        audit["outcome_patient_counts"] = draw_outcomes(fig, data["paired"])
        audit["pathways"] = draw_pathways(fig, data["pathways"])
        audit["invivo_effects"] = draw_invivo(fig, data["invivo"])
        audit["source_tables"] = [p.name for p in sorted(DATA.iterdir()) if p.is_file()]
        fig.canvas.draw()
        (HERE / "Figure_6_REDRAW_NONCANONICAL_manifest.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
        save_callback(fig, 6)
