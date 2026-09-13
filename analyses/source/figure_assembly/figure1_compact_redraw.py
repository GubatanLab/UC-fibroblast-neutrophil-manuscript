"""Native, compact Figure 1 redraw. Original coordinates and result estimates only."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "source_data" / "figure1"
WIDTH, HEIGHT = 183.0, 170.0
INK, MUTED, GRID = "#25282C", "#5F666D", "#E4E7EB"
FB, NEU = "#D9782D", "#6C4C9A"
NONINFL, INFL = "#2C7F7B", "#C73E4D"
FIB_COLORS = {
    "Activated crypt-top fibroblast": "#8E63A9",
    "Crypt-bottom fibroblast": "#62A878",
    "Inflammatory fibroblast": "#C73E4D",
    "LP fibroblast": "#D9903D",
    "Resting crypt-top fibroblast": "#4D89A8",
}
NEU_COLORS = {
    "Neutrophil CXCR4": "#6C4C9A", "Neutrophil MX1": "#367FA3",
    "Neutrophil OSM": "#C73E4D", "Neutrophil PADI4": "#D9903D",
    "Neutrophil lowRNA": "#A5A5A5",
}


def text(fig, x, y, value, size=5.5, weight="normal", color=INK, **kwargs):
    return fig.text(x / WIDTH, 1 - y / HEIGHT, value, fontsize=size,
                    fontweight=weight, color=color, va=kwargs.pop("va", "top"), **kwargs)


def axis(fig, x, y, width, height):
    return fig.add_axes([x / WIDTH, 1 - (y + height) / HEIGHT,
                         width / WIDTH, height / HEIGHT])


def heading(fig, x, y, letter, caption):
    text(fig, x, y - .2, letter, size=7, weight="bold")
    text(fig, x + 4.5, y, caption, size=6.8, weight="bold")


def line(fig, x1, y1, x2, y2, color=GRID, lw=.5, **kwargs):
    fig.add_artist(Line2D([x1 / WIDTH, x2 / WIDTH],
                          [1 - y1 / HEIGHT, 1 - y2 / HEIGHT],
                          transform=fig.transFigure, color=color, linewidth=lw, **kwargs))


def style_axis(ax, grid="x"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["bottom"].set_color("#60666D")
    ax.spines["left"].set_color("#60666D")
    ax.tick_params(length=2, width=.5, pad=1.5, labelsize=5.2, colors=INK)
    ax.set_axisbelow(True)
    ax.grid(axis=grid, color=GRID, linewidth=.35)


def mini_umap_axes(ax):
    for end in [(.26, .10), (.10, .30)]:
        ax.annotate("", xy=end, xytext=(.10, .10), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=.55,
                                    mutation_scale=3.8, shrinkA=0, shrinkB=0))
    ax.text(.20, .015, "UMAP1", transform=ax.transAxes, fontsize=5,
            ha="center", va="bottom", color=MUTED)
    ax.text(.015, .23, "UMAP2", transform=ax.transAxes, fontsize=5,
            ha="left", va="center", rotation=90, color=MUTED)


def state_legend(fig, left, items, colors):
    # Two compact rows; the exact same state colors are used in the forest plot.
    for index, (key, caption) in enumerate(items):
        row, col = divmod(index, 3)
        x, y = left + col * 28.5, 40.5 + row * 3.0
        fig.add_artist(Line2D([x / WIDTH], [1 - y / HEIGHT],
                              marker="o", markersize=2.6, linestyle="", color=colors[key],
                              transform=fig.transFigure))
        text(fig, x + 1.8, y, caption, size=5.0, va="center")


def draw_umaps(fig):
    heading(fig, 3.0, 2.7, "a", "Human UC fibroblast and neutrophil state landscape")
    conditions = [("Healthy", "Control"), ("Uninflamed UC", "UC noninflamed"),
                  ("Inflamed UC", "UC inflamed")]
    audit = {}
    for lineage, left, colors, size in [("fibroblast", 7.5, FIB_COLORS, .9),
                                       ("neutrophil", 97.0, NEU_COLORS, .32)]:
        data = pd.read_csv(DATA / f"{lineage}_existing_umap.csv")
        assert len(data) == (5187 if lineage == "fibroblast" else 38684)
        assert not data.cell.duplicated().any() and set(data.state) == set(colors)
        assert set(data.condition) == {c[0] for c in conditions}
        text(fig, left, 7.2, lineage.capitalize() + "s", size=6.2, weight="bold")
        text(fig, left + 82, 7.2, f"{len(data):,} cells", size=5.1, ha="right", color=MUTED)
        x0, x1 = data.UMAP1.min(), data.UMAP1.max()
        y0, y1 = data.UMAP2.min(), data.UMAP2.max()
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        span_y = max((y1 - y0) * 1.16, (x1 - x0) * 1.16 / (27.5 / 22.21))
        span_x = span_y * (27.5 / 22.21)
        cx -= span_x * .05
        cy -= span_y * .05
        # Shared limits and equal distance scaling across every condition within each lineage.
        for index, (condition, caption) in enumerate(conditions):
            subset = data.loc[data.condition.eq(condition)].sample(frac=1, random_state=8201)
            x = left + index * 28.2
            text(fig, x + 13.75, 10.9, caption, size=5.4, ha="center", weight="bold")
            text(fig, x + 13.75, 13.6, f"n = {len(subset):,}", size=5.0, ha="center", color=MUTED)
            ax = axis(fig, x, 16.1, 27.5, 22.21)
            ax.scatter(subset.UMAP1, subset.UMAP2, c=subset.state.map(colors), s=size,
                       alpha=.88, linewidths=0, rasterized=True)
            ax.set_xlim(cx - span_x / 2, cx + span_x / 2)
            ax.set_ylim(cy - span_y / 2, cy + span_y / 2)
            ax.set_aspect("equal", adjustable="box")
            ax.axis("off")
            mini_umap_axes(ax)
        audit[lineage] = data.groupby("condition").size().to_dict()
    state_legend(fig, 8.0, [
        ("Activated crypt-top fibroblast", "Activated crypt-top"),
        ("Crypt-bottom fibroblast", "Crypt-bottom"),
        ("Inflammatory fibroblast", "Inflammatory"),
        ("LP fibroblast", "Lamina propria"),
        ("Resting crypt-top fibroblast", "Resting crypt-top"),
    ], FIB_COLORS)
    state_legend(fig, 97.5, [(s, s.replace("Neutrophil ", "").replace("lowRNA", "Low-RNA"))
                            for s in NEU_COLORS], NEU_COLORS)
    return audit


def draw_abundance(fig):
    heading(fig, 3.0, 47.0, "b", "Paired cell-state remodeling")
    text(fig, 7.5, 50.7, "Inflamed minus noninflamed; median and 95% CI", size=5.1, color=MUTED)
    data = pd.read_csv(DATA / "paired_cell_state_abundance_tests.csv").set_index("state")
    ax = axis(fig, 34, 54.5, 44, 34.5)
    ax.set_xlim(-44, 54)
    ax.set_ylim(89, 54.5)
    ax.set_xticks([-40, -20, 0, 20, 40])
    ax.set_yticks([])
    style_axis(ax)
    ax.spines["left"].set_visible(False)
    ax.axvline(0, color="#777F88", linewidth=.7)
    text(fig, 7.5, 53.6, "Fibroblasts | 13 patient pairs", size=5.1, weight="bold")
    text(fig, 7.5, 71.9, "Neutrophils | 17 patient pairs", size=5.1, weight="bold")
    text(fig, 84.9, 53.6, "q", size=5.2, ha="center", weight="bold")
    rows = [
        ("Inflammatory fibroblast", "Inflammatory", 57.2),
        ("Resting crypt-top fibroblast", "Resting crypt-top", 60.3),
        ("Activated crypt-top fibroblast", "Activated crypt-top", 63.4),
        ("Crypt-bottom fibroblast", "Crypt-bottom", 66.5),
        ("LP fibroblast", "Lamina propria", 69.6),
        ("Neutrophil OSM", "OSM", 75.6),
        ("Neutrophil CXCR4", "CXCR4", 78.7),
        ("Neutrophil MX1", "MX1", 81.8),
        ("Neutrophil PADI4", "PADI4", 84.9),
        ("Neutrophil lowRNA", "Low-RNA", 88.0),
    ]
    colors = FIB_COLORS | NEU_COLORS
    for state, caption, y in rows:
        row = data.loc[state]
        assert row.ci_low <= row.effect_pp <= row.ci_high
        text(fig, 8.5, y, caption, size=5.15, va="center")
        ax.hlines(y, row.ci_low, row.ci_high, color=colors[state], linewidth=.85)
        ax.scatter([row.effect_pp], [y], s=16, color=colors[state], edgecolors="white", linewidths=.35, zorder=3)
        text(fig, 84.9, y, f"{row.q_value:.3f}", size=5.0, ha="center", va="center",
             weight="bold" if row.q_value < .05 else "normal")
    text(fig, 56, 93.0, "Change in abundance (percentage points)", size=5.5, ha="center")


def draw_activation(fig):
    heading(fig, 96.0, 47.0, "c", "Inflammatory-fibroblast activation")
    text(fig, 100.5, 50.7, "Inflammatory vs other fibroblasts; 10 biopsies", size=5.1, color=MUTED)
    data = pd.read_csv(DATA / "inflammatory_fibroblast_target_activation_tests.csv").set_index("feature_ascii")
    rows = [("FAP", "FAP"), ("ITGA5", "ITGA5"), ("ITGB1", "ITGB1"),
            ("alpha5beta1 module", "α5β1 module"), ("FAP-alpha5beta1 module", "FAP-α5β1 module")]
    ax = axis(fig, 126.0, 54.5, 42.0, 24.0)
    ax.set_xlim(-.05, 1.25)
    ax.set_ylim(78.5, 54.5)
    ax.set_xticks([0, .4, .8, 1.2])
    ax.set_yticks([])
    style_axis(ax)
    ax.spines["left"].set_visible(False)
    ax.axvline(0, color="#777F88", linewidth=.7)
    text(fig, 176.0, 53.6, "q", size=5.2, ha="center", weight="bold")
    for y, (key, caption) in zip([57.2, 61.8, 66.4, 71.0, 75.6], rows):
        row = data.loc[key]
        assert row.ci_low <= row.effect <= row.ci_high and row.n_biopsies == 10
        text(fig, 101, y, caption, size=5.3, va="center")
        ax.hlines(y, row.ci_low, row.ci_high, color=INFL, linewidth=.9)
        ax.scatter([row.effect], [y], s=23, color=INFL, edgecolors="white", linewidths=.35, zorder=3)
        text(fig, 176, y, f"{row.q_value:.3f}", size=5.1, ha="center", va="center")
    text(fig, 145.0, 82.5, "Within-biopsy expression difference", size=5.5, ha="center")


def draw_priorities(fig):
    heading(fig, 3.0, 98.0, "d", "Prioritized reciprocal signals")
    data = pd.read_csv(DATA / "figure1_prioritized_interactions.csv")
    text(fig, 7.5, 102.2, "Fibroblast → neutrophil", size=5.2, weight="bold")
    text(fig, 7.5, 117.1, "Neutrophil → fibroblast", size=5.2, weight="bold")
    ax = axis(fig, 35.0, 104.5, 44, 25.1)
    ax.set_xlim(.70, .93)
    ax.set_ylim(130, 104.5)
    ax.set_xticks([.70, .80, .90])
    ax.set_yticks([])
    style_axis(ax)
    ax.spines["left"].set_visible(False)
    items = [("HGF", "CD44", 106.8, FB), ("CXCL1", "CXCR2", 110.3, FB),
             ("CSF3", "CSF3R", 113.8, FB), ("OSM", "OSMR", 120.5, NEU),
             ("IL1B", "IL1R1", 124.0, NEU), ("NAMPT", "ITGA5", 127.5, NEU),
             ("NAMPT", "ITGB1", 127.5, NEU)]
    plotted = []
    for ligand, receptor, y, color in items:
        subset = data.loc[data.ligand.eq(ligand) & data.receptor.eq(receptor)]
        assert len(subset) == 1
        row = subset.iloc[0]
        ax.hlines(y, .70, row.prioritization_score, color="#C8CDD3", linewidth=.5)
        ax.scatter([row.prioritization_score], [y], s=31 * row.fraction_expressing_ligand_receptor,
                   color=color, edgecolors="white", linewidths=.35, zorder=3)
        if ligand != "NAMPT":
            text(fig, 8.5, y, f"{ligand} → {receptor}", size=5.2, va="center")
        else:
            ax.text(row.prioritization_score, y - 1.4, receptor, ha="center", va="center", fontsize=5.0)
        plotted.append({"ligand": ligand, "receptor": receptor, "score": float(row.prioritization_score)})
    text(fig, 8.5, 127.5, "NAMPT → α5β1*", size=5.2, va="center")
    text(fig, 56.5, 133.6, "MultiNicheNet prioritization score", size=5.4, ha="center")
    text(fig, 8.5, 137.4, "LR-expressing fraction", size=5.0, va="center", color=MUTED)
    for x, fraction in [(48, .56), (61, .78), (74, 1.0)]:
        fig.add_artist(Line2D([x / WIDTH], [1 - 137.4 / HEIGHT], marker="o", linestyle="",
                              markersize=np.sqrt(31 * fraction), color="#68717B", transform=fig.transFigure))
        text(fig, x + 2, 137.4, f"{fraction:.0%}", size=5.0, va="center")
    return plotted


def draw_pairs(fig):
    heading(fig, 96.0, 88.0, "e", "Paired patient-level circuit scores")
    text(fig, 100.5, 92.2, "9 patient pairs; median Δ and BH-adjusted q", size=5.1, color=MUTED)
    scores = pd.read_csv(DATA / "paired_lr_interaction_scores.csv")
    tests = pd.read_csv(DATA / "paired_lr_interaction_tests.csv").set_index("interaction")
    pairs = [("CXCL1 -> CXCR2", "CXCL1 → CXCR2"),
             ("NAMPT -> ITGA5+ITGB1", "NAMPT → α5β1*"), ("OSM -> OSMR", "OSM → OSMR")]
    for index, (key, caption) in enumerate(pairs):
        x = 107.5 + index * 24.6
        text(fig, x + 10.8, 95.3, caption, size=5.1, weight="bold", ha="center")
        stats = tests.loc[key]
        text(fig, x + 10.8, 98.0, f"Δ {stats.median_delta:+.2f}; q={stats.q_value:.3f}",
             size=5.0, ha="center", color=MUTED)
        ax = axis(fig, x, 101.7, 21.6, 29.0)
        wide = scores.loc[scores.interaction.eq(key)].pivot(index="patient", columns="group", values="score")
        assert len(wide) == 9 and not wide.isna().any().any()
        assert np.isclose(np.median(wide["Inflamed.UC"] - wide["Uninflamed.UC"]), stats.median_delta)
        for _, row in wide.iterrows():
            vals = [row["Uninflamed.UC"], row["Inflamed.UC"]]
            ax.plot([0, 1], vals, color="#AEB5BD", linewidth=.5, alpha=.85, zorder=1)
            ax.scatter([0, 1], vals, c=[NONINFL, INFL], s=9.5, edgecolors="white", linewidths=.25, zorder=3)
        ax.set_xlim(-.25, 1.25)
        low, high = wide.to_numpy().min(), wide.to_numpy().max()
        ax.set_ylim(low - (high - low) * .1, high + (high - low) * .1)
        ax.set_xticks([0, 1], ["Noninfl.", "Inflamed"])
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        style_axis(ax, grid="y")
        ax.tick_params(labelsize=5.0)
    text(fig, 101.0, 116.2, "LR co-expression score", size=5.3, rotation=90, ha="center", va="center")
    text(fig, 100.5, 137.0, "*α5β1 score: min(ITGA5, ITGB1); components shown in d.", size=5.0, color=MUTED)


def model_box(fig, x, y, width, height, fill, edge, title_value, subtitle):
    patch = FancyBboxPatch((x / WIDTH, 1 - (y + height) / HEIGHT), width / WIDTH, height / HEIGHT,
                           boxstyle="round,pad=0.0018,rounding_size=0.003",
                           transform=fig.transFigure, facecolor=fill, edgecolor=edge, linewidth=.65)
    fig.add_artist(patch)
    text(fig, x + width / 2, y + 1.8, title_value, size=5.5, weight="bold", ha="center")
    text(fig, x + width / 2, y + 5.2, subtitle, size=5.1, ha="center")


def model_arrow(fig, start, end, color, rad=0):
    fig.add_artist(FancyArrowPatch((start[0] / WIDTH, 1 - start[1] / HEIGHT),
                                   (end[0] / WIDTH, 1 - end[1] / HEIGHT),
                                   transform=fig.transFigure, arrowstyle="-|>",
                                   connectionstyle=f"arc3,rad={rad}", mutation_scale=5.2,
                                   linewidth=.9, color=color, shrinkA=0, shrinkB=0))


def draw_model(fig):
    heading(fig, 3.0, 143.0, "f", "Proposed circuit and perturbation tests")
    model_box(fig, 7.5, 149.0, 25, 8.7, "#F4F5F6", "#AEB5BD", "Inflamed UC", "ECM + inflammation")
    model_box(fig, 46, 149.0, 44, 8.7, "#FCEDDF", FB, "FAP+ / α5β1-high", "Inflammatory fibroblasts")
    model_box(fig, 130, 149.0, 48.5, 8.7, "#EEE7F5", NEU, "OSM · PADI4 · CXCR4 · MX1", "Neutrophil states")
    model_arrow(fig, (33, 153.3), (45, 153.3), "#8B929B")
    text(fig, 39, 149.8, "Remodeling", size=5.0, ha="center")
    model_arrow(fig, (90.7, 151.2), (129.3, 151.2), FB, rad=-.14)
    model_arrow(fig, (129.3, 155.5), (90.7, 155.5), NEU, rad=-.14)
    text(fig, 110, 146.0, "CXCL1 · CSF3 · HGF", size=5.0, ha="center")
    text(fig, 110, 159.4, "NAMPT* · OSM · IL1B", size=5.0, ha="center")
    for x, width, title_value, subtitle, color in [
        (8, 44, "FAP-cell ablation", "Niche-cell necessity", FB),
        (62, 44, "α5β1 blockade", "State transmission", NEU),
        (116, 62, "Extracellular NAMPT blockade", "Ligand dependence", NEU),
    ]:
        line(fig, x, 162.0, x + width, 162.0, color=color, lw=.7)
        text(fig, x + width / 2, 162.5, title_value + " | " + subtitle, size=5.0, ha="center")


def build_figure1(save_callback):
    with mpl.rc_context({"font.family": "Arial", "font.size": 5.5,
                         "axes.linewidth": .55, "pdf.fonttype": 42, "ps.fonttype": 42}):
        fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor="white")
        audit = {"page_mm": [WIDTH, HEIGHT], "canonical_outputs_modified": False,
                 "umap_recomputed": False, "statistics_recomputed": False,
                 "revision": "Larger UMAPs; compressed five-row activation panel; 58% taller paired plots",
                 "umap_panel_mm": [27.5, 22.21], "paired_plot_mm": [21.6, 29.0]}
        audit["umap_cells"] = draw_umaps(fig)
        draw_abundance(fig)
        draw_activation(fig)
        audit["prioritized_interactions"] = draw_priorities(fig)
        draw_pairs(fig)
        draw_model(fig)
        fig.canvas.draw()
        # Scientific labels are native PDF text, not flattened source-panel labels.
        audit["scientific_text_count"] = len(fig.texts) + sum(len(ax.texts) for ax in fig.axes)
        audit["axes_count"] = len(fig.axes)
        audit["source_tables"] = [p.name for p in sorted(DATA.glob("*.csv"))]
        (HERE / "Figure_1_REDRAW_NONCANONICAL_manifest.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
        save_callback(fig, 1)
