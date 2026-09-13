"""Compact, source-faithful Figure 3; all artwork and data marks are native vectors."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE / "source_data" / "figure3"
WIDTH, HEIGHT = 183.0, 168.0
INK, MUTED, GRID = "#25282C", "#5F666D", "#E1E5E8"
UC, BLOCKADE = "#C84C4C", "#168C88"
CMAP = LinearSegmentedColormap.from_list("figure3_effect", ["#286599", "#F7F7F5", "#B63848"])
STATES = ["CXCR4", "MX1/ISG", "OSM", "PADI4"]
PROGRAMS = ["OSM inflammation", "CXCR4 aging/retention", "Degranulation", "Immature granulopoiesis"]
FACTORIAL_CONTRASTS = ["alpha5beta1_fibroblast_interaction", "NAMPT_fibroblast_interaction", "alpha5beta1_NAMPT_synergy_UC"]


def text(fig, x, y, value, size=5.5, weight="normal", color=INK, **kwargs):
    return fig.text(x / WIDTH, 1 - y / HEIGHT, value, fontsize=size,
                    fontweight=weight, color=color, va=kwargs.pop("va", "top"), **kwargs)


def axis(fig, x, y, width, height, label):
    return fig.add_axes([x / WIDTH, 1 - (y + height) / HEIGHT,
                         width / WIDTH, height / HEIGHT], label=label)


def heading(fig, x, y, letter, caption):
    text(fig, x, y - .2, letter, size=7, weight="bold")
    text(fig, x + 4.5, y, caption, size=6.7, weight="bold")


def line(fig, x1, y1, x2, y2, color=GRID, lw=.5, **kwargs):
    fig.add_artist(Line2D([x1 / WIDTH, x2 / WIDTH], [1 - y1 / HEIGHT, 1 - y2 / HEIGHT],
                          transform=fig.transFigure, color=color, linewidth=lw, **kwargs))


def arrow(fig, start, end, color=MUTED, **kwargs):
    fig.add_artist(FancyArrowPatch((start[0] / WIDTH, 1 - start[1] / HEIGHT),
                                   (end[0] / WIDTH, 1 - end[1] / HEIGHT),
                                   transform=fig.transFigure, arrowstyle="-|>",
                                   mutation_scale=4.5, linewidth=.6, color=color,
                                   shrinkA=0, shrinkB=0, **kwargs))


def box(fig, x, y, width, height, fill="#F5F6F7", edge="#BAC1C7"):
    fig.add_artist(FancyBboxPatch((x / WIDTH, 1 - (y + height) / HEIGHT),
                                  width / WIDTH, height / HEIGHT,
                                  boxstyle="round,pad=0,rounding_size=0.004",
                                  transform=fig.transFigure, facecolor=fill,
                                  edgecolor=edge, linewidth=.55))


def draw_design(fig):
    heading(fig, 3, 3.0, "a", "Matched culture and blockade design")
    box(fig, 8, 10, 32, 10.0)
    text(fig, 24, 12.0, "Peripheral-blood\nneutrophils", size=6.1,
         weight="bold", ha="center", linespacing=1.3)
    box(fig, 52, 10, 34, 10.0, fill="#FFFFFF")
    text(fig, 69, 12.1, "Pre-culture reference\nIndependent donors", size=5.3,
         ha="center", linespacing=1.4)
    arrow(fig, (40.5, 15), (50.7, 15), linestyle=(0, (2, 2)))
    line(fig, 24, 20.0, 24, 22.1, color=MUTED, lw=.6)
    line(fig, 6.0, 22.1, 24, 22.1, color=MUTED, lw=.6)
    line(fig, 6.0, 22.1, 6.0, 45.0, color=MUTED, lw=.6)
    arms = [
        (24, 6.0, "Neutrophils only", "Untreated / NAMPTi / α5β1i", "#F4F2F7"),
        (32, 6.0, "Control FB", "Untreated coculture", "#F2F6F7"),
        (40, 10.0, "UC FB", "Untreated / α5β1i\nNAMPTi / dual blockade", "#FBF0EE"),
    ]
    for y, h, group, treatments, fill in arms:
        box(fig, 11, y, 75, h, fill=fill)
        arrow(fig, (6.0, y + h / 2), (10.0, y + h / 2))
        text(fig, 13, y + h / 2, group, size=5.8, weight="bold", va="center")
        line(fig, 39, y + 1.1, 39, y + h - 1.1, color="#CED3D7", lw=.45)
        text(fig, 41, y + h / 2, treatments, size=5.4, va="center", linespacing=1.25)
    text(fig, 8, 53.0, "Recovered neutrophils → transcriptomic profiling", size=5.4,
         weight="bold")


def matrix(fig, x, y, w, h, values, limit, label, significant=None):
    ax = axis(fig, x, y, w, h, label)
    arr = np.asarray(values, dtype=float)
    assert np.isfinite(arr).all() and np.abs(arr).max() <= limit
    im = ax.pcolormesh(np.arange(arr.shape[1] + 1) - .5,
                       np.arange(arr.shape[0] + 1) - .5, arr,
                       shading="flat", cmap=CMAP, norm=Normalize(-limit, limit),
                       edgecolors="white", linewidth=.45, rasterized=False)
    ax.set_xlim(-.5, arr.shape[1] - .5)
    ax.set_ylim(arr.shape[0] - .5, -.5)
    ax.set_axis_off()
    for row in range(arr.shape[0]):
        for col in range(arr.shape[1]):
            val = arr[row, col]
            star = "*" if significant is not None and significant[row, col] else ""
            rgba = im.cmap(im.norm(val))
            rgb = np.asarray(rgba[:3])
            linear = np.where(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)
            luminance = float(linear @ np.array([.2126, .7152, .0722]))
            ax.text(col, row, f"{val:+.2f}{star}", ha="center", va="center",
                    fontsize=6.0, color="white" if luminance < .179 else "black")
    return im


def colorbar(fig, im, x, y, w, limit, label):
    cax = axis(fig, x, y, w, 1.7, label)
    cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=[-limit, 0, limit])
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor("face")
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=1.5, width=.45, pad=1.0, labelsize=5.0, colors=MUTED)
    return cb


def draw_factorial(fig):
    heading(fig, 95, 3.0, "b", "Factorial effects on neutrophil states")
    text(fig, 99.5, 7.2, "6 complete donors; *FDR < 0.05", size=5.2, color=MUTED)
    data = pd.read_csv(DATA / "Figure_3F_factorial_state_abundance.csv")
    assert data.n_complete_donors.eq(6).all() and len(data) == 12
    p = data.pivot(index="state_label", columns="contrast_id", values="effect").loc[STATES, FACTORIAL_CONTRASTS]
    q = data.pivot(index="state_label", columns="contrast_id", values="FDR").loc[STATES, FACTORIAL_CONTRASTS]
    for i, caption in enumerate(["FB-dependent\nα5β1i", "FB-dependent\nNAMPTi", "α5β1i × NAMPTi\ninteraction"]):
        text(fig, 121 + i * 23, 12.0, caption, size=5.5, ha="center", linespacing=1.25)
    for i, state in enumerate(STATES):
        text(fig, 108, 19.0 + (i + .5) * 7.2, state, size=5.8, ha="right", va="center")
    im = matrix(fig, 109.5, 19, 69, 28.8, p, 1.3, "factorial", q.to_numpy() < .05)
    text(fig, 109.5, 51, "CLR effect", size=5.2, color=MUTED)
    colorbar(fig, im, 144, 51, 33, 1.3, "factorial_colorbar")


def draw_programs(fig):
    heading(fig, 3, 61.0, "c", "α5β1-responsive functional programs")
    text(fig, 7.5, 65.2, "α5β1i versus UC; 5 donor pairs; *FDR < 0.05", size=5.1, color=MUTED)
    data = pd.read_csv(DATA / "Figure_3E_restored_functional_programs.csv")
    selected = data.loc[data.contrast_id.eq("UA5_vs_UF") & data.program_label.isin(PROGRAMS)]
    assert len(selected) == 16 and selected.n_pairs.eq(5).all()
    p = selected.pivot(index="program_label", columns="state_label", values="logFC").loc[PROGRAMS, STATES]
    q = selected.pivot(index="program_label", columns="state_label", values="adj.P.Val").loc[PROGRAMS, STATES]
    for i, state in enumerate(STATES):
        text(fig, 43.5 + i * 12.5, 70.9, state, size=5.5, ha="center")
    labels = ["OSM inflammation", "CXCR4 aging/\nretention", "Degranulation", "Immature\ngranulopoiesis"]
    for i, caption in enumerate(labels):
        text(fig, 35.7, 75 + (i + .5) * 7.6, caption, size=5.6, ha="right",
             va="center", linespacing=1.15)
    im = matrix(fig, 37.25, 75, 50, 30.4, p, .35, "programs", q.to_numpy() < .05)
    text(fig, 7.5, 109.0, "Program-score effect", size=5.1, color=MUTED)
    colorbar(fig, im, 53, 109, 33, .35, "programs_colorbar")


def draw_reversal(fig):
    heading(fig, 95, 61.0, "d", "Blockade-associated signalling reversal")
    text(fig, 99.5, 65.2, "Predicted fibroblast-to-neutrophil routes", size=5.1, color=MUTED)
    data = pd.read_csv(DATA / "Figure_3J_fibroblast_neutrophil_signaling_reversals.csv")
    assert len(data) == 12
    p = data.pivot(index="route_label", columns="contrast_id", values="consensus_change")
    p = p.loc[:, ["UF_vs_CF", "UA5_vs_UF"]].sort_values("UF_vs_CF")
    routes = data.drop_duplicates("route_label").set_index("route_label").loc[p.index]
    for i, caption in enumerate(["UC versus\ncontrol", "α5β1i\nversus UC"]):
        text(fig, 143.225 + i * 23.45, 69.2, caption, size=5.4, ha="center", linespacing=1.1)
    for i, (_, route) in enumerate(routes.iterrows()):
        cy = 75 + (i + .5) * (30.4 / 6)
        group = "ECM-FB" if "[ECM-FB]" in route.name else "CCN-FB"
        state = next(s for s in STATES if s in route.neutrophil_state)
        text(fig, 99.5, cy - 2.0, f"{route.ligand} → {route.receptor}",
             size=5.7, fontstyle="italic")
        text(fig, 99.5, cy + .1, f"{group} → {state}", size=5.0, color=MUTED)
    im = matrix(fig, 131.5, 75, 46.9, 30.4, p, .9, "reversal")
    text(fig, 99.5, 109.0, "LR consensus change", size=5.1, color=MUTED)
    colorbar(fig, im, 145, 109, 32, .9, "reversal_colorbar")


def draw_trajectories(fig):
    heading(fig, 3, 116.0, "e", "α5β1-responsive gene expression along pseudotime")
    # An external, shared legend prevents condition labels covering the curves.
    for x, color, caption in [(115, UC, "UC fibroblasts"), (145.5, BLOCKADE, "UC + α5β1i")]:
        line(fig, x, 117.2, x + 4, 117.2, color=color, lw=1)
        text(fig, x + 5.1, 117.2, caption, size=5.3, va="center")
    text(fig, 7.5, 120.2, "Mean donor-by-bin expression ± s.e.", size=5.1, color=MUTED)
    data = pd.read_csv(DATA / "Figure_3I_alpha5beta1_gene_trajectories.csv")
    assert len(data) == 48
    for index, gene in enumerate(["PRIM1", "GRIP1", "COMMD1"]):
        left = 15 + index * 57
        text(fig, left + 23.8, 123.5, gene, size=6.6, weight="bold", fontstyle="italic", ha="center")
        ax = axis(fig, left, 128, 47.6, 30.3, f"trajectory_{gene}")
        q = data.loc[data.gene.eq(gene)]
        lower, upper = np.inf, -np.inf
        for condition, color in [("UC fibroblasts", UC), ("UC + α5β1i", BLOCKADE)]:
            z = q.loc[q.condition_label.eq(condition)].sort_values("pt_bin")
            assert len(z) == 8 and z.pt_bin.tolist() == list(range(1, 9))
            x, y, se = z.pt_bin.to_numpy(), z.mean_expression.to_numpy(), z.se_expression.to_numpy()
            ax.fill_between(x, y - se, y + se, color=color, alpha=.15, linewidth=0)
            ax.plot(x, y, color=color, linewidth=1.0, marker="o", markersize=1.6,
                    markeredgewidth=0, label=condition)
            lower = min(lower, float(np.min(y - se)))
            upper = max(upper, float(np.max(y + se)))
        margin = (upper - lower) * .075
        ax.set_ylim(max(0, lower - margin), upper + margin)
        ax.set_xlim(.85, 8.15)
        ax.set_xticks([1, 2, 4, 6, 8])
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10], min_n_ticks=3))
        ax.tick_params(length=2, width=.5, pad=1.5, labelsize=5.1, colors=INK)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_linewidth(.55)
        ax.spines[["left", "bottom"]].set_color("#60666D")
        ax.set_axisbelow(True)
        ax.grid(axis="y", color=GRID, linewidth=.35)
        ax.set_xlabel("Pseudotime bin", fontsize=5.4, labelpad=1.6)
        if index == 0:
            ax.set_ylabel("Log-normalized expression", fontsize=5.4, labelpad=2)


def build_figure3(save):
    mpl.rcParams.update({"font.family": "Arial", "font.size": 5.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42,
                         "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(WIDTH / 25.4, HEIGHT / 25.4), facecolor="white")
    draw_design(fig)
    draw_factorial(fig)
    draw_programs(fig)
    draw_reversal(fig)
    draw_trajectories(fig)
    save(fig, 3)
