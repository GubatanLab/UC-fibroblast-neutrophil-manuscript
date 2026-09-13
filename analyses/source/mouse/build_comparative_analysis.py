"""Comparative FAP-ablation versus alpha5beta1-blockade analysis.

This is a noncanonical, replicate-aware comparison built from existing mouse
colitis differential-abundance summaries and matched MultiNicheNet 2.1.0
outputs. It does not overwrite manuscript figures.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import mannwhitneyu, spearmanr


HERE = Path(__file__).resolve().parent
WORK = Path(__file__).resolve().parents[3]
SOURCE = WORK / "source_data/mouse_comparison/inputs"
TABLES = WORK / "results/mouse_comparison/tables"
FIGURES = WORK / "results/mouse_comparison/figures"

INK = "#25282C"
MUTED = "#646B73"
GRID = "#E2E6E9"
ABLATION = "#7651A6"
BLOCKADE = "#168C88"
EPITHELIAL = "#D97724"
IMMUNE = "#3B6FB6"
STROMAL = "#C84C4C"
COMPARTMENT_COLORS = {"Epithelial": EPITHELIAL, "Immune": IMMUNE, "Stromal": STROMAL}
INTERVENTION_ORDER = ["FAP ablation", "alpha5beta1 blockade"]
INTERVENTION_LABELS = {
    "FAP ablation": "FAP ablation",
    "alpha5beta1 blockade": "α5β1 blockade",
}
INTERVENTION_COLORS = {"FAP ablation": ABLATION, "alpha5beta1 blockade": BLOCKADE}
FOCUSED_ROUTES = [
    ("Fibroblast -> Neutrophil", "Il1b", "Il1rap"),
    ("Fibroblast -> Neutrophil", "Il1b", "Il1r2"),
    ("Fibroblast -> Neutrophil", "Il33", "Il1rap"),
    ("Fibroblast -> Neutrophil", "Cxcl1", "Cxcr2"),
    ("Fibroblast -> Neutrophil", "Cxcl2", "Cxcr2"),
    ("Fibroblast -> Neutrophil", "Cxcl5", "Cxcr2"),
    ("Fibroblast -> Neutrophil", "Csf1", "Csf1r"),
    ("Fibroblast -> Neutrophil", "Fn1", "Itga5"),
    ("Fibroblast -> Neutrophil", "Icam1", "Itgb2"),
    ("Neutrophil -> Fibroblast", "Osm", "Osmr"),
    ("Neutrophil -> Fibroblast", "Osm", "Il6st"),
    ("Neutrophil -> Fibroblast", "Il1b", "Il1r1"),
    ("Neutrophil -> Fibroblast", "Vegfa", "Nrp1"),
    ("Neutrophil -> Fibroblast", "Nampt", "Itga5"),
    ("Neutrophil -> Fibroblast", "Thbs1", "Itga4"),
    ("Neutrophil -> Fibroblast", "Tgm2", "Itga4"),
    ("Neutrophil -> Fibroblast", "S100a9", "Alcam"),
    ("Neutrophil -> Fibroblast", "Col1a1", "Itga5"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_hashes():
    paths = []
    for number in range(1, 7):
        folder = WORK / f"Figure {number}"
        for suffix in ("*.pdf", "*.png", "*.tiff"):
            paths.extend(sorted(folder.glob(suffix)))
    return {str(path.relative_to(WORK)): sha256(path) for path in paths}


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    result = np.full(values.shape, np.nan)
    keep = np.isfinite(values)
    p = values[keep]
    if not len(p):
        return result
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[keep] = restored
    return result


def bootstrap_mean_difference(reference, treatment, rng, draws=10000):
    reference = np.asarray(reference, dtype=float)
    treatment = np.asarray(treatment, dtype=float)
    observed = 100 * (treatment.mean() - reference.mean())
    boot = np.empty(draws)
    for index in range(draws):
        a = reference[rng.integers(0, len(reference), len(reference))]
        b = treatment[rng.integers(0, len(treatment), len(treatment))]
        boot[index] = 100 * (b.mean() - a.mean())
    return observed, *np.quantile(boot, [.025, .975])


def bootstrap_treatment_difference(ablation, blockade, rng, draws=10000):
    ablation = np.asarray(ablation, dtype=float)
    blockade = np.asarray(blockade, dtype=float)
    observed = 100 * (ablation.mean() - blockade.mean())
    boot = np.empty(draws)
    for index in range(draws):
        a = ablation[rng.integers(0, len(ablation), len(ablation))]
        b = blockade[rng.integers(0, len(blockade), len(blockade))]
        boot[index] = 100 * (a.mean() - b.mean())
    return observed, *np.quantile(boot, [.025, .975])


def clean_axis(ax, grid="x"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#32363A")
    ax.spines[["left", "bottom"]].set_linewidth(.55)
    ax.tick_params(labelsize=5.3, width=.5, length=2.2, pad=1.8, colors=INK)
    if grid in ("x", "both"):
        ax.grid(axis="x", color=GRID, linewidth=.5, zorder=0)
    if grid in ("y", "both"):
        ax.grid(axis="y", color=GRID, linewidth=.5, zorder=0)
    ax.set_axisbelow(True)


def panel_label(fig, x, y, letter, title):
    fig.text(x, y, letter, fontsize=8.0, fontweight="bold", color=INK, va="top")
    fig.text(x + .025, y, title, fontsize=7.0, fontweight="bold", color=INK, va="top")


def calculate_da(data, minimum_total_cells=None, seed=20260904):
    local = data.copy()
    if minimum_total_cells is not None:
        local = local.loc[local.total_cells >= minimum_total_cells]
    rng = np.random.default_rng(seed)
    records = []
    scopes = ["Compartment recovery", "Level 3 fraction within compartment"]
    for intervention in INTERVENTION_ORDER:
        block = local.loc[local.contrast_short.eq(intervention)]
        for scope in scopes:
            scoped = block.loc[block.feature_scope.eq(scope)]
            keys = scoped[["compartment", "feature"]].drop_duplicates()
            for key in keys.itertuples(index=False):
                q = scoped.loc[scoped.compartment.eq(key.compartment) & scoped.feature.eq(key.feature)]
                reference = q.loc[q.treatment_status.eq("Reference"), "proportion"].dropna().to_numpy(float)
                treatment = q.loc[q.treatment_status.eq("Treatment"), "proportion"].dropna().to_numpy(float)
                if len(reference) < 2 or len(treatment) < 2:
                    continue
                effect, low, high = bootstrap_mean_difference(reference, treatment, rng)
                test = mannwhitneyu(treatment, reference, alternative="two-sided", method="asymptotic")
                records.append({
                    "intervention": intervention,
                    "feature_scope": scope,
                    "compartment": key.compartment,
                    "feature": key.feature,
                    "reference_n": len(reference),
                    "treatment_n": len(treatment),
                    "reference_mean": reference.mean(),
                    "treatment_mean": treatment.mean(),
                    "effect_percentage_points": effect,
                    "ci_low": low,
                    "ci_high": high,
                    "mann_whitney_u": float(test.statistic),
                    "p_value": float(test.pvalue),
                    "minimum_total_cells": minimum_total_cells or 0,
                })
    result = pd.DataFrame(records)
    result["fdr_global_scope"] = np.nan
    for _, idx in result.groupby(["intervention", "feature_scope"]).groups.items():
        result.loc[idx, "fdr_global_scope"] = bh_adjust(result.loc[idx, "p_value"])
    return result


def comparative_da_table(primary, raw):
    left = primary.loc[primary.intervention.eq("FAP ablation")].copy()
    right = primary.loc[primary.intervention.eq("alpha5beta1 blockade")].copy()
    keys = ["feature_scope", "compartment", "feature"]
    keep = keys + ["effect_percentage_points", "ci_low", "ci_high", "fdr_global_scope",
                   "reference_n", "treatment_n"]
    merged = left[keep].merge(right[keep], on=keys, suffixes=("_ablation", "_blockade"))
    rng = np.random.default_rng(20260905)
    direct = []
    for row in merged.itertuples(index=False):
        block = raw.loc[
            raw.feature_scope.eq(row.feature_scope) & raw.compartment.eq(row.compartment) &
            raw.feature.eq(row.feature) & raw.treatment_status.eq("Treatment")
        ]
        a = block.loc[block.contrast_short.eq("FAP ablation"), "proportion"].dropna().to_numpy(float)
        b = block.loc[block.contrast_short.eq("alpha5beta1 blockade"), "proportion"].dropna().to_numpy(float)
        if len(a) >= 2 and len(b) >= 2:
            direct.append(bootstrap_treatment_difference(a, b, rng))
        else:
            direct.append((np.nan, np.nan, np.nan))
    direct = np.asarray(direct)
    merged["ablation_minus_blockade_percentage_points"] = direct[:, 0]
    merged["difference_ci_low"] = direct[:, 1]
    merged["difference_ci_high"] = direct[:, 2]
    merged["effect_concordant"] = np.sign(merged.effect_percentage_points_ablation) == np.sign(
        merged.effect_percentage_points_blockade)
    merged["max_abs_effect"] = merged[["effect_percentage_points_ablation",
                                       "effect_percentage_points_blockade"]].abs().max(axis=1)
    return merged


def milo_overview():
    ablation = pd.read_csv(SOURCE / "FAP_ablation_miloR_summary.csv")
    ablation = ablation.loc[ablation.contrast.eq("FAP_ablation")].copy()
    ablation["intervention"] = "FAP ablation"
    ablation = ablation.rename(columns={"bh_fdr_lt_0_05": "significant"})
    blockade = pd.read_csv(SOURCE / "a5B1_blockade_miloR_summary.csv")
    blockade = blockade.loc[blockade.contrast.eq("Blockade - severe DSS colitis")].copy()
    blockade["intervention"] = "alpha5beta1 blockade"
    blockade = blockade.rename(columns={"fdr_lt_0_05": "significant"})
    columns = ["intervention", "compartment", "neighborhoods", "significant", "enriched", "depleted"]
    result = pd.concat([ablation[columns], blockade[columns]], ignore_index=True)
    result["enriched_fraction"] = result.enriched / result.neighborhoods
    result["depleted_fraction"] = result.depleted / result.neighborhoods
    result["significant_fraction"] = result.significant / result.neighborhoods
    return result


def short_state(value):
    replacements = {
        " epithelial cell": "",
        " fibroblast": " fibro.",
        " neutrophil": " neut.",
        " monocyte/macrophage": " mono./macro.",
        " dendritic cell": " DC",
        "Conventional ": "Conv. ",
        "Inflammatory ": "Inflam. ",
        "transitional": "trans.",
        "responsive": "resp.",
        "secretory": "secr.",
    }
    for old, new in replacements.items():
        value = value.replace(old, new)
    return value


def draw_da_figure(primary, comparison, milo):
    fig = plt.figure(figsize=(183 / 25.4, 168 / 25.4), facecolor="white")
    fig.text(.025, .975, "FAP ablation versus α5β1 blockade: differential abundance",
             fontsize=9.0, fontweight="bold", color=INK, va="top")
    fig.text(.025, .945,
             "Primary inference uses biological mice",
             fontsize=5.4, color=MUTED, va="top")

    # a: descriptive neighborhood overview from separately constructed MiloR graphs.
    panel_label(fig, .025, .902, "a", "Neighborhood-level remodeling")
    ax = fig.add_axes([.145, .695, .315, .17])
    y_positions, y_labels = [], []
    cursor = 0
    for compartment in ["Epithelial", "Immune", "Stromal"]:
        for intervention in INTERVENTION_ORDER:
            row = milo.loc[milo.compartment.eq(compartment) & milo.intervention.eq(intervention)].iloc[0]
            y_positions.append(cursor)
            y_labels.append(f"{compartment} · {INTERVENTION_LABELS[intervention]}")
            ax.barh(cursor, -100 * row.depleted_fraction, color="#4C78A8", height=.62)
            ax.barh(cursor, 100 * row.enriched_fraction, color="#D85C5C", height=.62)
            cursor += 1
        cursor += .35
    ax.axvline(0, color="#687078", lw=.55)
    ax.set_yticks(y_positions, y_labels)
    ax.set_xlim(-100, 100)
    ax.set_xticks([-100, -50, 0, 50, 100], ["100", "50", "0", "50", "100"])
    ax.set_xlabel("Significant neighborhoods (%)", fontsize=5.2, labelpad=2)
    clean_axis(ax, grid="x")
    ax.tick_params(axis="y", length=0, labelsize=4.7)
    ax.spines["left"].set_visible(False)
    fig.legend(handles=[Line2D([0], [0], color="#4C78A8", lw=5, label="Lower after intervention"),
                        Line2D([0], [0], color="#D85C5C", lw=5, label="Higher after intervention")],
               frameon=False, fontsize=4.7, ncol=2, loc="center",
               bbox_to_anchor=(.30, .875))
    fig.text(.145, .655, "Descriptive only: neighborhood graphs and model eligibility differ by series.",
             fontsize=4.4, color=MUTED)

    # b: biological-mouse compartment recovery effects.
    panel_label(fig, .52, .902, "b", "Compartment recovery")
    ax = fig.add_axes([.61, .695, .35, .17])
    recovery = primary.loc[primary.feature_scope.eq("Compartment recovery")]
    y = np.arange(3)[::-1]
    for offset, intervention, marker in [(.12, "FAP ablation", "o"), (-.12, "alpha5beta1 blockade", "^")]:
        q = recovery.loc[recovery.intervention.eq(intervention)].set_index("compartment")
        order = ["Epithelial", "Immune", "Stromal"]
        values = q.loc[order, "effect_percentage_points"].to_numpy()
        lows = q.loc[order, "ci_low"].to_numpy()
        highs = q.loc[order, "ci_high"].to_numpy()
        ax.errorbar(values, y + offset, xerr=[values - lows, highs - values], fmt=marker,
                    color=INTERVENTION_COLORS[intervention], markersize=4.1, elinewidth=1,
                    capsize=1.7, label=INTERVENTION_LABELS[intervention], zorder=3)
    ax.axvline(0, color="#777E85", lw=.6)
    ax.set_yticks(y, ["Epithelial", "Immune", "Stromal"])
    ax.set_xlabel("Change after intervention (percentage points)", fontsize=5.2, labelpad=2)
    clean_axis(ax, grid="x")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.legend(frameon=False, fontsize=4.8, loc="lower right")

    # c: top unbiased state effects (top eight by maximum absolute effect per compartment).
    panel_label(fig, .025, .625, "c", "Level 3 state remodeling")
    states = comparison.loc[comparison.feature_scope.eq("Level 3 fraction within compartment")].copy()
    selected = (states.sort_values("max_abs_effect", ascending=False)
                .groupby("compartment", sort=False).head(8))
    selected["compartment_order"] = selected.compartment.map({"Epithelial": 0, "Immune": 1, "Stromal": 2})
    selected = selected.sort_values(["compartment_order", "effect_percentage_points_ablation"],
                                    ascending=[True, False]).reset_index(drop=True)
    matrix = selected[["effect_percentage_points_ablation", "effect_percentage_points_blockade"]].to_numpy()
    limit = max(10, np.nanmax(np.abs(matrix)))
    ax = fig.add_axes([.16, .14, .31, .43])
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
    ax.set_xticks([0, 1], ["FAP ablation", "α5β1 blockade"])
    ax.set_yticks(np.arange(len(selected)), [short_state(value) for value in selected.feature])
    ax.tick_params(length=0, labelsize=4.45, pad=1.5,
                   top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.spines[:].set_visible(False)
    for i, row in selected.iterrows():
        for j, suffix in enumerate(["ablation", "blockade"]):
            q_value = row[f"fdr_global_scope_{suffix}"]
            mark = "••" if q_value < .01 else ("•" if q_value < .05 else "")
            if mark:
                ax.text(j, i, mark, ha="center", va="center", fontsize=5.4, fontweight="bold",
                        color="white" if abs(matrix[i, j]) > limit * .35 else INK)
    last = None
    for i, compartment in enumerate(selected.compartment):
        if last is not None and compartment != last:
            ax.axhline(i - .5, color="white", lw=2.2)
        last = compartment
    cax = fig.add_axes([.19, .070, .25, .012])
    cb = fig.colorbar(image, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=4.5, length=1.5, pad=1)
    cb.set_label("Change within compartment (percentage points)", fontsize=4.8, labelpad=1.5)
    fig.text(.16, .112, "Dots: global-scope FDR <0.05 (two dots, <0.01)", fontsize=4.3, color=MUTED)

    # d: cross-intervention concordance of all shared Level 3 effects.
    panel_label(fig, .52, .625, "d", "Cross-intervention effect concordance")
    ax = fig.add_axes([.59, .13, .37, .42])
    for compartment, q in states.groupby("compartment"):
        ax.scatter(q.effect_percentage_points_blockade, q.effect_percentage_points_ablation,
                   s=15, color=COMPARTMENT_COLORS[compartment], alpha=.78,
                   edgecolor="white", linewidth=.35, label=compartment, zorder=3)
    bound = max(10, np.nanmax(np.abs(states[["effect_percentage_points_ablation",
                                            "effect_percentage_points_blockade"]].to_numpy()))) * 1.10
    ax.plot([-bound, bound], [-bound, bound], color="#A3AAB0", lw=.65, ls="--", zorder=1)
    ax.axhline(0, color="#A3AAB0", lw=.45)
    ax.axvline(0, color="#A3AAB0", lw=.45)
    ax.set_xlim(-bound, bound)
    ax.set_ylim(-bound, bound)
    ax.set_xlabel("α5β1-blockade effect (percentage points)", fontsize=5.2)
    ax.set_ylabel("FAP-ablation effect (percentage points)", fontsize=5.2)
    clean_axis(ax, grid="both")
    rho, p_value = spearmanr(states.effect_percentage_points_blockade,
                             states.effect_percentage_points_ablation)
    ax.text(.03, .97, f"Spearman ρ={rho:.2f}; P={p_value:.2g}\nn={len(states)} shared states",
            transform=ax.transAxes, ha="left", va="top", fontsize=4.8, color=MUTED)
    divergent = (states.assign(divergence=(states.effect_percentage_points_ablation -
                                           states.effect_percentage_points_blockade).abs())
                 .sort_values("divergence", ascending=False)
                 .groupby("compartment", sort=False).head(2))
    for side, q in divergent.groupby(divergent.effect_percentage_points_blockade.ge(0)):
        q = q.sort_values("effect_percentage_points_ablation").copy()
        desired = q.effect_percentage_points_ablation.to_numpy(float)
        minimum_gap = bound * .075
        for i in range(1, len(desired)):
            desired[i] = max(desired[i], desired[i - 1] + minimum_gap)
        desired -= np.mean(desired - q.effect_percentage_points_ablation.to_numpy(float))
        text_x = bound * (.50 if side else -.50)
        for y_text, row in zip(desired, q.itertuples(index=False)):
            ax.annotate(short_state(row.feature),
                        (row.effect_percentage_points_blockade, row.effect_percentage_points_ablation),
                        xytext=(text_x, y_text), textcoords="data",
                        fontsize=3.85, color=INK, va="center",
                        ha="left" if side else "right",
                        arrowprops={"arrowstyle": "-", "color": "#9BA2A8", "lw": .4})
    ax.legend(frameon=False, fontsize=4.7, loc="lower right")

    fig.text(.025, .018,
             "Bootstrap intervals resample biological mice equally. Direct head-to-head interpretation is limited because intervention and source batch are aligned.",
             fontsize=4.4, color=MUTED)
    target_pdf = FIGURES / "Figure_1_Differential_Abundance_Comparison.pdf"
    target_png = FIGURES / "Figure_1_Differential_Abundance_Comparison.png"
    target_svg = FIGURES / "Figure_1_Differential_Abundance_Comparison.svg"
    fig.savefig(target_pdf, dpi=450, metadata={"Title": "FAP ablation versus alpha5beta1 blockade differential abundance"})
    fig.savefig(target_png, dpi=450)
    fig.savefig(target_svg)
    plt.close(fig)
    return {"pdf": target_pdf, "png": target_png, "svg": target_svg,
            "state_rho": float(rho), "state_p": float(p_value), "states": int(len(states))}


def prepare_multinichenet():
    keys = ["direction", "ligand", "receptor"]
    fap = pd.read_csv(SOURCE / "FAP_ablation_MultiNicheNet_nonambient.csv")
    blockade = pd.read_csv(SOURCE / "a5B1_blockade_MultiNicheNet_nonambient.csv")

    def take(data, comparison, name):
        q = data.loc[data.comparison.eq(comparison), keys + ["score_delta"]].copy()
        assert not q.duplicated(keys).any()
        return q.rename(columns={"score_delta": name})

    merged = take(fap, "PBS_DSS_vs_PBS", "disease_score_fap_series")
    merged = merged.merge(take(fap, "PBS_DSS_vs_GCV_DSS", "attenuation_fap_ablation"), on=keys)
    merged = merged.merge(take(blockade, "DSS_vs_Control", "disease_score_blockade_series"), on=keys)
    merged = merged.merge(take(blockade, "DSS_vs_Blockade", "attenuation_a5b1_blockade"), on=keys)
    merged["disease_induced_both_series"] = ((merged.disease_score_fap_series > 0) &
                                               (merged.disease_score_blockade_series > 0))
    conditions = [
        (merged.attenuation_fap_ablation > 0) & (merged.attenuation_a5b1_blockade > 0),
        (merged.attenuation_fap_ablation > 0) & (merged.attenuation_a5b1_blockade <= 0),
        (merged.attenuation_fap_ablation <= 0) & (merged.attenuation_a5b1_blockade > 0),
    ]
    merged["attenuation_class"] = np.select(
        conditions,
        ["Attenuated by both", "FAP ablation only", "α5β1 blockade only"],
        default="Not attenuated by either",
    )
    merged["consensus_strength"] = np.where(
        merged.disease_induced_both_series & merged.attenuation_class.eq("Attenuated by both"),
        np.minimum(merged.attenuation_fap_ablation, merged.attenuation_a5b1_blockade), np.nan)
    merged["interaction"] = merged.ligand + " → " + merged.receptor
    merged["attenuation_difference_ablation_minus_blockade"] = (
        merged.attenuation_fap_ablation - merged.attenuation_a5b1_blockade)
    focused_index = {route: index for index, route in enumerate(FOCUSED_ROUTES)}
    merged["prespecified_focus"] = [tuple(row) in focused_index for row in merged[keys].to_numpy()]
    merged["focus_order"] = [focused_index.get(tuple(row), np.nan) for row in merged[keys].to_numpy()]
    return merged


def draw_mnn_figure(data):
    fig = plt.figure(figsize=(183 / 25.4, 155 / 25.4), facecolor="white")
    fig.text(.025, .975, "FAP ablation versus α5β1 blockade: reciprocal signaling",
             fontsize=9.0, fontweight="bold", color=INK, va="top")
    fig.text(.025, .945,
             "Comparative re-analysis of matched nonambient neutrophil–fibroblast MultiNicheNet 2.1.0 outputs",
             fontsize=5.4, color=MUTED, va="top")

    core = data.loc[data.disease_induced_both_series].copy()
    panel_label(fig, .025, .90, "a", "Interaction-level concordance")
    ax = fig.add_axes([.09, .52, .43, .33])
    direction_colors = {"Fibroblast -> Neutrophil": EPITHELIAL,
                        "Neutrophil -> Fibroblast": IMMUNE}
    for direction, q in data.groupby("direction"):
        ax.scatter(q.attenuation_a5b1_blockade, q.attenuation_fap_ablation,
                   s=9.5, alpha=.42, color=direction_colors[direction], linewidth=0,
                   label=direction.replace(" -> ", " → "), rasterized=True)
    ax.axhline(0, color="#8E969D", lw=.55)
    ax.axvline(0, color="#8E969D", lw=.55)
    rho, p_value = spearmanr(data.attenuation_a5b1_blockade, data.attenuation_fap_ablation)
    ax.text(.03, .97, f"All common routes\nSpearman ρ={rho:.2f}; P={p_value:.2g}; n={len(data)}",
            transform=ax.transAxes, va="top", fontsize=4.8, color=MUTED)
    clean_axis(ax, grid="both")
    ax.set_xlabel("DSS − α5β1 blockade priority-score difference", fontsize=5.1)
    ax.set_ylabel("DSS − FAP ablation priority-score difference", fontsize=5.1)
    ax.legend(frameon=False, fontsize=4.5, loc="lower right")

    panel_label(fig, .56, .90, "b", "Disease-induced route classification")
    ax = fig.add_axes([.63, .59, .33, .23])
    categories = ["Attenuated by both", "FAP ablation only", "α5β1 blockade only",
                  "Not attenuated by either"]
    category_colors = ["#3E8F72", ABLATION, BLOCKADE, "#B9BEC3"]
    counts = (core.groupby(["direction", "attenuation_class"]).size().unstack(fill_value=0)
              .reindex(columns=categories, fill_value=0))
    left = np.zeros(len(counts))
    y = np.arange(len(counts))[::-1]
    for category, color in zip(categories, category_colors):
        values = counts[category].to_numpy()
        ax.barh(y, values, left=left, height=.52, color=color, label=category)
        for yy, xx, ll in zip(y, values, left):
            if xx >= 20:
                ax.text(ll + xx / 2, yy, str(int(xx)), ha="center", va="center",
                        fontsize=4.5, color="white" if color != "#B9BEC3" else INK)
        left += values
    ax.set_yticks(y, [value.replace(" -> ", " → ") for value in counts.index])
    ax.set_xlabel("Routes induced by DSS in both source series (count)", fontsize=5.0)
    clean_axis(ax, grid="x")
    ax.tick_params(axis="y", length=0, labelsize=4.8)
    ax.spines["left"].set_visible(False)
    fig.legend(handles=[Line2D([0], [0], color=color, lw=5, label=category)
                        for category, color in zip(categories, category_colors)],
               frameon=False, fontsize=4.25, ncol=2, loc="center",
               bbox_to_anchor=(.79, .858))

    panel_label(fig, .025, .455, "c", "Prespecified disease-anchored reciprocal routes")
    selected = data.loc[data.prespecified_focus].sort_values("focus_order")
    directions = ["Fibroblast -> Neutrophil", "Neutrophil -> Fibroblast"]
    for block, direction in enumerate(directions):
        q = selected.loc[selected.direction.eq(direction)].sort_values("focus_order", ascending=False)
        ax = fig.add_axes([.15 + block * .49, .105, .31, .285])
        y = np.arange(len(q))
        for yy, row in enumerate(q.itertuples(index=False)):
            values = [row.attenuation_fap_ablation, row.attenuation_a5b1_blockade]
            ax.plot(values, [yy, yy], color="#BCC2C7", lw=.9, zorder=1)
            ax.scatter(values[0], yy, s=17, color=ABLATION, marker="o", zorder=3)
            ax.scatter(values[1], yy, s=21, color=BLOCKADE, marker="^", zorder=4)
        ax.set_yticks(y, q.interaction)
        ax.set_xlim(left=0)
        clean_axis(ax, grid="x")
        ax.tick_params(axis="y", length=0, labelsize=4.65)
        ax.spines["left"].set_visible(False)
        for tick in ax.get_yticklabels():
            tick.set_fontstyle("italic")
        ax.set_xlabel("DSS − intervention priority-score difference", fontsize=4.9)
        fig.text(.15 + block * .49, .405, direction.replace(" -> ", " → "),
                 fontsize=5.5, fontweight="bold", color=INK)
    fig.legend(handles=[Line2D([0], [0], marker="o", color=ABLATION, lw=0, markersize=4,
                               label="FAP ablation"),
                        Line2D([0], [0], marker="^", color=BLOCKADE, lw=0, markersize=4,
                               label="α5β1 blockade")],
               frameon=False, ncol=2, fontsize=4.8, loc="lower center", bbox_to_anchor=(.5, .035))
    fig.text(.025, .012,
             "Priority-score differences are ranking quantities, not p-values or signaling-flux estimates. Ablation: 4 neutrophils across 3 eligible units; blockade: 3 eligible mice, two with ≤11 fibroblasts.",
             fontsize=4.15, color=MUTED)

    target_pdf = FIGURES / "Figure_2_MultiNicheNet_Interaction_Comparison.pdf"
    target_png = FIGURES / "Figure_2_MultiNicheNet_Interaction_Comparison.png"
    target_svg = FIGURES / "Figure_2_MultiNicheNet_Interaction_Comparison.svg"
    fig.savefig(target_pdf, dpi=450, metadata={"Title": "FAP ablation versus alpha5beta1 blockade MultiNicheNet comparison"})
    fig.savefig(target_png, dpi=450)
    fig.savefig(target_svg)
    plt.close(fig)
    return {"pdf": target_pdf, "png": target_png, "svg": target_svg,
            "rho": float(rho), "p": float(p_value), "common_routes": int(len(data)),
            "disease_core": int(len(core)),
            "attenuated_by_both": int(core.attenuation_class.eq("Attenuated by both").sum())}


def build_report(primary, sensitivity, comparison, mnn, da_meta, mnn_meta):
    recovery = primary.loc[primary.feature_scope.eq("Compartment recovery")]
    recovery_lines = []
    for intervention in INTERVENTION_ORDER:
        q = recovery.loc[recovery.intervention.eq(intervention)].set_index("compartment")
        values = "; ".join(
            f"{compartment} {q.loc[compartment, 'effect_percentage_points']:+.1f} pp "
            f"(95% CI {q.loc[compartment, 'ci_low']:+.1f} to {q.loc[compartment, 'ci_high']:+.1f})"
            for compartment in ["Epithelial", "Immune", "Stromal"]
        )
        recovery_lines.append(f"- **{INTERVENTION_LABELS[intervention]}:** {values}.")

    states = comparison.loc[comparison.feature_scope.eq("Level 3 fraction within compartment")]
    divergent = states.assign(divergence=(states.effect_percentage_points_ablation -
                                          states.effect_percentage_points_blockade).abs()).nlargest(8, "divergence")
    divergent_lines = [
        f"- {row.feature} ({row.compartment}): ablation {row.effect_percentage_points_ablation:+.1f} pp; "
        f"blockade {row.effect_percentage_points_blockade:+.1f} pp."
        for row in divergent.itertuples(index=False)
    ]

    core = mnn.loc[mnn.disease_induced_both_series]
    class_counts = core.attenuation_class.value_counts()
    top = (mnn.loc[mnn.prespecified_focus & mnn.consensus_strength.notna()]
           .nlargest(12, "consensus_strength"))
    top_lines = [
        f"- {row.direction.replace(' -> ', ' → ')}: {row.interaction} "
        f"(ablation {row.attenuation_fap_ablation:.3f}; blockade {row.attenuation_a5b1_blockade:.3f})."
        for row in top.itertuples(index=False)
    ]
    sensitivity_recovery = sensitivity.loc[sensitivity.feature_scope.eq("Compartment recovery")]
    sensitivity_delta = primary.merge(
        sensitivity,
        on=["intervention", "feature_scope", "compartment", "feature"],
        suffixes=("_primary", "_sensitivity"),
    )
    max_sensitivity_change = (sensitivity_delta.effect_percentage_points_primary -
                              sensitivity_delta.effect_percentage_points_sensitivity).abs().max()

    report = f"""# FAP ablation versus α5β1 blockade comparative analysis

## Scope and design

This noncanonical analysis compares two interventions during DSS colitis using existing finalized datasets. The six DSS reference mice (D1–D6) are shared. The primary differential-abundance analysis gives each included biological mouse equal weight. It includes 6 DSS reference mice, 5 FAP-ablation mice and 12 α5β1-blockade mice. A prespecified sensitivity analysis excludes samples with fewer than 200 recovered cells, leaving 4 FAP-ablation mice; the largest change in any reported effect estimate was {max_sensitivity_change:.1f} percentage points.

Treatment and source batch are aligned (DSS=M1, ablation=M3, blockade=I1). Accordingly, the within-intervention estimates are interpretable as cohort-specific treatment-associated effects, while the direct ablation-versus-blockade contrast is descriptive and not an unbiased causal head-to-head comparison.

## Differential abundance methods

- Compartment recovery is expressed as a fraction of all recovered epithelial, immune and stromal cells.
- Level 3 abundance is expressed within its parent compartment.
- Effects are differences in equally weighted biological-mouse mean proportions, reported in percentage points.
- Confidence intervals are percentile intervals from 10,000 biological-mouse bootstrap resamples.
- Two-sided Mann–Whitney tests compare each intervention with the shared DSS reference group; Benjamini–Hochberg correction is applied across all features within each intervention and feature scope.
- The MiloR panel summarizes existing neighborhood models. Each intervention used a separately constructed graph and its source-model eligibility rules; neighborhood fractions are descriptive and are not directly pooled.

## Differential abundance results

### Compartment recovery

{chr(10).join(recovery_lines)}

Across {da_meta['states']} shared Level 3 states, the intervention-effect vectors had Spearman ρ={da_meta['state_rho']:.2f} (P={da_meta['state_p']:.2g}). This correlation measures pattern similarity, not equivalence of intervention magnitude.

### Largest differences between intervention-associated state effects

{chr(10).join(divergent_lines)}

Positive values denote a larger within-compartment fraction after intervention. A positive relative fraction does not necessarily imply greater absolute cell recovery when the parent compartment is depleted.

## MultiNicheNet comparative methods

Matched nonambient fibroblast–neutrophil outputs from separate MultiNicheNet 2.1.0 runs were joined by direction, ligand and receptor. Four score differences were retained: DSS versus control and DSS versus intervention in each source series. A disease-anchored route required a positive DSS-versus-control score difference in both series. A concordantly attenuated route additionally required positive DSS-versus-ablation and DSS-versus-blockade score differences. `consensus_strength` is the smaller of the two intervention-associated score differences, so highly ranked routes must be supported by both interventions.

MultiNicheNet priority-score differences integrate expression, specificity, prevalence, differential expression and predicted ligand activity. They are ranking quantities, not p-values, direct ligand–receptor binding measurements or estimates of signaling flux.

## MultiNicheNet results

{mnn_meta['common_routes']} ligand–receptor routes were evaluable in all four matched comparisons. {mnn_meta['disease_core']} were DSS-induced in both source series; among these, {mnn_meta['attenuated_by_both']} were lower after both interventions. The full interaction-level attenuation vectors had Spearman ρ={mnn_meta['rho']:.2f} (P={mnn_meta['p']:.2g}).

Disease-anchored classifications: {', '.join(f'{key}: {int(class_counts.get(key, 0))}' for key in ['Attenuated by both', 'FAP ablation only', 'α5β1 blockade only', 'Not attenuated by either'])}.

### Highest prespecified consensus routes

{chr(10).join(top_lines)}

## Critical limitations

- FAP-ablation MultiNicheNet retained only four neutrophils across three source-model units. This is an extreme-depletion, hypothesis-generating comparison rather than reliable mouse-level signaling inference.
- The blockade MultiNicheNet model retained only M16, M21 and M23; M21 and M23 had 2 and 11 fibroblasts, respectively.
- The two treatment cohorts are batch-aligned, and the interaction models were fit separately. Apparent differences between interventions may reflect cohort, recovery, batch or normalization differences.
- FAP ablation changes sender/receiver occupancy profoundly. Reduced MultiNicheNet priority may therefore reflect depletion as well as within-cell-state regulation.
- No FACS re-analysis or new experiment was included.

## Deliverables

- `figures/Figure_1_Differential_Abundance_Comparison.*`
- `figures/Figure_2_MultiNicheNet_Interaction_Comparison.*`
- `Comparative_Analysis_Figures_Combined.pdf`
- Complete and display-level tables in `tables/`
- Copied, immutable analysis inputs in `source_data/`
"""
    (HERE / "ANALYSIS_REPORT.md").write_text(report, encoding="utf-8")


def combine_pdfs(paths):
    from pypdf import PdfReader, PdfWriter
    writer = PdfWriter()
    for path in paths:
        for page in PdfReader(str(path)).pages:
            writer.add_page(page)
    target = HERE / "Comparative_Analysis_Figures_Combined.pdf"
    with target.open("wb") as handle:
        writer.write(handle)
    return target


def make_contact_sheet(paths):
    images = [Image.open(path).convert("RGB") for path in paths]
    width = max(image.width for image in images)
    resized = []
    for image in images:
        scale = width / image.width
        resized.append(image.resize((width, round(image.height * scale)), Image.Resampling.LANCZOS))
    gutter = 50
    canvas = Image.new("RGB", (width, sum(image.height for image in resized) + gutter * 3), "white")
    y = gutter
    for image in resized:
        canvas.paste(image, (0, y))
        y += image.height + gutter
    target = HERE / "Comparative_Analysis_Figures_contact_sheet.png"
    canvas.save(target, dpi=(200, 200), optimize=True)
    return target


def main():
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    before = canonical_hashes()
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": 5.3, "axes.linewidth": .55,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white",
    })

    raw = pd.read_csv(SOURCE / "sample_level_remodeling_fractions.csv")
    biological = pd.read_csv(WORK / "config/biological_mouse_units.csv")
    raw = raw.loc[raw.contrast_short.isin(INTERVENTION_ORDER) & raw.sample_uid.isin(biological.sample_uid)].copy()
    primary = calculate_da(raw, minimum_total_cells=None, seed=20260904)
    sensitivity = calculate_da(raw, minimum_total_cells=200, seed=20260906)
    comparison = comparative_da_table(primary, raw)
    milo = milo_overview()
    primary.to_csv(TABLES / "differential_abundance_mouse_level_primary.csv", index=False)
    sensitivity.to_csv(TABLES / "differential_abundance_sensitivity_min200_cells.csv", index=False)
    comparison.to_csv(TABLES / "differential_abundance_intervention_comparison.csv", index=False)
    milo.to_csv(TABLES / "miloR_neighborhood_remodeling_overview.csv", index=False)
    da_meta = draw_da_figure(primary, comparison, milo)

    mnn = prepare_multinichenet()
    mnn.to_csv(TABLES / "MultiNicheNet_all_matched_interactions.csv", index=False)
    mnn.loc[mnn.disease_induced_both_series].to_csv(
        TABLES / "MultiNicheNet_disease_induced_core.csv", index=False)
    mnn.loc[mnn.consensus_strength.notna()].sort_values("consensus_strength", ascending=False).to_csv(
        TABLES / "MultiNicheNet_concordantly_attenuated_ranked.csv", index=False)
    mnn_meta = draw_mnn_figure(mnn)

    combined = combine_pdfs([da_meta["pdf"], mnn_meta["pdf"]])
    contact = make_contact_sheet([da_meta["png"], mnn_meta["png"]])
    build_report(primary, sensitivity, comparison, mnn, da_meta, mnn_meta)

    after = canonical_hashes()
    manifest = {
        "status": "NONCANONICAL COMPARATIVE ANALYSIS",
        "date": "2026-09-04",
        "canonical_outputs_modified": before != after,
        "primary_da": {key: value for key, value in da_meta.items() if key not in {"pdf", "png", "svg"}},
        "multinichenet": {key: value for key, value in mnn_meta.items() if key not in {"pdf", "png", "svg"}},
        "inputs": {path.name: {"sha256": sha256(path), "bytes": path.stat().st_size}
                   for path in sorted(SOURCE.glob("*"))},
        "outputs": [str(path.relative_to(HERE)) for path in sorted(TABLES.glob("*"))] +
                   [str(path.relative_to(HERE)) for path in sorted(FIGURES.glob("*"))] +
                   [combined.name, contact.name, "ANALYSIS_REPORT.md"],
    }
    (HERE / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    if before != after:
        raise RuntimeError("Canonical figure files changed during analysis")
    print(json.dumps({
        "built": True,
        "canonical_outputs_modified": False,
        "differential_abundance_states": da_meta["states"],
        "multinichenet_common_routes": mnn_meta["common_routes"],
        "multinichenet_disease_core": mnn_meta["disease_core"],
        "multinichenet_attenuated_by_both": mnn_meta["attenuated_by_both"],
        "combined_pdf": str(combined),
    }, indent=2))


if __name__ == "__main__":
    main()
