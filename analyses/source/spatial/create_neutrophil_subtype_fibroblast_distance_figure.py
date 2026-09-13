from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore")

ROOT = Path("giotto_codex_results")
OUT = ROOT / "fap_high_a5b1_high_rerun"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GROUPS = ("UC_Noninflamed", "UC_Inflamed")
GROUP_LABEL = {"UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}
GROUP_COLOR = {"UC_Noninflamed": "#D4AD35", "UC_Inflamed": "#D95F5F"}
TARGETS = ("FAP-high", "a5B1-high", "FAP-high/a5B1-high")
SUBTYPES = (
    "Neutrophil", "Neutrophil_CXCR4", "Neutrophil_MX1",
    "Neutrophil_OSM", "Neutrophil_PADI4",
)
SUBTYPE_LABEL = {
    "Neutrophil": "Unassigned",
    "Neutrophil_CXCR4": "CXCR4+",
    "Neutrophil_MX1": "MX1+",
    "Neutrophil_OSM": "OSM+",
    "Neutrophil_PADI4": "PADI4+",
}
N_RANDOM = 200
SEED = 20260807
FDR_ALPHA = 0.10


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    x = p[ok]
    order = np.argsort(x)
    ranked = x[order] * len(x) / np.arange(1, len(x) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1)
    out[ok] = adjusted
    return out


rng = np.random.default_rng(SEED)
usecols = [
    "cell_ID", "PatientID", "Diagnosis2", "cell_type", "x", "y",
    "FAPa_cell_clr", "a5B1_cell_clr",
]
d = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv", usecols=usecols)
d = d[d.Diagnosis2.isin(GROUPS)].copy()
d["is_fibroblast"] = d.cell_type.eq("Fibroblast")

thresholds = pd.read_csv(TAB / "00_patient_high_thresholds.csv")
d = d.merge(
    thresholds[["PatientID", "FAP_q75_clr", "a5B1_q75_clr"]],
    on="PatientID", how="left", validate="many_to_one",
)
d["FAP-high"] = d.is_fibroblast & d.FAPa_cell_clr.ge(d.FAP_q75_clr)
d["a5B1-high"] = d.is_fibroblast & d.a5B1_cell_clr.ge(d.a5B1_q75_clr)
d["FAP-high/a5B1-high"] = d["FAP-high"] & d["a5B1-high"]

rows = []
for (patient, diagnosis), z in d.groupby(["PatientID", "Diagnosis2"], sort=True):
    fibroblasts = z[z.is_fibroblast]
    fibro_xy = fibroblasts[["x", "y"]].to_numpy()
    subtype_xy = {
        subtype: z.loc[z.cell_type.eq(subtype), ["x", "y"]].to_numpy()
        for subtype in SUBTYPES
    }
    for target in TARGETS:
        target_xy = z.loc[z[target], ["x", "y"]].to_numpy()
        if len(target_xy) == 0 or len(fibro_xy) < len(target_xy):
            continue
        observed = {}
        for subtype, xy in subtype_xy.items():
            if len(xy):
                observed[subtype] = float(np.median(cKDTree(target_xy).query(xy)[0]))
        expected = {subtype: [] for subtype in observed}
        for _ in range(N_RANDOM):
            random_idx = rng.choice(len(fibro_xy), size=len(target_xy), replace=False)
            random_tree = cKDTree(fibro_xy[random_idx])
            for subtype, xy in subtype_xy.items():
                if subtype in observed:
                    expected[subtype].append(float(np.median(random_tree.query(xy)[0])))
        for subtype, observed_median in observed.items():
            expected_values = np.asarray(expected[subtype])
            ratios = np.log2((observed_median + 1e-9) / (expected_values + 1e-9))
            rows.append({
                "PatientID": patient,
                "Diagnosis2": diagnosis,
                "target": target,
                "neutrophil_subtype": subtype,
                "n_subtype_cells": len(subtype_xy[subtype]),
                "n_target_fibroblasts": len(target_xy),
                "observed_median_distance_um": observed_median,
                "expected_median_random_distance_um": float(np.median(expected_values)),
                "median_log2_observed_random_distance": float(np.median(ratios)),
            })

patient_results = pd.DataFrame(rows)
patient_results.to_csv(TAB / "17_neutrophil_subtype_proximity_by_patient.csv", index=False)

stats_rows = []
for (target, subtype), q in patient_results.groupby(["target", "neutrophil_subtype"]):
    noninflamed = q.loc[q.Diagnosis2.eq("UC_Noninflamed"), "median_log2_observed_random_distance"].to_numpy()
    inflamed = q.loc[q.Diagnosis2.eq("UC_Inflamed"), "median_log2_observed_random_distance"].to_numpy()
    if len(noninflamed) and len(inflamed):
        test = mannwhitneyu(inflamed, noninflamed, alternative="two-sided")
        stats_rows.append({
            "target": target,
            "neutrophil_subtype": subtype,
            "analysis": "UC inflamed versus UC noninflamed",
            "group": "",
            "n_UC_noninflamed": len(noninflamed),
            "n_UC_inflamed": len(inflamed),
            "median_UC_noninflamed": np.median(noninflamed),
            "median_UC_inflamed": np.median(inflamed),
            "effect_inflamed_minus_noninflamed": np.median(inflamed) - np.median(noninflamed),
            "p_value": test.pvalue,
        })
    for group, values in (("UC_Noninflamed", noninflamed), ("UC_Inflamed", inflamed)):
        if len(values) >= 3 and np.any(values != 0):
            test = wilcoxon(values, alternative="two-sided", zero_method="wilcox")
            stats_rows.append({
                "target": target,
                "neutrophil_subtype": subtype,
                "analysis": "Effect versus zero",
                "group": group,
                "n_UC_noninflamed": len(noninflamed),
                "n_UC_inflamed": len(inflamed),
                "median_UC_noninflamed": np.median(noninflamed),
                "median_UC_inflamed": np.median(inflamed),
                "effect_inflamed_minus_noninflamed": np.median(values),
                "p_value": test.pvalue,
            })

stats = pd.DataFrame(stats_rows)
stats["fdr_within_analysis"] = stats.groupby("analysis").p_value.transform(bh_adjust)
stats.to_csv(TAB / "18_neutrophil_subtype_proximity_statistics.csv", index=False)

# Companion figure to panel D: patient-level values, equal patient weighting.
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})
order = [SUBTYPE_LABEL[s] for s in SUBTYPES]
hue_order = [GROUP_LABEL[g] for g in GROUPS]
palette = {GROUP_LABEL[g]: GROUP_COLOR[g] for g in GROUPS}
plot_data = patient_results.copy()
plot_data["Subtype"] = plot_data.neutrophil_subtype.map(SUBTYPE_LABEL)
plot_data["Group"] = plot_data.Diagnosis2.map(GROUP_LABEL)

between = stats[stats.analysis.eq("UC inflamed versus UC noninflamed")].copy()
global_min = plot_data.median_log2_observed_random_distance.min()
global_max = plot_data.median_log2_observed_random_distance.max()
span = max(global_max - global_min, 1.0)
annotation_y = global_max + 0.11 * span

fig, axes = plt.subplots(1, 3, figsize=(19, 6.8), sharey=True)
for i, (ax, target) in enumerate(zip(axes, TARGETS)):
    q = plot_data[plot_data.target.eq(target)]
    sns.boxplot(
        data=q, x="Subtype", y="median_log2_observed_random_distance",
        hue="Group", order=order, hue_order=hue_order, palette=palette,
        showfliers=False, width=.68, linewidth=1.2, ax=ax,
    )
    sns.stripplot(
        data=q, x="Subtype", y="median_log2_observed_random_distance",
        hue="Group", order=order, hue_order=hue_order, palette=palette,
        dodge=True, jitter=.10, size=4.2, alpha=.72, linewidth=.35,
        edgecolor="black", legend=False, ax=ax,
    )
    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ax.set_title(target)
    ax.set_xlabel("")
    ax.set_ylabel("log2 observed / random nearest distance" if i == 0 else "")
    ax.tick_params(axis="x", rotation=28, labelsize=9)
    ax.set_ylim(global_min - .08 * span, global_max + .24 * span)
    for xpos, subtype in enumerate(SUBTYPES):
        r = between[
            between.target.eq(target) & between.neutrophil_subtype.eq(subtype)
        ].iloc[0]
        ax.text(
            xpos, annotation_y,
            f"P={r.p_value:.3g}\nq={r.fdr_within_analysis:.3g}",
            ha="center", va="bottom", fontsize=7.5,
        )
    handles, labels = ax.get_legend_handles_labels()
    ax.get_legend().remove()

fig.tight_layout(rect=(0, 0, 1, .79))
fig.suptitle("Neutrophil subtype proximity to high-state fibroblasts in UC", fontweight="bold", y=.985)
fig.text(
    .5, .925,
    "Each point is one patient; positive values indicate greater distance (lower proximity). P and q compare UC inflamed with UC noninflamed.",
    ha="center", fontsize=11,
)
fig.legend(handles[:2], labels[:2], title="", frameon=False, fontsize=9,
           loc="upper center", bbox_to_anchor=(.5, .885), ncol=2)
fig.savefig(FIG / "Figure_Neutrophil_Subtype_Fibroblast_Distance.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(FIG / "Figure_Neutrophil_Subtype_Fibroblast_Distance.pdf", bbox_inches="tight", facecolor="white")
plt.close(fig)

sig = between[between.fdr_within_analysis.lt(FDR_ALPHA)].sort_values(
    ["target", "fdr_within_analysis", "p_value"]
)
if len(sig):
    sig_lines = []
    for r in sig.itertuples():
        direction = "greater" if r.effect_inflamed_minus_noninflamed < 0 else "lower"
        sig_lines.append(
            f"- {SUBTYPE_LABEL[r.neutrophil_subtype]} neutrophils to {r.target} fibroblasts: "
            f"{direction} proximity in UC inflamed (inflamed-minus-noninflamed normalized distance "
            f"{r.effect_inflamed_minus_noninflamed:+.3f}, raw P={r.p_value:.3g}, "
            f"FDR={r.fdr_within_analysis:.3g})."
        )
else:
    sig_lines = [f"- No subtype-specific UC-state differences met FDR < {FDR_ALPHA:.2f}."]

n_note = patient_results.groupby(["neutrophil_subtype", "Diagnosis2"]).PatientID.nunique().unstack(fill_value=0)
n_lines = [
    f"- {SUBTYPE_LABEL[s]}: {int(n_note.loc[s, 'UC_Noninflamed'])} UC noninflamed and "
    f"{int(n_note.loc[s, 'UC_Inflamed'])} UC inflamed patients."
    for s in SUBTYPES
]

summary = f"""# Neutrophil subtype-to-fibroblast proximity

The analysis includes UC noninflamed and UC inflamed samples only. For each patient, neutrophil subtype, and high-state fibroblast target, the median nearest-cell distance was divided by the median distance to 200 equal-sized random fibroblast subsets from the same patient. Patient-level values were compared so samples contributed equally regardless of cell count. Positive log2 values indicate lower proximity; negative inflamed-minus-noninflamed effects indicate greater proximity in inflamed UC. Raw two-sided P values and Benjamini-Hochberg FDR values are reported, with FDR < {FDR_ALPHA:.2f} as the discovery threshold.

## UC-state differences meeting FDR < {FDR_ALPHA:.2f}

{chr(10).join(sig_lines)}

## Patient coverage

{chr(10).join(n_lines)}

No patient with an observed subtype was excluded. P21 had no PADI4 subtype cells and therefore has no PADI4-specific distance estimate; all other patient-subtype combinations are represented.
"""
(OUT / "NEUTROPHIL_SUBTYPE_PROXIMITY_SUMMARY.md").write_text(summary, encoding="utf-8")

print("Created neutrophil subtype-to-fibroblast distance figure and statistics.")
