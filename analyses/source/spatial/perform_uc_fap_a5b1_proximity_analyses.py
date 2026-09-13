from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore")

ROOT = Path("giotto_codex_results")
SOURCE = ROOT / "fap_a5b1_neutrophil_reprogramming"
OUT = ROOT / "uc_fap_a5b1_proximity_analyses"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GROUPS = ("UC_Noninflamed", "UC_Inflamed")
GROUP_LABEL = {"UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}
GROUP_COLOR = {"UC_Noninflamed": "#D4AD35", "UC_Inflamed": "#D95F5F"}
TARGETS = {
    "FAP+ fibroblasts": ("FAP+/a5B1+", "FAP+/a5B1-"),
    "a5B1+ fibroblasts": ("FAP+/a5B1+", "FAP-/a5B1+"),
    "FAP+/a5B1+ fibroblasts": ("FAP+/a5B1+",),
}
RADII = (25, 50, 100)
N_NEUTROPHILS = 30
N_BOOT = 200
SEED = 20260805


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


def rank_biserial(x, y):
    u = mannwhitneyu(x, y, alternative="two-sided").statistic
    return 2 * u / (len(x) * len(y)) - 1


def uc_tests(frame, value, endpoint, zero_tests=False):
    a = frame.loc[frame.Diagnosis2.eq("UC_Noninflamed"), value].dropna().to_numpy()
    b = frame.loc[frame.Diagnosis2.eq("UC_Inflamed"), value].dropna().to_numpy()
    rows = []
    if len(a) and len(b):
        test = mannwhitneyu(b, a, alternative="two-sided")
        rows.append({
            "endpoint": endpoint, "analysis": "UC inflamed versus UC noninflamed",
            "group": "", "n_UC_noninflamed": len(a), "n_UC_inflamed": len(b),
            "median_UC_noninflamed": np.median(a), "median_UC_inflamed": np.median(b),
            "effect_inflamed_minus_noninflamed": np.median(b) - np.median(a),
            "rank_biserial": rank_biserial(b, a), "p_value": test.pvalue,
        })
    if zero_tests:
        for group, x in (("UC_Noninflamed", a), ("UC_Inflamed", b)):
            if len(x) >= 3 and np.any(x != 0):
                test = wilcoxon(x, alternative="two-sided", zero_method="wilcox")
                rows.append({
                    "endpoint": endpoint, "analysis": "Normalized proximity versus random expectation",
                    "group": group, "n_UC_noninflamed": len(a), "n_UC_inflamed": len(b),
                    "median_UC_noninflamed": np.median(a), "median_UC_inflamed": np.median(b),
                    "effect_inflamed_minus_noninflamed": np.median(x),
                    "rank_biserial": np.nan, "p_value": test.pvalue,
                })
    return rows


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


cells = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv", low_memory=False)
thresholds = pd.read_csv(SOURCE / "tables" / "00_fibroblast_state_thresholds.csv").set_index("marker")
fap_threshold = float(thresholds.loc["FAP", "primary_threshold_clr"])
a5_threshold = float(thresholds.loc["a5B1", "primary_threshold_clr"])

cells["is_fibroblast"] = cells.cell_type.eq("Fibroblast")
cells["is_neutrophil"] = cells.cell_type.str.startswith("Neutrophil", na=False)
fap = cells.FAPa_cell_clr.ge(fap_threshold)
a5 = cells.a5B1_cell_clr.ge(a5_threshold)
cells["fibroblast_state"] = "Non-fibroblast"
cells.loc[cells.is_fibroblast & fap & a5, "fibroblast_state"] = "FAP+/a5B1+"
cells.loc[cells.is_fibroblast & fap & ~a5, "fibroblast_state"] = "FAP+/a5B1-"
cells.loc[cells.is_fibroblast & ~fap & a5, "fibroblast_state"] = "FAP-/a5B1+"
cells.loc[cells.is_fibroblast & ~fap & ~a5, "fibroblast_state"] = "FAP-/a5B1-"
cells = cells[cells.Diagnosis2.isin(GROUPS)].copy()

# Exact cell-level nearest distances in both directions.
neut_distance_rows = []
fib_distance_rows = []
patient_rows = []
for pid, z in cells.groupby("PatientID", sort=True):
    fib = z[z.is_fibroblast]
    neut = z[z.is_neutrophil]
    neut_xy = neut[["x", "y"]].to_numpy()
    if fib.empty or neut.empty:
        continue
    neut_tree = cKDTree(neut_xy)
    for target, states in TARGETS.items():
        target_fib = fib[fib.fibroblast_state.isin(states)]
        if target_fib.empty:
            continue
        target_xy = target_fib[["x", "y"]].to_numpy()
        n_to_f = cKDTree(target_xy).query(neut_xy)[0]
        f_to_n = neut_tree.query(target_xy)[0]
        neut_distance_rows.extend({
            "cell_ID": cid, "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
            "target": target, "neutrophil_to_target_fibroblast_um": dist,
        } for cid, dist in zip(neut.cell_ID, n_to_f))
        fib_distance_rows.extend({
            "cell_ID": cid, "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
            "target": target, "target_fibroblast_to_neutrophil_um": dist,
        } for cid, dist in zip(target_fib.cell_ID, f_to_n))
        row = {
            "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0], "target": target,
            "n_neutrophils": len(neut), "n_target_fibroblasts": len(target_fib),
            "median_neutrophil_to_target_um": np.median(n_to_f),
            "q25_neutrophil_to_target_um": np.quantile(n_to_f, .25),
            "q75_neutrophil_to_target_um": np.quantile(n_to_f, .75),
            "median_target_to_neutrophil_um": np.median(f_to_n),
        }
        for radius in RADII:
            row[f"fraction_neutrophils_within_{radius}um"] = np.mean(n_to_f <= radius)
            row[f"fraction_target_fibroblasts_within_{radius}um_of_neutrophil"] = np.mean(f_to_n <= radius)
        patient_rows.append(row)

neut_distances = pd.DataFrame(neut_distance_rows)
fib_distances = pd.DataFrame(fib_distance_rows)
patient_metrics = pd.DataFrame(patient_rows)
neut_distances.to_csv(TAB / "01_neutrophil_to_target_distances_by_cell.csv", index=False)
fib_distances.to_csv(TAB / "02_target_fibroblast_to_neutrophil_distances_by_cell.csv", index=False)
patient_metrics.to_csv(TAB / "03_raw_proximity_metrics_by_patient.csv", index=False)

raw_stats_rows = []
for target, q in patient_metrics.groupby("target"):
    raw_stats_rows.extend(uc_tests(q, "median_neutrophil_to_target_um", f"{target}; median neutrophil-to-target distance"))
    raw_stats_rows.extend(uc_tests(q, "median_target_to_neutrophil_um", f"{target}; median target-to-neutrophil distance"))
    for radius in RADII:
        raw_stats_rows.extend(uc_tests(
            q, f"fraction_neutrophils_within_{radius}um",
            f"{target}; fraction neutrophils within {radius} um",
        ))
        raw_stats_rows.extend(uc_tests(
            q, f"fraction_target_fibroblasts_within_{radius}um_of_neutrophil",
            f"{target}; fraction target fibroblasts within {radius} um of neutrophil",
        ))
raw_stats = pd.DataFrame(raw_stats_rows)
raw_stats["fdr"] = bh_adjust(raw_stats.p_value)
raw_stats.to_csv(TAB / "04_raw_proximity_statistics.csv", index=False)

# Equal-neutrophil, abundance-adjusted label-permutation bootstrap.
rng = np.random.default_rng(SEED)
boot_rows = []
eligibility_rows = []
for pid, z in cells.groupby("PatientID", sort=True):
    fib = z[z.is_fibroblast].reset_index(drop=True)
    neut = z[z.is_neutrophil].reset_index(drop=True)
    eligible = len(neut) >= N_NEUTROPHILS
    eligibility_rows.append({
        "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
        "n_fibroblasts": len(fib), "n_neutrophils": len(neut), "eligible": eligible,
    })
    if not eligible:
        continue
    fib_xy = fib[["x", "y"]].to_numpy()
    neut_xy = neut[["x", "y"]].to_numpy()
    target_masks = {target: fib.fibroblast_state.isin(states).to_numpy() for target, states in TARGETS.items()}
    observed_distances = {
        target: cKDTree(fib_xy[mask]).query(neut_xy)[0]
        for target, mask in target_masks.items()
    }
    for iteration in range(1, N_BOOT + 1):
        ni = rng.choice(len(neut), N_NEUTROPHILS, replace=False)
        selected_neut_xy = neut_xy[ni]
        selected_neut_tree = cKDTree(selected_neut_xy)
        for target, mask in target_masks.items():
            n_target = int(mask.sum())
            observed_n_to_f = observed_distances[target][ni]
            observed_f_to_n = selected_neut_tree.query(fib_xy[mask])[0]
            random_idx = rng.choice(len(fib), n_target, replace=False)
            random_xy = fib_xy[random_idx]
            expected_n_to_f = cKDTree(random_xy).query(selected_neut_xy)[0]
            expected_f_to_n = selected_neut_tree.query(random_xy)[0]
            row = {
                "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
                "iteration": iteration, "target": target,
                "sampled_neutrophils": N_NEUTROPHILS, "n_target_fibroblasts": n_target,
                "observed_median_neutrophil_to_target_um": np.median(observed_n_to_f),
                "expected_median_neutrophil_to_random_fibroblast_um": np.median(expected_n_to_f),
                "log2_neutrophil_to_target_distance_ratio": np.log2(
                    (np.median(observed_n_to_f) + .5) / (np.median(expected_n_to_f) + .5)
                ),
                "observed_median_target_to_neutrophil_um": np.median(observed_f_to_n),
                "expected_median_random_fibroblast_to_neutrophil_um": np.median(expected_f_to_n),
                "log2_target_to_neutrophil_distance_ratio": np.log2(
                    (np.median(observed_f_to_n) + .5) / (np.median(expected_f_to_n) + .5)
                ),
            }
            for radius in RADII:
                observed_fraction = np.mean(observed_n_to_f <= radius)
                expected_fraction = np.mean(expected_n_to_f <= radius)
                row[f"observed_fraction_neutrophils_within_{radius}um"] = observed_fraction
                row[f"expected_fraction_neutrophils_within_{radius}um"] = expected_fraction
                row[f"difference_fraction_neutrophils_within_{radius}um"] = observed_fraction - expected_fraction
            boot_rows.append(row)

eligibility = pd.DataFrame(eligibility_rows)
eligibility.to_csv(TAB / "05_equal_count_bootstrap_eligibility.csv", index=False)
boot = pd.DataFrame(boot_rows)
boot.to_csv(TAB / "06_equal_count_proximity_bootstrap_iterations.csv", index=False)

agg_spec = {
    "n_iterations": ("iteration", "nunique"),
    "median_log2_neutrophil_to_target_distance_ratio": ("log2_neutrophil_to_target_distance_ratio", "median"),
    "median_log2_target_to_neutrophil_distance_ratio": ("log2_target_to_neutrophil_distance_ratio", "median"),
}
for radius in RADII:
    agg_spec[f"median_difference_fraction_within_{radius}um"] = (
        f"difference_fraction_neutrophils_within_{radius}um", "median"
    )
patient_boot = boot.groupby(["PatientID", "Diagnosis2", "target"]).agg(**agg_spec).reset_index()
patient_boot.to_csv(TAB / "07_equal_count_abundance_adjusted_proximity_by_patient.csv", index=False)

normalized_stats_rows = []
for target, q in patient_boot.groupby("target"):
    normalized_stats_rows.extend(uc_tests(
        q, "median_log2_neutrophil_to_target_distance_ratio",
        f"{target}; normalized neutrophil-to-target distance ratio", zero_tests=True,
    ))
    normalized_stats_rows.extend(uc_tests(
        q, "median_log2_target_to_neutrophil_distance_ratio",
        f"{target}; normalized target-to-neutrophil distance ratio", zero_tests=True,
    ))
    for radius in RADII:
        normalized_stats_rows.extend(uc_tests(
            q, f"median_difference_fraction_within_{radius}um",
            f"{target}; observed-minus-random fraction within {radius} um", zero_tests=True,
        ))
normalized_stats = pd.DataFrame(normalized_stats_rows)
normalized_stats["fdr_within_analysis"] = normalized_stats.groupby("analysis").p_value.transform(bh_adjust)
normalized_stats.to_csv(TAB / "08_equal_count_abundance_adjusted_proximity_statistics.csv", index=False)

# Paired FAP+ versus a5B1+ normalized distance comparison.
paired = patient_boot[patient_boot.target.isin(("FAP+ fibroblasts", "a5B1+ fibroblasts"))].pivot(
    index=["PatientID", "Diagnosis2"], columns="target",
    values="median_log2_neutrophil_to_target_distance_ratio",
).reset_index()
paired["FAP_minus_a5B1_normalized_distance"] = paired["FAP+ fibroblasts"] - paired["a5B1+ fibroblasts"]
paired.to_csv(TAB / "09_paired_FAP_vs_a5B1_normalized_distance.csv", index=False)
paired_rows = []
for group in GROUPS:
    x = paired.loc[paired.Diagnosis2.eq(group), "FAP_minus_a5B1_normalized_distance"].to_numpy()
    test = wilcoxon(x, alternative="two-sided", zero_method="wilcox")
    paired_rows.append({
        "analysis": "Within-patient FAP+ minus a5B1+ normalized distance",
        "group": group, "n_patients": len(x), "median_effect": np.median(x),
        "p_value": test.pvalue,
    })
x = paired.loc[paired.Diagnosis2.eq("UC_Inflamed"), "FAP_minus_a5B1_normalized_distance"].to_numpy()
y = paired.loc[paired.Diagnosis2.eq("UC_Noninflamed"), "FAP_minus_a5B1_normalized_distance"].to_numpy()
test = mannwhitneyu(x, y, alternative="two-sided")
paired_rows.append({
    "analysis": "Inflammation difference in paired FAP-versus-a5B1 distance",
    "group": "UC inflamed vs UC noninflamed", "n_patients": len(x) + len(y),
    "median_effect": np.median(x) - np.median(y), "p_value": test.pvalue,
})
paired_stats = pd.DataFrame(paired_rows)
paired_stats["fdr"] = bh_adjust(paired_stats.p_value)
paired_stats.to_csv(TAB / "10_paired_FAP_vs_a5B1_normalized_distance_statistics.csv", index=False)

# Figures.
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})
order = [GROUP_LABEL[g] for g in GROUPS]
palette = {GROUP_LABEL[g]: GROUP_COLOR[g] for g in GROUPS}
fig, axes = plt.subplots(2, 3, figsize=(21, 13))

q = patient_metrics.copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="target", y="median_neutrophil_to_target_um", hue="Group",
            hue_order=order, palette=palette, showfliers=False, ax=axes[0, 0])
axes[0, 0].set_title("A  Raw nearest-target distance")
axes[0, 0].set_xlabel("")
axes[0, 0].set_ylabel("Patient median distance (um)")
axes[0, 0].tick_params(axis="x", rotation=25, labelsize=9)
axes[0, 0].legend(title="", frameon=False, fontsize=9)

long_fraction = patient_metrics.melt(
    id_vars=["PatientID", "Diagnosis2", "target"],
    value_vars=[f"fraction_neutrophils_within_{r}um" for r in RADII],
    var_name="radius", value_name="fraction",
)
long_fraction["radius_um"] = long_fraction.radius.str.extract(r"(\d+)").astype(int)
long_fraction["Group"] = long_fraction.Diagnosis2.map(GROUP_LABEL)
sns.lineplot(data=long_fraction, x="radius_um", y="fraction", hue="Group", style="target",
             hue_order=order, palette=palette, markers=True, dashes=False, errorbar=None, ax=axes[0, 1])
axes[0, 1].set_title("B  Neutrophil proximity curves")
axes[0, 1].set_xlabel("Radius (um)")
axes[0, 1].set_ylabel("Fraction of neutrophils within radius")
axes[0, 1].legend(frameon=False, fontsize=8)

q = patient_boot.copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="target", y="median_log2_neutrophil_to_target_distance_ratio",
            hue="Group", hue_order=order, palette=palette, showfliers=False, ax=axes[0, 2])
axes[0, 2].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0, 2].set_title("C  Abundance-adjusted nearest distance")
axes[0, 2].set_xlabel("")
axes[0, 2].set_ylabel("log2 observed/random distance")
axes[0, 2].tick_params(axis="x", rotation=25, labelsize=9)
axes[0, 2].legend(title="", frameon=False, fontsize=9)

contact_heat = patient_boot.groupby(["Diagnosis2", "target"]).median_difference_fraction_within_50um.median().unstack()
contact_heat = contact_heat.reindex(index=GROUPS, columns=TARGETS)
sns.heatmap(contact_heat, annot=True, fmt=".2f", cmap="vlag", center=0,
            cbar_kws={"label": "Observed - random fraction"}, ax=axes[1, 0])
axes[1, 0].set_yticklabels(order, rotation=0)
axes[1, 0].set_title("D  Normalized 50 um contact fraction")
axes[1, 0].set_xlabel("")
axes[1, 0].set_ylabel("")

reciprocal_heat = patient_boot.groupby(["Diagnosis2", "target"]).median_log2_target_to_neutrophil_distance_ratio.median().unstack()
reciprocal_heat = reciprocal_heat.reindex(index=GROUPS, columns=TARGETS)
sns.heatmap(reciprocal_heat, annot=True, fmt=".2f", cmap="vlag", center=0,
            cbar_kws={"label": "log2 observed/random distance"}, ax=axes[1, 1])
axes[1, 1].set_yticklabels(order, rotation=0)
axes[1, 1].set_title("E  Fibroblast-to-neutrophil distance")
axes[1, 1].set_xlabel("")
axes[1, 1].set_ylabel("")

for _, row in paired.iterrows():
    axes[1, 2].plot(
        [0, 1], [row["FAP+ fibroblasts"], row["a5B1+ fibroblasts"]],
        color=GROUP_COLOR[row.Diagnosis2], alpha=.55, marker="o",
    )
axes[1, 2].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1, 2].set_xticks([0, 1], ["FAP+", "a5B1+"])
axes[1, 2].set_ylabel("log2 observed/random nearest distance")
axes[1, 2].set_title("F  Within-patient target proximity")

fig.suptitle("Neutrophil proximity to FAP+ and a5B1+ fibroblasts in UC", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_UC_FAP_a5B1_Proximity_Analyses")


def fmt(x):
    return f"{x:.3g}"


primary = normalized_stats[
    normalized_stats.endpoint.str.contains("normalized neutrophil-to-target distance ratio", regex=False)
]
between = primary[primary.analysis.eq("UC inflamed versus UC noninflamed")].copy()
between["target"] = between.endpoint.str.split(";").str[0]
between = between.set_index("target")
within = primary[primary.analysis.eq("Normalized proximity versus random expectation")].copy()
within["target"] = within.endpoint.str.split(";").str[0]

proximity_lines = []
for target in TARGETS:
    row = between.loc[target]
    proximity_lines.append(
        f"- {target}: UC noninflamed {row.median_UC_noninflamed:.3f}; "
        f"UC inflamed {row.median_UC_inflamed:.3f}; inflamed-minus-noninflamed "
        f"{row.effect_inflamed_minus_noninflamed:.3f}, P={fmt(row.p_value)}, "
        f"FDR={fmt(row.fdr_within_analysis)}."
    )
within_lines = []
for target in TARGETS:
    for group in GROUPS:
        row = within[within.target.eq(target) & within.group.eq(group)].iloc[0]
        direction = "farther than random" if row.effect_inflamed_minus_noninflamed > 0 else "closer than random"
        within_lines.append(
            f"- {target}, {GROUP_LABEL[group]}: {row.effect_inflamed_minus_noninflamed:.3f} "
            f"({direction}), P={fmt(row.p_value)}, FDR={fmt(row.fdr_within_analysis)}."
        )

raw_distance = raw_stats[
    raw_stats.endpoint.str.contains("median neutrophil-to-target distance", regex=False)
].set_index(raw_stats[
    raw_stats.endpoint.str.contains("median neutrophil-to-target distance", regex=False)
].endpoint.str.split(";").str[0])
raw_lines = []
for target in TARGETS:
    row = raw_distance.loc[target]
    raw_lines.append(
        f"- {target}: UC noninflamed {row.median_UC_noninflamed:.1f} um; "
        f"UC inflamed {row.median_UC_inflamed:.1f} um; P={fmt(row.p_value)}, FDR={fmt(row.fdr)}."
    )

paired_lines = [
    f"- {r.group}: median FAP-minus-a5B1 normalized distance {r.median_effect:.3f}, "
    f"P={fmt(r.p_value)}, FDR={fmt(r.fdr)}."
    for r in paired_stats.itertuples()
]

contact50 = normalized_stats[
    normalized_stats.endpoint.str.contains("observed-minus-random fraction within 50 um", regex=False)
    & normalized_stats.analysis.eq("Normalized proximity versus random expectation")
].copy()
contact50["target"] = contact50.endpoint.str.split(";").str[0]
contact_lines = []
for target in TARGETS:
    for group in GROUPS:
        row = contact50[contact50.target.eq(target) & contact50.group.eq(group)].iloc[0]
        contact_lines.append(
            f"- {target}, {GROUP_LABEL[group]}: observed-minus-random fraction "
            f"{row.effect_inflamed_minus_noninflamed:.3f}, P={fmt(row.p_value)}, "
            f"FDR={fmt(row.fdr_within_analysis)}."
        )

excluded = eligibility[~eligibility.eligible]
excluded_text = ", ".join(f"{r.PatientID} ({int(r.n_neutrophils)} neutrophils)" for r in excluded.itertuples()) or "none"

report = f"""# Neutrophil proximity to FAP+ and a5B1+ fibroblasts in UC

## Methods

FAP+ and a5B1+ were evaluated as marginal targets, with the overlapping FAP+/a5B1+ subset reported separately. Raw nearest-neighbor distances and fractions within 25, 50, and 100 um were calculated in both directions. The primary abundance-adjusted analysis used {N_BOOT} iterations of exactly {N_NEUTROPHILS} neutrophils per sample and compared each observed target to an equal-sized random fibroblast subset from the same patient. Positive normalized distance indicates avoidance; negative values indicate closer-than-random proximity. Patient was the biological replicate. Excluded from equal-count analysis: {excluded_text}.

## Raw nearest-neighbor distances

{chr(10).join(raw_lines)}

## Equal-count, abundance-adjusted proximity

{chr(10).join(proximity_lines)}

Within-group comparisons against the random-label expectation were:

{chr(10).join(within_lines)}

At 50 um, normalized contact-fraction deficits were:

{chr(10).join(contact_lines)}

## Paired FAP+ versus a5B1+ proximity

{chr(10).join(paired_lines)}

Complete bidirectional-distance and 25/50/100 um contact-fraction statistics are included in the accompanying tables. These analyses quantify spatial proximity, not direct contact, signaling, or causality.
"""
(OUT / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

manifest = {
    "targets": TARGETS, "radii_um": list(RADII),
    "equal_neutrophils_per_iteration": N_NEUTROPHILS,
    "bootstrap_permutation_iterations": N_BOOT,
    "randomization": "equal-sized random fibroblast subset within patient",
    "excluded_equal_count_samples": excluded.PatientID.tolist(),
    "seed": SEED, "biological_replicate": "patient",
}
(OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print("Completed UC FAP/a5B1 proximity analyses.")
print(f"Saved outputs to {OUT}")
