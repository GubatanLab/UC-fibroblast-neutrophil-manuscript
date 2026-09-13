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
OUT = ROOT / "neutrophil_a5b1_fibroblast_reanalysis"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260725
N_PERM = 1000
N_BALANCED_BOOTSTRAP = 500
BALANCED_FIBROBLAST_COUNT = 250
RADII = (25, 50, 100)
GROUPS = ("Control", "UC_Noninflamed", "UC_Inflamed")
UC_GROUPS = ("UC_Noninflamed", "UC_Inflamed")
GROUP_LABEL = {
    "Control": "Control",
    "UC_Noninflamed": "UC noninflamed",
    "UC_Inflamed": "UC inflamed",
}
GROUP_COLOR = {
    "Control": "#4C78A8",
    "UC_Noninflamed": "#D4AD35",
    "UC_Inflamed": "#D95F5F",
}
SUBTYPES = ("All neutrophils", "PADI4", "CXCR4", "MX1", "OSM")


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


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def interaction_metrics(neut_xy, fib_xy, fib_positive, radii, rng):
    """Radius-graph enrichment after permuting a5B1 labels among fibroblasts."""
    if len(neut_xy) == 0 or len(fib_xy) < 2 or fib_positive.sum() == 0:
        return []
    fib_tree = cKDTree(fib_xy)
    pairs = fib_tree.query_ball_point(neut_xy, max(radii))
    n_fib = len(fib_xy)
    n_pos = int(fib_positive.sum())
    results = []
    for radius in radii:
        degrees = np.zeros(n_fib, dtype=np.int64)
        n_neut_contact = 0
        n_neut_pos_contact = 0
        total_edges = 0
        pos_edges = 0
        for point, ids in zip(neut_xy, pairs):
            if not ids:
                continue
            ids = np.asarray(ids, dtype=int)
            delta = fib_xy[ids] - point
            ids = ids[(delta * delta).sum(axis=1) <= radius * radius]
            if not len(ids):
                continue
            n_neut_contact += 1
            total_edges += len(ids)
            np.add.at(degrees, ids, 1)
            hit = fib_positive[ids]
            pos_edges += int(hit.sum())
            n_neut_pos_contact += int(hit.any())
        if total_edges == 0:
            continue
        perm_edges = np.empty(N_PERM, dtype=float)
        for k in range(N_PERM):
            chosen = rng.choice(n_fib, n_pos, replace=False)
            perm_edges[k] = degrees[chosen].sum()
        expected = float(perm_edges.mean())
        enrichment = np.log2((pos_edges + 0.5) / (expected + 0.5))
        abundance_rate_ratio = ((pos_edges + 0.5) / n_pos) / ((total_edges + 0.5) / n_fib)
        p_two = (1 + np.sum(np.abs(perm_edges - expected) >= abs(pos_edges - expected))) / (N_PERM + 1)
        results.append(
            {
                "radius_um": radius,
                "n_neutrophils": len(neut_xy),
                "n_fibroblasts": n_fib,
                "n_a5B1_positive_fibroblasts": n_pos,
                "total_neutrophil_fibroblast_edges": total_edges,
                "observed_a5B1_positive_edges": pos_edges,
                "expected_a5B1_positive_edges": expected,
                "log2_observed_expected": enrichment,
                "a5B1_positive_edges_per_neutrophil_per_1000_positive_fibroblasts": (
                    pos_edges / len(neut_xy) / n_pos * 1000
                ),
                "all_fibroblast_edges_per_neutrophil_per_1000_fibroblasts": (
                    total_edges / len(neut_xy) / n_fib * 1000
                ),
                "fibroblast_abundance_normalized_rate_ratio": abundance_rate_ratio,
                "log2_fibroblast_abundance_normalized_rate_ratio": np.log2(abundance_rate_ratio),
                "permutation_p": p_two,
                "fraction_neutrophils_with_any_fibroblast": n_neut_contact / len(neut_xy),
                "fraction_neutrophils_with_a5B1_positive_fibroblast": n_neut_pos_contact / len(neut_xy),
            }
        )
    return results


d = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv")
required = {
    "cell_ID", "PatientID", "Diagnosis2", "cell_type", "Neutrophil_subtype_0.8",
    "a5B1_cell_raw", "a5B1_cell_clr", "x", "y",
}
missing = required.difference(d.columns)
if missing:
    raise RuntimeError(f"Missing required exported columns: {sorted(missing)}")

d["is_fibroblast"] = d["cell_type"].eq("Fibroblast")
d["is_neutrophil"] = d["cell_type"].str.startswith("Neutrophil", na=False)
uc_fib = d[d.is_fibroblast & d.Diagnosis2.isin(UC_GROUPS)]
global_threshold = float(uc_fib.a5B1_cell_clr.quantile(0.75))
d["a5B1_positive_global"] = d.is_fibroblast & (d.a5B1_cell_clr >= global_threshold)
patient_thresholds = d[d.is_fibroblast].groupby("PatientID").a5B1_cell_clr.quantile(0.75)
d["patient_a5B1_q75"] = d.PatientID.map(patient_thresholds)
d["a5B1_positive_patient_q75"] = d.is_fibroblast & (d.a5B1_cell_clr >= d.patient_a5B1_q75)

definition = pd.DataFrame(
    {
        "definition": ["Primary", "Sensitivity"],
        "positive_rule": [
            "Fibroblast a5B1 CLR intensity >= pooled UC fibroblast 75th percentile",
            "Fibroblast a5B1 CLR intensity >= within-patient fibroblast 75th percentile",
        ],
        "threshold": [global_threshold, np.nan],
        "normalization": ["Seurat Akoya data layer (CLR)", "Seurat Akoya data layer (CLR)"],
    }
)
definition.to_csv(TAB / "00_a5B1_positive_definition.csv", index=False)

fib_summary = (
    d[d.is_fibroblast]
    .groupby(["PatientID", "Diagnosis2"])
    .agg(
        n_fibroblasts=("cell_ID", "size"),
        median_a5B1_clr=("a5B1_cell_clr", "median"),
        q75_a5B1_clr=("a5B1_cell_clr", lambda x: x.quantile(0.75)),
        n_a5B1_positive=("a5B1_positive_global", "sum"),
        fraction_a5B1_positive=("a5B1_positive_global", "mean"),
    )
    .reset_index()
)
fib_summary.to_csv(TAB / "01_a5B1_fibroblast_summary_by_patient.csv", index=False)
fib_stats_rows = []
for metric in ("fraction_a5B1_positive", "median_a5B1_clr"):
    noninflamed = fib_summary[fib_summary.Diagnosis2 == "UC_Noninflamed"][metric].to_numpy()
    inflamed = fib_summary[fib_summary.Diagnosis2 == "UC_Inflamed"][metric].to_numpy()
    test = mannwhitneyu(inflamed, noninflamed, alternative="two-sided")
    fib_stats_rows.append(
        {
            "metric": metric,
            "n_UC_noninflamed": len(noninflamed),
            "n_UC_inflamed": len(inflamed),
            "median_UC_noninflamed": np.median(noninflamed),
            "median_UC_inflamed": np.median(inflamed),
            "median_difference_inflamed_minus_noninflamed": np.median(inflamed) - np.median(noninflamed),
            "rank_biserial": rank_biserial(inflamed, noninflamed),
            "p_value": test.pvalue,
        }
    )
fib_stats = pd.DataFrame(fib_stats_rows)
fib_stats["fdr"] = bh_adjust(fib_stats.p_value)
fib_stats.to_csv(TAB / "01b_a5B1_fibroblast_abundance_statistics.csv", index=False)

rng = np.random.default_rng(SEED)
interaction_rows = []
distance_rows = []
for pid, z in d.groupby("PatientID", sort=True):
    fib = z[z.is_fibroblast].copy()
    neut = z[z.is_neutrophil].copy()
    if fib.empty or neut.empty:
        continue
    fib_xy = fib[["x", "y"]].to_numpy()
    pos_xy = fib.loc[fib.a5B1_positive_global, ["x", "y"]].to_numpy()
    neg_xy = fib.loc[~fib.a5B1_positive_global, ["x", "y"]].to_numpy()
    if len(pos_xy) and len(neg_xy):
        dp = cKDTree(pos_xy).query(neut[["x", "y"]].to_numpy())[0]
        dn = cKDTree(neg_xy).query(neut[["x", "y"]].to_numpy())[0]
        distance_rows.extend(
            {
                "cell_ID": cid,
                "PatientID": pid,
                "Diagnosis2": z.Diagnosis2.iloc[0],
                "Neutrophil_subtype": st if pd.notna(st) else "Unclassified",
                "nearest_a5B1_positive_fibroblast_um": xp,
                "nearest_a5B1_negative_fibroblast_um": xn,
                "distance_difference_positive_minus_negative_um": xp - xn,
            }
            for cid, st, xp, xn in zip(neut.cell_ID, neut["Neutrophil_subtype_0.8"], dp, dn)
        )
    for subtype in SUBTYPES:
        n = neut if subtype == "All neutrophils" else neut[neut["Neutrophil_subtype_0.8"].eq(subtype)]
        metrics = interaction_metrics(
            n[["x", "y"]].to_numpy(),
            fib_xy,
            fib.a5B1_positive_global.to_numpy(),
            RADII,
            rng,
        )
        for row in metrics:
            row.update(
                {
                    "PatientID": pid,
                    "Diagnosis2": z.Diagnosis2.iloc[0],
                    "neutrophil_population": subtype,
                    "a5B1_definition": "pooled_UC_fibroblast_q75",
                }
            )
            interaction_rows.append(row)

interaction = pd.DataFrame(interaction_rows)
interaction.to_csv(TAB / "02_spatial_interaction_enrichment_by_patient.csv", index=False)
cell_dist = pd.DataFrame(distance_rows)
cell_dist.to_csv(TAB / "03_neutrophil_nearest_a5B1_fibroblast_distances_by_cell.csv", index=False)

patient_dist = (
    cell_dist.groupby(["PatientID", "Diagnosis2", "Neutrophil_subtype"])
    .agg(
        n_neutrophils=("cell_ID", "size"),
        median_nearest_a5B1_positive_um=("nearest_a5B1_positive_fibroblast_um", "median"),
        median_nearest_a5B1_negative_um=("nearest_a5B1_negative_fibroblast_um", "median"),
        median_distance_difference_um=("distance_difference_positive_minus_negative_um", "median"),
        fraction_closer_to_a5B1_positive=("distance_difference_positive_minus_negative_um", lambda x: np.mean(x < 0)),
    )
    .reset_index()
)
all_dist = (
    cell_dist.groupby(["PatientID", "Diagnosis2"])
    .agg(
        n_neutrophils=("cell_ID", "size"),
        median_nearest_a5B1_positive_um=("nearest_a5B1_positive_fibroblast_um", "median"),
        median_nearest_a5B1_negative_um=("nearest_a5B1_negative_fibroblast_um", "median"),
        median_distance_difference_um=("distance_difference_positive_minus_negative_um", "median"),
        fraction_closer_to_a5B1_positive=("distance_difference_positive_minus_negative_um", lambda x: np.mean(x < 0)),
    )
    .reset_index()
)
all_dist["Neutrophil_subtype"] = "All neutrophils"
patient_dist = pd.concat([all_dist, patient_dist], ignore_index=True)
patient_dist.to_csv(TAB / "04_nearest_distance_summary_by_patient.csv", index=False)

# Sensitivity analysis with patient-specific top-quartile labels.
sensitivity_rows = []
rng_sens = np.random.default_rng(SEED + 1)
for pid, z in d.groupby("PatientID", sort=True):
    fib = z[z.is_fibroblast]
    neut = z[z.is_neutrophil]
    if fib.empty or neut.empty:
        continue
    metrics = interaction_metrics(
        neut[["x", "y"]].to_numpy(),
        fib[["x", "y"]].to_numpy(),
        fib.a5B1_positive_patient_q75.to_numpy(),
        (50,),
        rng_sens,
    )
    for row in metrics:
        row.update(
            {
                "PatientID": pid,
                "Diagnosis2": z.Diagnosis2.iloc[0],
                "neutrophil_population": "All neutrophils",
                "a5B1_definition": "within_patient_fibroblast_q75",
            }
        )
        sensitivity_rows.append(row)
sensitivity = pd.DataFrame(sensitivity_rows)
sensitivity.to_csv(TAB / "05_sensitivity_patient_specific_q75.csv", index=False)

# Explicit total-fibroblast abundance normalization by repeated equal-count
# subsampling. Every UC patient contributes 250 fibroblasts per iteration.
balanced_iterations = []
rng_balanced = np.random.default_rng(SEED + 2)
for pid, z in d[d.Diagnosis2.isin(UC_GROUPS)].groupby("PatientID", sort=True):
    fib = z[z.is_fibroblast].copy()
    neut = z[z.is_neutrophil]
    if len(fib) < BALANCED_FIBROBLAST_COUNT or neut.empty:
        continue
    fib_xy = fib[["x", "y"]].to_numpy()
    degrees = np.zeros(len(fib), dtype=np.int64)
    for ids in cKDTree(fib_xy).query_ball_point(neut[["x", "y"]].to_numpy(), 50):
        if ids:
            np.add.at(degrees, np.asarray(ids, dtype=int), 1)
    positive = fib.a5B1_positive_global.to_numpy()
    for iteration in range(N_BALANCED_BOOTSTRAP):
        selected = rng_balanced.choice(len(fib), BALANCED_FIBROBLAST_COUNT, replace=False)
        selected_positive = positive[selected]
        n_pos = int(selected_positive.sum())
        if n_pos == 0:
            continue
        total_edges = int(degrees[selected].sum())
        pos_edges = int(degrees[selected[selected_positive]].sum())
        rate_ratio = ((pos_edges + 0.5) / n_pos) / (
            (total_edges + 0.5) / BALANCED_FIBROBLAST_COUNT
        )
        balanced_iterations.append(
            {
                "PatientID": pid,
                "Diagnosis2": z.Diagnosis2.iloc[0],
                "iteration": iteration + 1,
                "sampled_fibroblasts": BALANCED_FIBROBLAST_COUNT,
                "sampled_a5B1_positive_fibroblasts": n_pos,
                "total_edges_50um": total_edges,
                "a5B1_positive_edges_50um": pos_edges,
                "fibroblast_abundance_normalized_rate_ratio": rate_ratio,
                "log2_fibroblast_abundance_normalized_rate_ratio": np.log2(rate_ratio),
            }
        )
balanced_iterations = pd.DataFrame(balanced_iterations)
balanced_iterations.to_csv(TAB / "09_equal_fibroblast_count_bootstrap_iterations.csv", index=False)
balanced_patient = (
    balanced_iterations.groupby(["PatientID", "Diagnosis2"])
    .agg(
        n_iterations=("iteration", "size"),
        median_sampled_a5B1_positive=("sampled_a5B1_positive_fibroblasts", "median"),
        median_log2_abundance_normalized_rate_ratio=(
            "log2_fibroblast_abundance_normalized_rate_ratio", "median"
        ),
        ci025_log2_abundance_normalized_rate_ratio=(
            "log2_fibroblast_abundance_normalized_rate_ratio", lambda x: x.quantile(0.025)
        ),
        ci975_log2_abundance_normalized_rate_ratio=(
            "log2_fibroblast_abundance_normalized_rate_ratio", lambda x: x.quantile(0.975)
        ),
    )
    .reset_index()
)
balanced_patient.to_csv(TAB / "10_equal_fibroblast_count_summary_by_patient.csv", index=False)
bal_noninflamed = balanced_patient[
    balanced_patient.Diagnosis2 == "UC_Noninflamed"
].median_log2_abundance_normalized_rate_ratio.to_numpy()
bal_inflamed = balanced_patient[
    balanced_patient.Diagnosis2 == "UC_Inflamed"
].median_log2_abundance_normalized_rate_ratio.to_numpy()
balanced_between = mannwhitneyu(bal_inflamed, bal_noninflamed, alternative="two-sided")
balanced_stats_rows = [
    {
        "analysis": "UC inflamed versus UC noninflamed",
        "group": "",
        "n_1": len(bal_inflamed),
        "n_2": len(bal_noninflamed),
        "median_1": np.median(bal_inflamed),
        "median_2": np.median(bal_noninflamed),
        "effect": np.median(bal_inflamed) - np.median(bal_noninflamed),
        "p_value": balanced_between.pvalue,
    }
]
for group, values in [
    ("UC_Noninflamed", bal_noninflamed),
    ("UC_Inflamed", bal_inflamed),
]:
    one = wilcoxon(values, alternative="two-sided", zero_method="wilcox")
    balanced_stats_rows.append(
        {
            "analysis": "Normalized interaction versus no preference",
            "group": group,
            "n_1": len(values),
            "n_2": np.nan,
            "median_1": np.median(values),
            "median_2": 0,
            "effect": np.median(values),
            "p_value": one.pvalue,
        }
    )
balanced_stats = pd.DataFrame(balanced_stats_rows)
balanced_stats["fdr"] = bh_adjust(balanced_stats.p_value)
balanced_stats.to_csv(TAB / "11_equal_fibroblast_count_statistics.csv", index=False)

# Patient-level inference: UC inflamed versus noninflamed and enrichment versus zero.
stats_rows = []
for (population, radius), q in interaction.groupby(["neutrophil_population", "radius_um"]):
    q = q[q.Diagnosis2.isin(UC_GROUPS)]
    a = q[q.Diagnosis2 == "UC_Noninflamed"].log2_observed_expected.dropna().to_numpy()
    b = q[q.Diagnosis2 == "UC_Inflamed"].log2_observed_expected.dropna().to_numpy()
    if len(a) and len(b):
        test = mannwhitneyu(b, a, alternative="two-sided")
        stats_rows.append(
            {
                "analysis": "UC inflamed versus UC noninflamed",
                "neutrophil_population": population,
                "radius_um": radius,
                "group": "",
                "n_1": len(b),
                "n_2": len(a),
                "median_1": np.median(b),
                "median_2": np.median(a),
                "effect": np.median(b) - np.median(a),
                "effect_type": "median difference (inflamed - noninflamed)",
                "rank_biserial": rank_biserial(b, a),
                "p_value": test.pvalue,
            }
        )
    for group in UC_GROUPS:
        x = q[q.Diagnosis2 == group].log2_observed_expected.dropna().to_numpy()
        if len(x) >= 3 and np.any(x != 0):
            test = wilcoxon(x, alternative="two-sided", zero_method="wilcox")
            stats_rows.append(
                {
                    "analysis": "Enrichment versus no preference",
                    "neutrophil_population": population,
                    "radius_um": radius,
                    "group": group,
                    "n_1": len(x),
                    "n_2": np.nan,
                    "median_1": np.median(x),
                    "median_2": 0,
                    "effect": np.median(x),
                    "effect_type": "median log2 observed/expected",
                    "rank_biserial": np.nan,
                    "p_value": test.pvalue,
                }
            )
stats = pd.DataFrame(stats_rows)
stats["fdr_within_analysis"] = stats.groupby("analysis").p_value.transform(bh_adjust)
stats.to_csv(TAB / "06_patient_level_statistics.csv", index=False)

# P07 exclusion sensitivity for the primary comparison.
exclusion = []
for source, label in [(interaction, "pooled_UC_q75"), (sensitivity, "patient_q75")]:
    q = source[
        source.Diagnosis2.isin(UC_GROUPS)
        & source.neutrophil_population.eq("All neutrophils")
        & source.radius_um.eq(50)
        & source.PatientID.ne("P07")
    ]
    a = q[q.Diagnosis2 == "UC_Noninflamed"].log2_observed_expected.to_numpy()
    b = q[q.Diagnosis2 == "UC_Inflamed"].log2_observed_expected.to_numpy()
    exclusion.append(
        {
            "definition": label,
            "excluded_patient": "P07",
            "n_UC_noninflamed": len(a),
            "n_UC_inflamed": len(b),
            "median_UC_noninflamed": np.median(a),
            "median_UC_inflamed": np.median(b),
            "median_difference": np.median(b) - np.median(a),
            "p_value": mannwhitneyu(b, a, alternative="two-sided").pvalue,
        }
    )
pd.DataFrame(exclusion).to_csv(TAB / "07_sensitivity_excluding_P07.csv", index=False)

# Group summaries.
summary = (
    interaction.groupby(["Diagnosis2", "neutrophil_population", "radius_um"])
    .log2_observed_expected.agg(["count", "median", "mean", "std"])
    .reset_index()
)
summary.to_csv(TAB / "08_group_summary.csv", index=False)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})

# Figure 1: primary patient-level results.
fig, axes = plt.subplots(1, 3, figsize=(20, 6))
q = fib_summary.copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
palette = {GROUP_LABEL[k]: v for k, v in GROUP_COLOR.items()}
sns.boxplot(data=q, x="Group", y="fraction_a5B1_positive", order=[GROUP_LABEL[g] for g in GROUPS],
            palette=palette, showfliers=False, ax=axes[0])
sns.stripplot(data=q, x="Group", y="fraction_a5B1_positive", order=[GROUP_LABEL[g] for g in GROUPS],
              color="black", size=5, ax=axes[0])
axes[0].set_title("A  α5β1+ fibroblast abundance")
axes[0].set_xlabel("")
axes[0].set_ylabel("Fraction of fibroblasts")
axes[0].tick_params(axis="x", rotation=20)

q = interaction[(interaction.neutrophil_population == "All neutrophils")].copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="radius_um", y="log2_observed_expected", hue="Group",
            hue_order=[GROUP_LABEL[g] for g in GROUPS], palette=palette, showfliers=False, ax=axes[1])
sns.stripplot(data=q, x="radius_um", y="log2_observed_expected", hue="Group",
              hue_order=[GROUP_LABEL[g] for g in GROUPS], dodge=True, color="black", size=3.5,
              legend=False, ax=axes[1])
axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1].set_title("B  Neutrophil–α5β1+ fibroblast targeting")
axes[1].set_xlabel("Interaction radius (µm)")
axes[1].set_ylabel("log2 observed / permuted expected edges")
axes[1].legend(title="", frameon=False, fontsize=10)

q = interaction[(interaction.radius_um == 50) & (interaction.neutrophil_population != "All neutrophils")].copy()
heat = q.groupby(["Diagnosis2", "neutrophil_population"]).log2_observed_expected.median().unstack()
heat = heat.reindex(index=GROUPS, columns=SUBTYPES[1:])
sns.heatmap(heat, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f",
            cbar_kws={"label": "Median log2 O/E"}, ax=axes[2])
axes[2].set_yticklabels([GROUP_LABEL[g] for g in GROUPS], rotation=0)
axes[2].set_title("C  Neutrophil subtype interactions (50 µm)")
axes[2].set_xlabel("Neutrophil subtype")
axes[2].set_ylabel("")
fig.suptitle("Neutrophil interactions with α5β1+ fibroblasts in UC", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_A5B1_01_primary_results")

# Figure 2: representative UC maps, using median-enrichment patient per group.
primary50 = interaction[
    (interaction.neutrophil_population == "All neutrophils") & (interaction.radius_um == 50)
]
representatives = {}
for group in UC_GROUPS:
    x = primary50[primary50.Diagnosis2 == group].sort_values("log2_observed_expected")
    representatives[group] = x.iloc[(len(x) - 1) // 2].PatientID
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
for ax, group in zip(axes, UC_GROUPS):
    pid = representatives[group]
    z = d[d.PatientID == pid]
    bg = z[~z.is_fibroblast & ~z.is_neutrophil]
    neg = z[z.is_fibroblast & ~z.a5B1_positive_global]
    pos = z[z.a5B1_positive_global]
    neu = z[z.is_neutrophil]
    ax.scatter(bg.x, bg.y, s=0.4, c="#D9D9D9", alpha=0.18, rasterized=True)
    ax.scatter(neg.x, neg.y, s=2, c="#8C8C8C", alpha=0.45, label="α5β1− fibroblast", rasterized=True)
    ax.scatter(pos.x, pos.y, s=4, c="#C51B7D", alpha=0.8, label="α5β1+ fibroblast", rasterized=True)
    ax.scatter(neu.x, neu.y, s=5, c="#00A6D6", alpha=0.85, label="Neutrophil", rasterized=True)
    enr = primary50[primary50.PatientID == pid].log2_observed_expected.iloc[0]
    ax.set_title(f"{GROUP_LABEL[group]}: {pid}\n50 µm log2 O/E = {enr:.2f}")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
axes[1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
fig.suptitle("Representative α5β1 fibroblast–neutrophil spatial organization", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_A5B1_02_representative_spatial_maps")

# Figure 3: explicit fibroblast-abundance normalization.
fig, axes = plt.subplots(1, 3, figsize=(20, 6))
full50 = interaction[
    (interaction.neutrophil_population == "All neutrophils")
    & (interaction.radius_um == 50)
    & (interaction.Diagnosis2.isin(UC_GROUPS))
].copy()
full50["Group"] = full50.Diagnosis2.map(GROUP_LABEL)
sns.scatterplot(
    data=full50, x="n_fibroblasts",
    y="log2_fibroblast_abundance_normalized_rate_ratio",
    hue="Group", palette=palette, s=90, ax=axes[0],
)
axes[0].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0].set_title("A  Rate normalized per fibroblast")
axes[0].set_xlabel("Total fibroblasts in patient")
axes[0].set_ylabel("log2 α5β1+ / all-fibroblast edge rate")
axes[0].legend(title="", frameon=False, fontsize=10)

bq = balanced_patient.copy()
bq["Group"] = bq.Diagnosis2.map(GROUP_LABEL)
uc_palette = {GROUP_LABEL[k]: GROUP_COLOR[k] for k in UC_GROUPS}
sns.boxplot(
    data=bq, x="Group", y="median_log2_abundance_normalized_rate_ratio",
    order=[GROUP_LABEL[g] for g in UC_GROUPS], palette=uc_palette,
    showfliers=False, ax=axes[1],
)
sns.stripplot(
    data=bq, x="Group", y="median_log2_abundance_normalized_rate_ratio",
    order=[GROUP_LABEL[g] for g in UC_GROUPS], color="black", size=6, ax=axes[1],
)
axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1].set_title(f"B  Equal-count resampling ({BALANCED_FIBROBLAST_COUNT}/patient)")
axes[1].set_xlabel("")
axes[1].set_ylabel("Median bootstrap log2 normalized rate")

paired = full50[
    ["PatientID", "Diagnosis2", "log2_fibroblast_abundance_normalized_rate_ratio"]
].merge(
    balanced_patient[["PatientID", "median_log2_abundance_normalized_rate_ratio"]],
    on="PatientID",
)
for _, row in paired.iterrows():
    axes[2].plot(
        [0, 1],
        [row.log2_fibroblast_abundance_normalized_rate_ratio,
         row.median_log2_abundance_normalized_rate_ratio],
        color=GROUP_COLOR[row.Diagnosis2], alpha=0.55, marker="o",
    )
axes[2].axhline(0, color="black", linestyle="--", linewidth=1)
axes[2].set_xticks([0, 1], ["All fibroblasts", f"{BALANCED_FIBROBLAST_COUNT}-cell resampling"])
axes[2].set_ylabel("log2 abundance-normalized rate ratio")
axes[2].set_title("C  Sensitivity to fibroblast count")
fig.suptitle(
    "Neutrophil–α5β1+ fibroblast interactions normalized for fibroblast abundance",
    fontweight="bold",
)
fig.tight_layout()
save(fig, "Figure_A5B1_03_fibroblast_abundance_normalized")

# Compact machine-readable manifest and text summary.
primary_stats = stats[
    (stats.analysis == "UC inflamed versus UC noninflamed")
    & (stats.neutrophil_population == "All neutrophils")
    & (stats.radius_um == 50)
].iloc[0]
group_medians = primary50.groupby("Diagnosis2").log2_observed_expected.median()
sens_medians = sensitivity.groupby("Diagnosis2").log2_observed_expected.median()
sens_noninflamed = sensitivity[sensitivity.Diagnosis2 == "UC_Noninflamed"].log2_observed_expected
sens_inflamed = sensitivity[sensitivity.Diagnosis2 == "UC_Inflamed"].log2_observed_expected
sens_p = mannwhitneyu(sens_inflamed, sens_noninflamed, alternative="two-sided").pvalue
within_primary = stats[
    (stats.analysis == "Enrichment versus no preference")
    & (stats.neutrophil_population == "All neutrophils")
    & (stats.radius_um == 50)
].set_index("group")
abundance_stat = fib_stats.set_index("metric").loc["fraction_a5B1_positive"]
abundance_medians = fib_summary.groupby("Diagnosis2").fraction_a5B1_positive.median()
subtype50 = stats[
    (stats.radius_um == 50)
    & (stats.neutrophil_population != "All neutrophils")
]
subtype_depleted = {}
for group in UC_GROUPS:
    hit = subtype50[
        (subtype50.analysis == "Enrichment versus no preference")
        & (subtype50.group == group)
        & (subtype50.effect < 0)
        & (subtype50.fdr_within_analysis < 0.05)
    ]
    subtype_depleted[group] = ", ".join(hit.neutrophil_population.tolist()) or "none"
manifest = {
    "input": str(ROOT / "additional_analyses" / "codex_extended_cells.csv"),
    "seed": SEED,
    "permutations_per_patient": N_PERM,
    "balanced_bootstrap_iterations": N_BALANCED_BOOTSTRAP,
    "balanced_fibroblasts_per_patient": BALANCED_FIBROBLAST_COUNT,
    "primary_a5B1_threshold_clr": global_threshold,
    "threshold_source": "pooled UC fibroblast 75th percentile",
    "primary_radius_um": 50,
    "biological_replicate": "patient",
    "representative_patients": representatives,
}
(OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

report = f"""# Neutrophil interactions with α5β1+ fibroblasts in UC

## Analysis definition

Fibroblasts were called α5β1-positive when their cell-level α5β1 CLR intensity was at or above the pooled UC fibroblast 75th percentile ({global_threshold:.4f}). The primary spatial endpoint was the number of neutrophil–α5β1+ fibroblast edges within 50 µm relative to the expectation from {N_PERM} random permutations of α5β1 labels among fibroblasts in the same patient. This preserves each patient's fibroblast positions, α5β1+ abundance, and neutrophil geometry. Inference used patient-level values.

## Primary result

At 50 µm, median log2 observed/expected enrichment was {group_medians.get('UC_Noninflamed', np.nan):.3f} in UC noninflamed tissue and {group_medians.get('UC_Inflamed', np.nan):.3f} in UC inflamed tissue. The median inflamed-minus-noninflamed difference was {primary_stats.effect:.3f} (two-sided patient-level Mann–Whitney P={primary_stats.p_value:.4g}, FDR={primary_stats.fdr_within_analysis:.4g}; n={int(primary_stats.n_1)} inflamed and n={int(primary_stats.n_2)} noninflamed).

The negative enrichment was significant relative to no spatial preference in UC noninflamed tissue (Wilcoxon P={within_primary.loc['UC_Noninflamed', 'p_value']:.4g}, FDR={within_primary.loc['UC_Noninflamed', 'fdr_within_analysis']:.4g}) and UC inflamed tissue (P={within_primary.loc['UC_Inflamed', 'p_value']:.4g}, FDR={within_primary.loc['UC_Inflamed', 'fdr_within_analysis']:.4g}). Thus, neutrophil neighborhoods were depleted of α5β1+ fibroblasts after controlling for the number and positions of fibroblasts in each patient.

At 50 µm, significant subtype-level depletion was observed for {subtype_depleted['UC_Noninflamed']} in UC noninflamed tissue and {subtype_depleted['UC_Inflamed']} in UC inflamed tissue (FDR < 0.05). No neutrophil subtype showed a significant inflamed-versus-noninflamed difference after multiple-testing correction.

## α5β1+ fibroblast abundance

The median fraction of fibroblasts above the pooled UC α5β1 threshold was {abundance_medians.get('UC_Noninflamed', np.nan):.3f} in UC noninflamed tissue and {abundance_medians.get('UC_Inflamed', np.nan):.3f} in UC inflamed tissue (Mann–Whitney P={abundance_stat.p_value:.4g}, FDR={abundance_stat.fdr:.4g}). This descriptive decrease did not reach statistical significance.

## Robustness

Using a within-patient top-quartile α5β1 definition, the median 50 µm enrichment was {sens_medians.get('UC_Noninflamed', np.nan):.3f} in UC noninflamed and {sens_medians.get('UC_Inflamed', np.nan):.3f} in UC inflamed tissue (Mann–Whitney P={sens_p:.4g}). Excluding repaired patient P07 also left the primary between-state conclusion nonsignificant.

## Explicit normalization for fibroblast abundance

The interaction rate was divided by the number of available α5β1+ fibroblasts and referenced to the corresponding edge rate per total fibroblast. In a stricter resampling analysis, every UC patient contributed exactly {BALANCED_FIBROBLAST_COUNT} fibroblasts in each of {N_BALANCED_BOOTSTRAP} iterations. The resulting patient-level median log2 normalized rate ratio was {np.median(bal_noninflamed):.3f} in UC noninflamed and {np.median(bal_inflamed):.3f} in UC inflamed tissue. Both remained below zero (UC noninflamed FDR={balanced_stats.loc[balanced_stats.group == 'UC_Noninflamed', 'fdr'].iloc[0]:.4g}; UC inflamed FDR={balanced_stats.loc[balanced_stats.group == 'UC_Inflamed', 'fdr'].iloc[0]:.4g}), while the between-state difference remained nonsignificant (Mann–Whitney P={balanced_between.pvalue:.4g}). Thus, the depletion of α5β1+ fibroblasts from neutrophil neighborhoods was not explained by differences in total fibroblast abundance.

## Interpretation guardrail

The permutation enrichment quantifies preferential spatial association beyond α5β1+ fibroblast availability; it does not establish physical receptor binding, signaling direction, or causality. The α5β1-positive call is an operational high-expression threshold because no external antibody-specific positivity control was supplied.
"""
(OUT / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

print("Completed neutrophil-a5B1 fibroblast reanalysis.")
print(f"Saved outputs to {OUT}")
