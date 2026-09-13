from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu, wilcoxon


ROOT = Path("giotto_codex_results")
SOURCE = ROOT / "neutrophil_a5b1_fibroblast_reanalysis"
OUT = ROOT / "a5b1_fibroblast_neutrophil_three_group_analysis"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GROUPS = ("Control", "UC_Noninflamed", "UC_Inflamed")
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
PAIRWISE = (
    ("UC_Inflamed", "Control"),
    ("UC_Noninflamed", "Control"),
    ("UC_Inflamed", "UC_Noninflamed"),
)
PRIMARY_RADIUS = 50
PRIMARY_POPULATION = "All neutrophils"


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


def epsilon_squared_kruskal(h, n, k):
    if n <= k:
        return np.nan
    return max(0.0, (h - k + 1) / (n - k))


def three_group_tests(data, metric, endpoint, allow_zero_tests=False):
    rows = []
    arrays = {
        group: data.loc[data.Diagnosis2.eq(group), metric].dropna().to_numpy()
        for group in GROUPS
    }
    if all(len(arrays[group]) for group in GROUPS):
        test = kruskal(*(arrays[group] for group in GROUPS))
        rows.append(
            {
                "endpoint": endpoint,
                "analysis": "Kruskal-Wallis omnibus",
                "comparison": "Control vs UC noninflamed vs UC inflamed",
                "group_1": "",
                "group_2": "",
                "n_1": len(arrays[GROUPS[0]]),
                "n_2": len(arrays[GROUPS[1]]),
                "n_3": len(arrays[GROUPS[2]]),
                "median_1": np.median(arrays[GROUPS[0]]),
                "median_2": np.median(arrays[GROUPS[1]]),
                "median_3": np.median(arrays[GROUPS[2]]),
                "effect": epsilon_squared_kruskal(test.statistic, sum(map(len, arrays.values())), 3),
                "effect_type": "Kruskal-Wallis epsilon-squared",
                "p_value": test.pvalue,
            }
        )
    for first, second in PAIRWISE:
        x, y = arrays[first], arrays[second]
        if not len(x) or not len(y):
            continue
        test = mannwhitneyu(x, y, alternative="two-sided")
        rows.append(
            {
                "endpoint": endpoint,
                "analysis": "Pairwise Mann-Whitney",
                "comparison": f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}",
                "group_1": first,
                "group_2": second,
                "n_1": len(x),
                "n_2": len(y),
                "n_3": np.nan,
                "median_1": np.median(x),
                "median_2": np.median(y),
                "median_3": np.nan,
                "effect": np.median(x) - np.median(y),
                "effect_type": f"median difference ({GROUP_LABEL[first]} - {GROUP_LABEL[second]})",
                "rank_biserial": rank_biserial(x, y),
                "p_value": test.pvalue,
            }
        )
    if allow_zero_tests:
        for group in GROUPS:
            x = arrays[group]
            if len(x) >= 3 and np.any(x != 0):
                test = wilcoxon(x, alternative="two-sided", zero_method="wilcox")
                rows.append(
                    {
                        "endpoint": endpoint,
                        "analysis": "One-sample Wilcoxon versus zero",
                        "comparison": f"{GROUP_LABEL[group]} vs no preference",
                        "group_1": group,
                        "group_2": "",
                        "n_1": len(x),
                        "n_2": np.nan,
                        "n_3": np.nan,
                        "median_1": np.median(x),
                        "median_2": 0.0,
                        "median_3": np.nan,
                        "effect": np.median(x),
                        "effect_type": "median log2 observed/expected",
                        "rank_biserial": np.nan,
                        "p_value": test.pvalue,
                    }
                )
    return rows


def add_fdr(stats):
    stats = stats.copy()
    stats["fdr_within_endpoint_and_analysis"] = stats.groupby(
        ["endpoint", "analysis"], dropna=False
    ).p_value.transform(bh_adjust)
    return stats


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def require_columns(data, names, label):
    missing = set(names).difference(data.columns)
    if missing:
        raise RuntimeError(f"{label} is missing columns: {sorted(missing)}")


fib = pd.read_csv(SOURCE / "tables" / "01_a5B1_fibroblast_summary_by_patient.csv")
interaction = pd.read_csv(SOURCE / "tables" / "02_spatial_interaction_enrichment_by_patient.csv")
distance = pd.read_csv(SOURCE / "tables" / "04_nearest_distance_summary_by_patient.csv")
sensitivity = pd.read_csv(SOURCE / "tables" / "05_sensitivity_patient_specific_q75.csv")
definition = pd.read_csv(SOURCE / "tables" / "00_a5B1_positive_definition.csv")
cells = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv", low_memory=False)

require_columns(fib, ["PatientID", "Diagnosis2", "fraction_a5B1_positive", "median_a5B1_clr"], "fibroblast summary")
require_columns(
    interaction,
    ["PatientID", "Diagnosis2", "neutrophil_population", "radius_um", "log2_observed_expected"],
    "interaction summary",
)
require_columns(
    distance,
    ["PatientID", "Diagnosis2", "Neutrophil_subtype", "median_nearest_a5B1_positive_um"],
    "distance summary",
)
require_columns(cells, ["PatientID", "Diagnosis2", "cell_type", "a5B1_cell_clr", "x", "y"], "cell export")

observed_groups = set(fib.Diagnosis2.dropna().unique())
if observed_groups != set(GROUPS):
    raise RuntimeError(f"Expected groups {GROUPS}; observed {sorted(observed_groups)}")

# A fixed UC-derived threshold is applied to all groups, so controls are compared
# on the same measurement scale. The within-patient q75 result is a sensitivity analysis.
abundance_rows = []
for metric, endpoint in (
    ("fraction_a5B1_positive", "Fraction of fibroblasts called a5B1-positive"),
    ("median_a5B1_clr", "Median fibroblast a5B1 CLR intensity"),
):
    abundance_rows.extend(three_group_tests(fib, metric, endpoint))
abundance_stats = add_fdr(pd.DataFrame(abundance_rows))
abundance_stats.to_csv(TAB / "01_a5B1_fibroblast_abundance_three_group_statistics.csv", index=False)

# Full interaction analysis at 25, 50, and 100 um for all neutrophils and
# prespecified neutrophil subtypes. FDR is controlled across the three pairwise
# group contrasts separately for every population/radius endpoint.
interaction_rows = []
for (population, radius), q in interaction.groupby(["neutrophil_population", "radius_um"]):
    endpoint = f"{population}; {int(radius)} um log2 observed/expected edges"
    interaction_rows.extend(
        three_group_tests(q, "log2_observed_expected", endpoint, allow_zero_tests=True)
    )
interaction_stats = add_fdr(pd.DataFrame(interaction_rows))
interaction_stats.to_csv(TAB / "02_spatial_interaction_three_group_statistics.csv", index=False)

# Nearest-distance results are secondary because their magnitude is affected by
# the abundance of a5B1-positive fibroblasts; permutation enrichment is primary.
distance_rows = []
for subtype, q in distance.groupby("Neutrophil_subtype"):
    for metric, label in (
        ("median_nearest_a5B1_positive_um", "median nearest a5B1-positive fibroblast distance"),
        ("median_distance_difference_um", "median positive-minus-negative fibroblast distance"),
        ("fraction_closer_to_a5B1_positive", "fraction neutrophils closer to a5B1-positive fibroblasts"),
    ):
        distance_rows.extend(three_group_tests(q, metric, f"{subtype}; {label}"))
distance_stats = add_fdr(pd.DataFrame(distance_rows))
distance_stats.to_csv(TAB / "03_nearest_distance_three_group_statistics.csv", index=False)

# Threshold sensitivity: each patient's top fibroblast quartile is positive.
sensitivity_rows = three_group_tests(
    sensitivity[
        sensitivity.neutrophil_population.eq(PRIMARY_POPULATION)
        & sensitivity.radius_um.eq(PRIMARY_RADIUS)
    ],
    "log2_observed_expected",
    "All neutrophils; 50 um; within-patient fibroblast q75 sensitivity",
    allow_zero_tests=True,
)
sensitivity_stats = add_fdr(pd.DataFrame(sensitivity_rows))
sensitivity_stats.to_csv(TAB / "04_within_patient_q75_sensitivity_three_group_statistics.csv", index=False)

primary = interaction[
    interaction.neutrophil_population.eq(PRIMARY_POPULATION)
    & interaction.radius_um.eq(PRIMARY_RADIUS)
].copy()
primary_summary = (
    primary.groupby("Diagnosis2")
    .agg(
        n_patients=("PatientID", "nunique"),
        median_log2_observed_expected=("log2_observed_expected", "median"),
        q25_log2_observed_expected=("log2_observed_expected", lambda x: x.quantile(0.25)),
        q75_log2_observed_expected=("log2_observed_expected", lambda x: x.quantile(0.75)),
        median_fraction_neutrophils_with_a5B1_positive_fibroblast=(
            "fraction_neutrophils_with_a5B1_positive_fibroblast", "median"
        ),
        median_normalized_rate_ratio=("fibroblast_abundance_normalized_rate_ratio", "median"),
    )
    .reindex(GROUPS)
    .reset_index()
)
primary_summary.to_csv(TAB / "05_primary_50um_group_summary.csv", index=False)

primary_stats = interaction_stats[
    interaction_stats.endpoint.eq("All neutrophils; 50 um log2 observed/expected edges")
].copy()
primary_stats.to_csv(TAB / "06_primary_50um_three_group_statistics.csv", index=False)

secondary_rows = []
for metric, endpoint in (
    (
        "fraction_neutrophils_with_a5B1_positive_fibroblast",
        "All neutrophils; fraction with an a5B1-positive fibroblast within 50 um",
    ),
    (
        "log2_fibroblast_abundance_normalized_rate_ratio",
        "All neutrophils; 50 um log2 a5B1-positive/all-fibroblast edge-rate ratio",
    ),
):
    secondary_rows.extend(three_group_tests(primary, metric, endpoint))
secondary_stats = add_fdr(pd.DataFrame(secondary_rows))
secondary_stats.to_csv(TAB / "07_secondary_50um_three_group_statistics.csv", index=False)

# Representative sample is the patient nearest the group median primary enrichment.
representatives = {}
for group in GROUPS:
    q = primary[primary.Diagnosis2.eq(group)].copy()
    median = q.log2_observed_expected.median()
    representatives[group] = q.loc[(q.log2_observed_expected - median).abs().idxmin(), "PatientID"]

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})
order = [GROUP_LABEL[group] for group in GROUPS]
palette = {GROUP_LABEL[group]: GROUP_COLOR[group] for group in GROUPS}

# Main six-panel quantitative figure.
fig, axes = plt.subplots(2, 3, figsize=(20, 13))

q = fib.assign(Group=fib.Diagnosis2.map(GROUP_LABEL))
sns.boxplot(data=q, x="Group", y="fraction_a5B1_positive", hue="Group", order=order,
            palette=palette, showfliers=False, legend=False, ax=axes[0, 0])
sns.stripplot(data=q, x="Group", y="fraction_a5B1_positive", order=order,
              color="black", size=5, ax=axes[0, 0])
axes[0, 0].set_title("A  a5B1+ fibroblast abundance")
axes[0, 0].set_xlabel("")
axes[0, 0].set_ylabel("Fraction of fibroblasts")

q = primary.assign(Group=primary.Diagnosis2.map(GROUP_LABEL))
sns.boxplot(data=q, x="Group", y="log2_observed_expected", hue="Group", order=order,
            palette=palette, showfliers=False, legend=False, ax=axes[0, 1])
sns.stripplot(data=q, x="Group", y="log2_observed_expected", order=order,
              color="black", size=5, ax=axes[0, 1])
axes[0, 1].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0, 1].set_title("B  Spatial targeting at 50 um")
axes[0, 1].set_xlabel("")
axes[0, 1].set_ylabel("log2 observed / expected edges")

q = distance[distance.Neutrophil_subtype.eq(PRIMARY_POPULATION)].assign(
    Group=lambda x: x.Diagnosis2.map(GROUP_LABEL)
)
sns.boxplot(data=q, x="Group", y="median_nearest_a5B1_positive_um", hue="Group", order=order,
            palette=palette, showfliers=False, legend=False, ax=axes[0, 2])
sns.stripplot(data=q, x="Group", y="median_nearest_a5B1_positive_um", order=order,
              color="black", size=5, ax=axes[0, 2])
axes[0, 2].set_title("C  Nearest a5B1+ fibroblast")
axes[0, 2].set_xlabel("")
axes[0, 2].set_ylabel("Patient median distance (um)")

q = primary.assign(Group=primary.Diagnosis2.map(GROUP_LABEL))
sns.boxplot(data=q, x="Group", y="fraction_neutrophils_with_a5B1_positive_fibroblast",
            hue="Group", order=order, palette=palette, showfliers=False, legend=False,
            ax=axes[1, 0])
sns.stripplot(data=q, x="Group", y="fraction_neutrophils_with_a5B1_positive_fibroblast",
              order=order, color="black", size=5, ax=axes[1, 0])
axes[1, 0].set_title("D  Neutrophils with an a5B1+ neighbor")
axes[1, 0].set_xlabel("")
axes[1, 0].set_ylabel("Fraction within 50 um")

q = interaction[
    interaction.radius_um.eq(PRIMARY_RADIUS)
    & ~interaction.neutrophil_population.eq(PRIMARY_POPULATION)
]
heat = q.groupby(["Diagnosis2", "neutrophil_population"]).log2_observed_expected.median().unstack()
heat = heat.reindex(GROUPS)
sns.heatmap(heat, cmap="vlag", center=0, annot=True, fmt=".2f",
            cbar_kws={"label": "Median log2 O/E"}, ax=axes[1, 1])
axes[1, 1].set_yticklabels(order, rotation=0)
axes[1, 1].set_title("E  Subtype interactions at 50 um")
axes[1, 1].set_xlabel("Neutrophil subtype")
axes[1, 1].set_ylabel("")

sens_primary = sensitivity[
    sensitivity.neutrophil_population.eq(PRIMARY_POPULATION)
    & sensitivity.radius_um.eq(PRIMARY_RADIUS)
][["PatientID", "Diagnosis2", "log2_observed_expected"]].rename(
    columns={"log2_observed_expected": "Within-patient q75"}
)
paired = primary[["PatientID", "Diagnosis2", "log2_observed_expected"]].rename(
    columns={"log2_observed_expected": "Fixed pooled-UC q75"}
).merge(sens_primary, on=["PatientID", "Diagnosis2"])
for _, row in paired.iterrows():
    axes[1, 2].plot(
        [0, 1], [row["Fixed pooled-UC q75"], row["Within-patient q75"]],
        color=GROUP_COLOR[row.Diagnosis2], alpha=0.55, marker="o",
    )
axes[1, 2].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1, 2].set_xticks([0, 1], ["Fixed pooled-UC q75", "Within-patient q75"], rotation=12)
axes[1, 2].set_ylabel("log2 observed / expected edges")
axes[1, 2].set_title("F  Threshold sensitivity")

for ax in axes.flat:
    if ax not in (axes[1, 1], axes[1, 2]):
        ax.tick_params(axis="x", rotation=15)
fig.suptitle("Neutrophil spatial relationships with a5B1+ fibroblasts", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_A5B1_Neutrophil_Three_Group_Quantitative")

# Three representative whole-section maps.
threshold = float(definition.loc[definition.definition.eq("Primary"), "threshold"].iloc[0])
cells["is_fibroblast"] = cells.cell_type.eq("Fibroblast")
cells["is_neutrophil"] = cells.cell_type.str.startswith("Neutrophil", na=False)
cells["is_a5B1_positive_fibroblast"] = cells.is_fibroblast & cells.a5B1_cell_clr.ge(threshold)
fig, axes = plt.subplots(1, 3, figsize=(22, 7))
for ax, group in zip(axes, GROUPS):
    pid = representatives[group]
    z = cells[cells.PatientID.eq(pid)]
    background = z[~z.is_fibroblast & ~z.is_neutrophil]
    negative = z[z.is_fibroblast & ~z.is_a5B1_positive_fibroblast]
    positive = z[z.is_a5B1_positive_fibroblast]
    neutrophils = z[z.is_neutrophil]
    ax.scatter(background.x, background.y, s=0.4, c="#D9D9D9", alpha=0.16, rasterized=True)
    ax.scatter(negative.x, negative.y, s=2, c="#7F7F7F", alpha=0.42,
               label="a5B1-negative fibroblast", rasterized=True)
    ax.scatter(positive.x, positive.y, s=5, c="#C51B7D", alpha=0.82,
               label="a5B1-positive fibroblast", rasterized=True)
    ax.scatter(neutrophils.x, neutrophils.y, s=6, c="#00A6D6", alpha=0.88,
               label="Neutrophil", rasterized=True)
    enrichment = primary.loc[primary.PatientID.eq(pid), "log2_observed_expected"].iloc[0]
    ax.set_title(f"{GROUP_LABEL[group]}: {pid}\n50 um log2 O/E = {enrichment:.2f}")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
axes[-1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
fig.suptitle("Representative a5B1 fibroblast-neutrophil spatial organization", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_A5B1_Neutrophil_Three_Group_Spatial_Maps")

def fmt_p(value):
    return f"{value:.3g}"


omnibus = primary_stats[primary_stats.analysis.eq("Kruskal-Wallis omnibus")].iloc[0]
pairwise = primary_stats[primary_stats.analysis.eq("Pairwise Mann-Whitney")].set_index("comparison")
one_sample = primary_stats[
    primary_stats.analysis.eq("One-sample Wilcoxon versus zero")
].set_index("group_1")
abundance_primary = abundance_stats[
    abundance_stats.endpoint.eq("Fraction of fibroblasts called a5B1-positive")
]
abundance_omnibus = abundance_primary[
    abundance_primary.analysis.eq("Kruskal-Wallis omnibus")
].iloc[0]
abundance_pairwise = abundance_primary[
    abundance_primary.analysis.eq("Pairwise Mann-Whitney")
].set_index("comparison")
sens_omnibus = sensitivity_stats[
    sensitivity_stats.analysis.eq("Kruskal-Wallis omnibus")
].iloc[0]

all_neut_distance = distance[
    distance.Neutrophil_subtype.eq(PRIMARY_POPULATION)
].copy()
distance_medians = all_neut_distance.groupby("Diagnosis2").median_nearest_a5B1_positive_um.median()
distance_endpoint = "All neutrophils; median nearest a5B1-positive fibroblast distance"
distance_primary_stats = distance_stats[distance_stats.endpoint.eq(distance_endpoint)]
distance_omnibus = distance_primary_stats[
    distance_primary_stats.analysis.eq("Kruskal-Wallis omnibus")
].iloc[0]
distance_pairwise = distance_primary_stats[
    distance_primary_stats.analysis.eq("Pairwise Mann-Whitney")
].set_index("comparison")

subtype_within = interaction_stats[
    interaction_stats.endpoint.str.contains("; 50 um log2 observed/expected edges", regex=False)
    & ~interaction_stats.endpoint.str.startswith(PRIMARY_POPULATION)
    & interaction_stats.analysis.eq("One-sample Wilcoxon versus zero")
    & interaction_stats.effect.lt(0)
    & interaction_stats.fdr_within_endpoint_and_analysis.lt(0.05)
].copy()
subtype_between = interaction_stats[
    interaction_stats.endpoint.str.contains("; 50 um log2 observed/expected edges", regex=False)
    & ~interaction_stats.endpoint.str.startswith(PRIMARY_POPULATION)
    & interaction_stats.analysis.eq("Pairwise Mann-Whitney")
    & interaction_stats.fdr_within_endpoint_and_analysis.lt(0.05)
]

group_summary = primary_summary.set_index("Diagnosis2")
fib_group = fib.groupby("Diagnosis2").fraction_a5B1_positive.median()

pair_lines = []
for first, second in PAIRWISE:
    label = f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}"
    row = pairwise.loc[label]
    pair_lines.append(
        f"- {label}: median difference {row.effect:.3f}, rank-biserial {row.rank_biserial:.3f}, "
        f"P={fmt_p(row.p_value)}, FDR={fmt_p(row.fdr_within_endpoint_and_analysis)}."
    )

within_lines = []
for group in GROUPS:
    row = one_sample.loc[group]
    direction = "enrichment" if row.effect > 0 else "depletion"
    within_lines.append(
        f"- {GROUP_LABEL[group]}: median {row.effect:.3f} ({direction}), "
        f"P={fmt_p(row.p_value)}, FDR={fmt_p(row.fdr_within_endpoint_and_analysis)}."
    )

abundance_lines = []
for first, second in PAIRWISE:
    label = f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}"
    row = abundance_pairwise.loc[label]
    abundance_lines.append(
        f"- {label}: median difference {row.effect:.3f}, P={fmt_p(row.p_value)}, "
        f"FDR={fmt_p(row.fdr_within_endpoint_and_analysis)}."
    )

distance_lines = []
for first, second in PAIRWISE:
    label = f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}"
    row = distance_pairwise.loc[label]
    distance_lines.append(
        f"- {label}: median difference {row.effect:.1f} um, P={fmt_p(row.p_value)}, "
        f"FDR={fmt_p(row.fdr_within_endpoint_and_analysis)}."
    )

subtype_lines = []
for group in GROUPS:
    hit = subtype_within[subtype_within.group_1.eq(group)]
    names = [endpoint.split(";")[0] for endpoint in hit.endpoint]
    subtype_lines.append(
        f"- {GROUP_LABEL[group]}: {', '.join(names) if names else 'none'}."
    )

report = f"""# a5B1-positive fibroblast-neutrophil spatial analysis: three-group comparison

## Design

The analysis included {int(group_summary.loc['Control', 'n_patients'])} controls, {int(group_summary.loc['UC_Noninflamed', 'n_patients'])} UC noninflamed samples, and {int(group_summary.loc['UC_Inflamed', 'n_patients'])} UC inflamed samples. Patient was the biological replicate. Fibroblasts were called a5B1-positive using the fixed pooled-UC fibroblast 75th-percentile CLR threshold ({threshold:.4f}) for all three groups. The primary endpoint was the within-patient log2 ratio of observed neutrophil-a5B1-positive fibroblast edges within 50 um to the expectation from 1,000 random permutations of a5B1 labels among that patient's fibroblasts.

## Primary 50 um spatial result

Median log2 observed/expected values were {group_summary.loc['Control', 'median_log2_observed_expected']:.3f} in controls, {group_summary.loc['UC_Noninflamed', 'median_log2_observed_expected']:.3f} in UC noninflamed, and {group_summary.loc['UC_Inflamed', 'median_log2_observed_expected']:.3f} in UC inflamed tissue. The omnibus three-group Kruskal-Wallis test was P={fmt_p(omnibus.p_value)} (epsilon-squared={omnibus.effect:.3f}). Pairwise patient-level results, with FDR correction across the three contrasts, were:

{chr(10).join(pair_lines)}

Within-group tests against no spatial preference (log2 observed/expected = 0), FDR-corrected across the three groups, were:

{chr(10).join(within_lines)}

## a5B1-positive fibroblast abundance

Median positive fractions were {fib_group['Control']:.3f} in controls, {fib_group['UC_Noninflamed']:.3f} in UC noninflamed, and {fib_group['UC_Inflamed']:.3f} in UC inflamed tissue. The abundance omnibus test was P={fmt_p(abundance_omnibus.p_value)} (epsilon-squared={abundance_omnibus.effect:.3f}). Pairwise results were:

{chr(10).join(abundance_lines)}

## Nearest-distance and neutrophil-subtype results

The patient-level median distance from a neutrophil to the nearest a5B1-positive fibroblast was {distance_medians['Control']:.1f} um in controls, {distance_medians['UC_Noninflamed']:.1f} um in UC noninflamed, and {distance_medians['UC_Inflamed']:.1f} um in UC inflamed tissue (omnibus P={fmt_p(distance_omnibus.p_value)}). Pairwise results were:

{chr(10).join(distance_lines)}

Because nearest-positive distance depends on positive-fibroblast abundance, it is a secondary endpoint; the label-permutation enrichment above is primary. At 50 um, neutrophil subtypes significantly depleted relative to no spatial preference (FDR < 0.05 within each subtype endpoint across groups) were:

{chr(10).join(subtype_lines)}

No subtype-level pairwise group contrast had FDR < 0.05 ({len(subtype_between)} significant contrasts).

## Sensitivity and interpretation

With a within-patient top-quartile fibroblast definition, the three-group spatial-enrichment omnibus test was P={fmt_p(sens_omnibus.p_value)}. Complete 25, 50, and 100 um results, neutrophil-subtype analyses, distance endpoints, exact P values, effect sizes, and FDR values are in the accompanying tables.

The primary permutation endpoint adjusts for each patient's fibroblast positions and number of a5B1-positive fibroblasts. It measures preferential spatial association, not receptor binding, signaling direction, or causality. The a5B1-positive call is an operational high-expression threshold because no external positivity control was supplied.
"""
(OUT / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

manifest = {
    "source_analysis": str(SOURCE),
    "cell_input": str(ROOT / "additional_analyses" / "codex_extended_cells.csv"),
    "groups": list(GROUPS),
    "patients_per_group": {
        group: int(primary.loc[primary.Diagnosis2.eq(group), "PatientID"].nunique())
        for group in GROUPS
    },
    "biological_replicate": "patient",
    "primary_endpoint": "log2 observed/permuted-expected neutrophil-a5B1-positive fibroblast edges",
    "primary_radius_um": PRIMARY_RADIUS,
    "a5B1_threshold_clr": threshold,
    "a5B1_threshold_source": "pooled UC fibroblast 75th percentile",
    "pairwise_multiplicity": "Benjamini-Hochberg across three prespecified group contrasts per endpoint",
    "representative_patients": representatives,
}
(OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print("Completed three-group a5B1 fibroblast-neutrophil spatial analysis.")
print(f"Saved outputs to {OUT}")
