from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore")

ROOT = Path("giotto_codex_results")
SOURCE = ROOT / "fap_a5b1_neutrophil_reprogramming"
NORMALIZED = SOURCE / "cell_count_normalized"
OUT = ROOT / "uc_inflamed_noninflamed_fap_a5b1_comparison"
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
PRIMARY_RADIUS = 50
N_BOOT = 200
N_PROXIMAL = 5
N_DISTANT = 10
SEED = 20260804
MARKERS = (
    "OSM", "CXCR4", "PDL1", "HLA_ABC", "a5B1", "PADI4", "MX1",
    "CD66b", "CD16", "CD11b", "Ki67",
)
PROGRAM_MARKERS = ("OSM", "CXCR4", "PDL1", "HLA_ABC", "a5B1")
PROGRAM_NAME = "NF-niche program"


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


def uc_tests(frame, value, endpoint):
    noninflamed = frame.loc[frame.Diagnosis2.eq("UC_Noninflamed"), value].dropna().to_numpy()
    inflamed = frame.loc[frame.Diagnosis2.eq("UC_Inflamed"), value].dropna().to_numpy()
    rows = []
    if len(noninflamed) and len(inflamed):
        test = mannwhitneyu(inflamed, noninflamed, alternative="two-sided")
        rows.append({
            "endpoint": endpoint, "analysis": "UC inflamed versus UC noninflamed",
            "group": "", "n_UC_noninflamed": len(noninflamed),
            "n_UC_inflamed": len(inflamed),
            "median_UC_noninflamed": np.median(noninflamed),
            "median_UC_inflamed": np.median(inflamed),
            "effect_inflamed_minus_noninflamed": np.median(inflamed) - np.median(noninflamed),
            "rank_biserial": rank_biserial(inflamed, noninflamed),
            "p_value": test.pvalue,
        })
    for group, values in (("UC_Noninflamed", noninflamed), ("UC_Inflamed", inflamed)):
        if len(values) >= 3 and np.any(values != 0):
            test = wilcoxon(values, alternative="two-sided", zero_method="wilcox")
            rows.append({
                "endpoint": endpoint, "analysis": "Enrichment/effect versus zero",
                "group": group, "n_UC_noninflamed": len(noninflamed),
                "n_UC_inflamed": len(inflamed),
                "median_UC_noninflamed": np.median(noninflamed),
                "median_UC_inflamed": np.median(inflamed),
                "effect_inflamed_minus_noninflamed": np.median(values),
                "rank_biserial": np.nan, "p_value": test.pvalue,
            })
    return rows


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


# Aggregate the mutually exclusive states into marginal FAP+ and a5B1+ targets.
state_spatial = pd.read_csv(SOURCE / "tables" / "03_spatial_enrichment_by_patient.csv")
state_spatial = state_spatial[state_spatial.Diagnosis2.isin(GROUPS)]
aggregate_rows = []
for target, states in TARGETS.items():
    q = state_spatial[state_spatial.fibroblast_state.isin(states)]
    a = q.groupby(["PatientID", "Diagnosis2", "radius_um"]).agg(
        n_neutrophils=("n_neutrophils", "first"),
        n_fibroblasts=("n_fibroblasts", "first"),
        n_target_fibroblasts=("n_state_fibroblasts", "sum"),
        total_edges=("total_edges", "first"),
        observed_target_edges=("observed_state_edges", "sum"),
        expected_target_edges=("expected_state_edges", "sum"),
    ).reset_index()
    a["target"] = target
    a["log2_observed_expected"] = np.log2(
        (a.observed_target_edges + .5) / (a.expected_target_edges + .5)
    )
    a["edges_per_million_possible_pairs"] = (
        a.observed_target_edges / (a.n_neutrophils * a.n_target_fibroblasts) * 1e6
    )
    aggregate_rows.append(a)
spatial = pd.concat(aggregate_rows, ignore_index=True)
spatial.to_csv(TAB / "01_multiscale_spatial_metrics_by_patient.csv", index=False)

spatial_stats_rows = []
for (target, radius), q in spatial.groupby(["target", "radius_um"]):
    spatial_stats_rows.extend(uc_tests(
        q, "log2_observed_expected", f"{target}; {int(radius)} um log2 observed/expected"
    ))
spatial_stats = pd.DataFrame(spatial_stats_rows)
spatial_stats["fdr_within_analysis"] = spatial_stats.groupby("analysis").p_value.transform(bh_adjust)
spatial_stats.to_csv(TAB / "02_multiscale_spatial_statistics.csv", index=False)

# Equal-count bootstrap aggregation: exactly 250 fibroblasts and 30 neutrophils
# were used in every source iteration.
boot = pd.read_csv(NORMALIZED / "tables" / "03_equal_count_spatial_bootstrap_iterations.csv")
boot = boot[boot.Diagnosis2.isin(GROUPS)]
normalized_rows = []
for target, states in TARGETS.items():
    q = boot[boot.fibroblast_state.isin(states)]
    iteration = q.groupby(["PatientID", "Diagnosis2", "iteration"]).agg(
        sampled_target_fibroblasts=("sampled_state_fibroblasts", "sum"),
        observed_target_edges=("state_edges", "sum"),
        expected_target_edges=("expected_state_edges", "sum"),
        total_edges=("total_edges", "first"),
    ).reset_index()
    patient = iteration.groupby(["PatientID", "Diagnosis2"]).agg(
        n_iterations=("iteration", "nunique"),
        median_sampled_target_fibroblasts=("sampled_target_fibroblasts", "median"),
        summed_observed_target_edges=("observed_target_edges", "sum"),
        summed_expected_target_edges=("expected_target_edges", "sum"),
    ).reset_index()
    patient["target"] = target
    patient["aggregate_log2_observed_expected"] = np.log2(
        (patient.summed_observed_target_edges + .5)
        / (patient.summed_expected_target_edges + .5)
    )
    normalized_rows.append(patient)
normalized = pd.concat(normalized_rows, ignore_index=True)
normalized.to_csv(TAB / "03_equal_count_spatial_metrics_by_patient.csv", index=False)

normalized_stats_rows = []
for target, q in normalized.groupby("target"):
    normalized_stats_rows.extend(uc_tests(
        q, "aggregate_log2_observed_expected",
        f"{target}; equal 250-fibroblast/30-neutrophil aggregate log2 O/E",
    ))
normalized_stats = pd.DataFrame(normalized_stats_rows)
normalized_stats["fdr_within_analysis"] = normalized_stats.groupby("analysis").p_value.transform(bh_adjust)
normalized_stats.to_csv(TAB / "04_equal_count_spatial_statistics.csv", index=False)

# Paired within-patient comparison of FAP+ versus a5B1+ targeting.
paired = normalized[normalized.target.isin(("FAP+ fibroblasts", "a5B1+ fibroblasts"))].pivot(
    index=["PatientID", "Diagnosis2"], columns="target", values="aggregate_log2_observed_expected"
).reset_index()
paired["FAP_minus_a5B1_log2_enrichment"] = paired["FAP+ fibroblasts"] - paired["a5B1+ fibroblasts"]
paired.to_csv(TAB / "05_paired_FAP_vs_a5B1_targeting_by_patient.csv", index=False)
paired_rows = []
for group in GROUPS:
    x = paired.loc[paired.Diagnosis2.eq(group), "FAP_minus_a5B1_log2_enrichment"].dropna().to_numpy()
    test = wilcoxon(x, alternative="two-sided", zero_method="wilcox")
    paired_rows.append({
        "analysis": "Within-patient FAP+ versus a5B1+ targeting",
        "group": group, "n_patients": len(x), "median_FAP_minus_a5B1": np.median(x),
        "p_value": test.pvalue,
    })
x = paired.loc[paired.Diagnosis2.eq("UC_Inflamed"), "FAP_minus_a5B1_log2_enrichment"].to_numpy()
y = paired.loc[paired.Diagnosis2.eq("UC_Noninflamed"), "FAP_minus_a5B1_log2_enrichment"].to_numpy()
test = mannwhitneyu(x, y, alternative="two-sided")
paired_rows.append({
    "analysis": "Inflammation difference in FAP-versus-a5B1 preference",
    "group": "UC inflamed vs UC noninflamed", "n_patients": len(x) + len(y),
    "median_FAP_minus_a5B1": np.median(x) - np.median(y), "p_value": test.pvalue,
})
paired_stats = pd.DataFrame(paired_rows)
paired_stats["fdr"] = bh_adjust(paired_stats.p_value)
paired_stats.to_csv(TAB / "06_paired_FAP_vs_a5B1_targeting_statistics.csv", index=False)

# Marginal target abundance per patient.
state_abundance = pd.read_csv(SOURCE / "tables" / "01_fibroblast_state_abundance_by_patient.csv")
state_abundance = state_abundance[state_abundance.Diagnosis2.isin(GROUPS)]
abundance_rows = []
for target, states in TARGETS.items():
    q = state_abundance[state_abundance.fibroblast_state.isin(states)].groupby(
        ["PatientID", "Diagnosis2"]
    ).agg(n_fibroblasts=("n_fibroblasts", "first"), n_target=("n_state", "sum")).reset_index()
    q["target"] = target
    q["fraction_fibroblasts"] = q.n_target / q.n_fibroblasts
    abundance_rows.append(q)
abundance = pd.concat(abundance_rows, ignore_index=True)
abundance.to_csv(TAB / "07_target_abundance_by_patient.csv", index=False)
abundance_stats_rows = []
for target, q in abundance.groupby("target"):
    abundance_stats_rows.extend(uc_tests(q, "fraction_fibroblasts", target))
abundance_stats = pd.DataFrame(abundance_stats_rows)
abundance_stats["fdr_within_analysis"] = abundance_stats.groupby("analysis").p_value.transform(bh_adjust)
abundance_stats.to_csv(TAB / "08_target_abundance_statistics.csv", index=False)

# Equal near/far neutrophil counts for marginal-target marker comparisons.
cells = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv", low_memory=False)
crc = pd.read_csv(ROOT / "neutrophil_fibroblast_crc_marker_map" / "crc_associated_marker_expression.csv")
d = cells.merge(crc, on="cell_ID", how="left", validate="one_to_one")
d["is_neutrophil"] = d.cell_type.str.startswith("Neutrophil", na=False)
dist = pd.read_csv(SOURCE / "tables" / "04_neutrophil_fibroblast_state_distances.csv")
wide = dist.pivot(index="cell_ID", columns="fibroblast_state", values="nearest_state_distance_um")
neut = d[d.is_neutrophil].merge(wide, left_on="cell_ID", right_index=True, how="left", validate="one_to_one")
neut = neut[neut.Diagnosis2.isin(GROUPS)].copy()
for target, states in TARGETS.items():
    neut[target] = neut[list(states)].min(axis=1)
marker_source = {
    "OSM": "raw_OSM", "CXCR4": "raw_CXCR4", "PDL1": "raw_PDL1",
    "HLA_ABC": "raw_HLA_ABC", "a5B1": "raw_a5B1", "PADI4": "PADI4",
    "MX1": "MX.1", "CD66b": "CD66b", "CD16": "CD16", "CD11b": "CD11b",
    "Ki67": "raw_Ki67",
}
for marker, source in marker_source.items():
    neut[f"log2_{marker}"] = np.log2(neut[source].astype(float).clip(lower=0) + .1)
for marker in PROGRAM_MARKERS:
    col = f"log2_{marker}"
    neut[f"z_{marker}"] = neut.groupby("PatientID")[col].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0.0
    )
neut[PROGRAM_NAME] = neut[[f"z_{m}" for m in PROGRAM_MARKERS]].mean(axis=1)

rng = np.random.default_rng(SEED)
features = list(MARKERS) + [PROGRAM_NAME]
marker_rows = []
for (pid, diagnosis), q in neut.groupby(["PatientID", "Diagnosis2"]):
    for target in TARGETS:
        proximal = q[q[target].le(50)]
        distant = q[q[target].gt(100)]
        if len(proximal) < N_PROXIMAL or len(distant) < N_DISTANT:
            continue
        pi = rng.integers(0, len(proximal), size=(N_BOOT, N_PROXIMAL))
        di = rng.integers(0, len(distant), size=(N_BOOT, N_DISTANT))
        for feature in features:
            p = (proximal[PROGRAM_NAME] if feature == PROGRAM_NAME else proximal[f"log2_{feature}"]).to_numpy()
            f = (distant[PROGRAM_NAME] if feature == PROGRAM_NAME else distant[f"log2_{feature}"]).to_numpy()
            effects = np.median(p[pi], axis=1) - np.median(f[di], axis=1)
            marker_rows.append({
                "PatientID": pid, "Diagnosis2": diagnosis, "target": target,
                "feature": feature, "n_available_proximal": len(proximal),
                "n_available_distant": len(distant),
                "bootstrap_median_effect": np.median(effects),
                "bootstrap_ci025": np.quantile(effects, .025),
                "bootstrap_ci975": np.quantile(effects, .975),
            })
marker_effects = pd.DataFrame(marker_rows)
marker_effects.to_csv(TAB / "09_equal_count_marker_effects_by_patient.csv", index=False)

marker_stats_rows = []
for (target, feature), q in marker_effects.groupby(["target", "feature"]):
    marker_stats_rows.extend(uc_tests(
        q, "bootstrap_median_effect", f"{target}; {feature}; equal 5-near/10-far effect"
    ))
marker_stats = pd.DataFrame(marker_stats_rows)
marker_stats["target"] = marker_stats.endpoint.str.split(";").str[0]
marker_stats["feature"] = marker_stats.endpoint.str.split(";").str[1].str.strip()
marker_stats["fdr_within_target_analysis"] = marker_stats.groupby(
    ["target", "analysis"]
).p_value.transform(bh_adjust)
marker_stats.to_csv(TAB / "10_equal_count_marker_statistics.csv", index=False)

# Figures.
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})
order = [GROUP_LABEL[g] for g in GROUPS]
palette = {GROUP_LABEL[g]: GROUP_COLOR[g] for g in GROUPS}

fig, axes = plt.subplots(2, 3, figsize=(21, 13))

q = abundance.copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="target", y="fraction_fibroblasts", hue="Group", hue_order=order,
            palette=palette, showfliers=False, ax=axes[0, 0])
axes[0, 0].set_title("A  Fibroblast target abundance")
axes[0, 0].set_xlabel("")
axes[0, 0].set_ylabel("Fraction of fibroblasts")
axes[0, 0].tick_params(axis="x", rotation=25, labelsize=9)
axes[0, 0].legend(title="", frameon=False, fontsize=9)

multi = spatial.groupby(["target", "radius_um", "Diagnosis2"]).log2_observed_expected.median().unstack("Diagnosis2")
multi["inflamed_minus_noninflamed"] = multi["UC_Inflamed"] - multi["UC_Noninflamed"]
multi_heat = multi.inflamed_minus_noninflamed.unstack("radius_um").reindex(index=TARGETS)
sns.heatmap(multi_heat, annot=True, fmt=".2f", cmap="vlag", center=0,
            cbar_kws={"label": "Inflamed - noninflamed log2 O/E"}, ax=axes[0, 1])
axes[0, 1].set_title("B  Multiscale UC-state difference")
axes[0, 1].set_xlabel("Radius (um)")
axes[0, 1].set_ylabel("")
axes[0, 1].tick_params(axis="y", labelsize=9, rotation=0)

q = normalized.copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="target", y="aggregate_log2_observed_expected", hue="Group",
            hue_order=order, palette=palette, showfliers=False, ax=axes[0, 2])
axes[0, 2].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0, 2].set_title("C  Equal-count targeting at 50 um")
axes[0, 2].set_xlabel("")
axes[0, 2].set_ylabel("Aggregate log2 observed / expected")
axes[0, 2].tick_params(axis="x", rotation=25, labelsize=9)
axes[0, 2].legend(title="", frameon=False, fontsize=9)

for _, row in paired.iterrows():
    axes[1, 0].plot(
        [0, 1], [row["FAP+ fibroblasts"], row["a5B1+ fibroblasts"]],
        color=GROUP_COLOR[row.Diagnosis2], alpha=.55, marker="o",
    )
axes[1, 0].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1, 0].set_xticks([0, 1], ["FAP+", "a5B1+"])
axes[1, 0].set_ylabel("Equal-count log2 observed / expected")
axes[1, 0].set_title("D  Within-patient target preference")

heat = marker_effects.groupby(["Diagnosis2", "target", "feature"]).bootstrap_median_effect.median().reset_index()
heat["row"] = heat.Diagnosis2.map(GROUP_LABEL) + " | " + heat.target
heat = heat.pivot(index="row", columns="feature", values="bootstrap_median_effect")
row_order = [f"{GROUP_LABEL[g]} | {t}" for g in GROUPS for t in TARGETS]
heat = heat.reindex(index=row_order, columns=features)
sns.heatmap(heat, cmap="vlag", center=0, cbar_kws={"label": "Median near-far effect"}, ax=axes[1, 1])
axes[1, 1].set_title("E  Equal-count neutrophil phenotype")
axes[1, 1].set_xlabel("")
axes[1, 1].set_ylabel("")
axes[1, 1].tick_params(axis="x", rotation=45, labelsize=8)
axes[1, 1].tick_params(axis="y", labelsize=8)

diff = marker_effects.groupby(["Diagnosis2", "target", "feature"]).bootstrap_median_effect.median().unstack("Diagnosis2")
diff["inflamed_minus_noninflamed"] = diff["UC_Inflamed"] - diff["UC_Noninflamed"]
diff_heat = diff["inflamed_minus_noninflamed"].unstack("feature").reindex(index=TARGETS, columns=features)
sns.heatmap(diff_heat, cmap="vlag", center=0, cbar_kws={"label": "Inflamed - noninflamed"}, ax=axes[1, 2])
axes[1, 2].set_title("F  Inflammation-associated marker difference")
axes[1, 2].set_xlabel("")
axes[1, 2].set_ylabel("")
axes[1, 2].tick_params(axis="x", rotation=45, labelsize=8)
axes[1, 2].tick_params(axis="y", labelsize=9)

fig.suptitle("Neutrophil interactions with FAP+ and a5B1+ fibroblasts in UC", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_UC_Inflamed_vs_Noninflamed_FAP_a5B1_Comparison")


def fmt(x):
    return f"{x:.3g}"


primary_between = normalized_stats[
    normalized_stats.analysis.eq("UC inflamed versus UC noninflamed")
].set_index(normalized_stats[normalized_stats.analysis.eq("UC inflamed versus UC noninflamed")].endpoint.str.split(";").str[0])
primary_within = normalized_stats[
    normalized_stats.analysis.eq("Enrichment/effect versus zero")
]
primary_medians = normalized.groupby(["target", "Diagnosis2"]).aggregate_log2_observed_expected.median()

spatial_lines = []
for target in TARGETS:
    row = primary_between.loc[target]
    spatial_lines.append(
        f"- {target}: UC noninflamed {primary_medians.loc[(target, 'UC_Noninflamed')]:.3f}; "
        f"UC inflamed {primary_medians.loc[(target, 'UC_Inflamed')]:.3f}; "
        f"inflamed-minus-noninflamed {row.effect_inflamed_minus_noninflamed:.3f}, "
        f"P={fmt(row.p_value)}, FDR={fmt(row.fdr_within_analysis)}."
    )

within_spatial_lines = []
for target in TARGETS:
    for group in GROUPS:
        row = primary_within[
            primary_within.endpoint.str.startswith(target + ";")
            & primary_within.group.eq(group)
        ].iloc[0]
        within_spatial_lines.append(
            f"- {target}, {GROUP_LABEL[group]}: median {row.effect_inflamed_minus_noninflamed:.3f}, "
            f"P={fmt(row.p_value)}, FDR={fmt(row.fdr_within_analysis)}."
        )

paired_text = []
for row in paired_stats.itertuples():
    paired_text.append(
        f"- {row.group}: median effect {row.median_FAP_minus_a5B1:.3f}, "
        f"P={fmt(row.p_value)}, FDR={fmt(row.fdr)}."
    )

marker_between = marker_stats[
    marker_stats.analysis.eq("UC inflamed versus UC noninflamed")
    & marker_stats.fdr_within_target_analysis.lt(.05)
]
marker_lines = []
for target in TARGETS:
    q = marker_between[marker_between.target.eq(target)].sort_values(
        "effect_inflamed_minus_noninflamed", ascending=False
    )
    hits = [f"{r.feature} ({r.effect_inflamed_minus_noninflamed:+.2f})" for r in q.itertuples()]
    marker_lines.append(f"- {target}: {', '.join(hits) if hits else 'no marker difference after FDR correction'}.")

marker_within = marker_stats[
    marker_stats.analysis.eq("Enrichment/effect versus zero")
    & marker_stats.fdr_within_target_analysis.lt(.05)
]
marker_within_lines = []
for target in TARGETS:
    for group in GROUPS:
        q = marker_within[
            marker_within.target.eq(target) & marker_within.group.eq(group)
        ].sort_values("effect_inflamed_minus_noninflamed", ascending=False)
        hits = [f"{r.feature} ({r.effect_inflamed_minus_noninflamed:+.2f})" for r in q.itertuples()]
        marker_within_lines.append(
            f"- {target}, {GROUP_LABEL[group]}: "
            f"{', '.join(hits) if hits else 'none after FDR correction'}."
        )

abundance_between = abundance_stats[
    abundance_stats.analysis.eq("UC inflamed versus UC noninflamed")
].set_index("endpoint")
abundance_lines = []
for target in TARGETS:
    row = abundance_between.loc[target]
    abundance_lines.append(
        f"- {target}: UC noninflamed {row.median_UC_noninflamed:.3f}; "
        f"UC inflamed {row.median_UC_inflamed:.3f}; P={fmt(row.p_value)}, "
        f"FDR={fmt(row.fdr_within_analysis)}."
    )

report = f"""# Neutrophil interactions with FAP+ and a5B1+ fibroblasts: UC inflamed versus noninflamed

## Definitions and normalization

FAP+ includes FAP+/a5B1+ and FAP+/a5B1- fibroblasts; a5B1+ includes FAP+/a5B1+ and FAP-/a5B1+ fibroblasts. The overlapping double-positive subset is also reported separately. Primary comparisons use the cell-count-normalized bootstrap with exactly 250 fibroblasts and 30 neutrophils per iteration. P21 was excluded because it contained only 10 neutrophils. Patient was the biological replicate.

## Equal-count 50 um spatial comparison

{chr(10).join(spatial_lines)}

Within-group tests against no spatial preference were:

{chr(10).join(within_spatial_lines)}

Within-patient FAP+ versus a5B1+ preference tests were:

{chr(10).join(paired_text)}

## Neutrophil phenotype differences between UC states

Markers with a significantly different near-versus-far effect in inflamed versus noninflamed UC were:

{chr(10).join(marker_lines)}

Significant marker effects within each UC state were:

{chr(10).join(marker_within_lines)}

## Fibroblast target abundance

{chr(10).join(abundance_lines)}

Complete 25, 50, and 100 um results and within-group tests against no spatial preference are provided in the tables. Spatial and marker differences are associations and do not establish causal fibroblast-to-neutrophil reprogramming.
"""
(OUT / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

manifest = {
    "targets": TARGETS, "primary_radius_um": PRIMARY_RADIUS,
    "primary_normalization": "200 iterations of 250 fibroblasts and 30 neutrophils per sample",
    "marker_normalization": "200 iterations of 5 proximal and 10 distant neutrophils",
    "excluded_primary_sample": "P21", "biological_replicate": "patient",
    "seed": SEED,
}
(OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print("Completed UC inflamed versus noninflamed FAP/a5B1 interaction comparison.")
print(f"Saved outputs to {OUT}")
