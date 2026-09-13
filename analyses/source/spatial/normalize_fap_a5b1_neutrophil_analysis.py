from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import kruskal, mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore")

ROOT = Path("giotto_codex_results")
SOURCE = ROOT / "fap_a5b1_neutrophil_reprogramming"
OUT = SOURCE / "cell_count_normalized"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GROUPS = ("Control", "UC_Noninflamed", "UC_Inflamed")
GROUP_LABEL = {
    "Control": "Control", "UC_Noninflamed": "UC noninflamed",
    "UC_Inflamed": "UC inflamed",
}
GROUP_COLOR = {
    "Control": "#4C78A8", "UC_Noninflamed": "#D4AD35",
    "UC_Inflamed": "#D95F5F",
}
STATES = ("FAP+/a5B1+", "FAP+/a5B1-", "FAP-/a5B1+", "FAP-/a5B1-")
PAIRWISE = (
    ("UC_Inflamed", "Control"),
    ("UC_Noninflamed", "Control"),
    ("UC_Inflamed", "UC_Noninflamed"),
)
RADIUS = 50
N_FIBROBLASTS = 250
N_NEUTROPHILS = 30
N_PROXIMAL = 5
N_DISTANT = 10
N_BOOT = 200
MIN_VALID_BOOTSTRAPS = 20
SEED = 20260803

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


def group_tests(frame, value, endpoint, zero_tests=True):
    arrays = {
        group: frame.loc[frame.Diagnosis2.eq(group), value].dropna().to_numpy()
        for group in GROUPS
    }
    rows = []
    if all(len(arrays[group]) >= 3 for group in GROUPS):
        test = kruskal(*(arrays[group] for group in GROUPS))
        rows.append({
            "endpoint": endpoint, "analysis": "Kruskal-Wallis omnibus",
            "comparison": "Control vs UC noninflamed vs UC inflamed",
            "group_1": "", "group_2": "",
            "n_1": len(arrays[GROUPS[0]]), "n_2": len(arrays[GROUPS[1]]),
            "n_3": len(arrays[GROUPS[2]]), "effect": np.nan,
            "p_value": test.pvalue,
        })
    for first, second in PAIRWISE:
        x, y = arrays[first], arrays[second]
        if not len(x) or not len(y):
            continue
        test = mannwhitneyu(x, y, alternative="two-sided")
        rows.append({
            "endpoint": endpoint, "analysis": "Pairwise Mann-Whitney",
            "comparison": f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}",
            "group_1": first, "group_2": second,
            "n_1": len(x), "n_2": len(y), "n_3": np.nan,
            "effect": np.median(x) - np.median(y), "p_value": test.pvalue,
            "rank_biserial": rank_biserial(x, y),
        })
    if zero_tests:
        for group in GROUPS:
            x = arrays[group]
            if len(x) >= 3 and np.any(x != 0):
                test = wilcoxon(x, alternative="two-sided", zero_method="wilcox")
                rows.append({
                    "endpoint": endpoint, "analysis": "One-sample Wilcoxon versus zero",
                    "comparison": f"{GROUP_LABEL[group]} vs zero",
                    "group_1": group, "group_2": "", "n_1": len(x),
                    "n_2": np.nan, "n_3": np.nan, "effect": np.median(x),
                    "p_value": test.pvalue, "rank_biserial": np.nan,
                })
    return rows


def feature_statistics(frame):
    rows = []
    for (state, feature), q in frame.groupby(["fibroblast_state", "feature"]):
        arrays = {}
        for group in GROUPS:
            x = q.loc[q.Diagnosis2.eq(group), "bootstrap_median_effect"].dropna().to_numpy()
            arrays[group] = x
            if len(x) >= 3 and np.any(x != 0):
                test = wilcoxon(x, alternative="two-sided", zero_method="wilcox")
                rows.append({
                    "fibroblast_state": state, "feature": feature,
                    "analysis": "One-sample Wilcoxon versus zero",
                    "stratum": GROUP_LABEL[group], "n_patients": len(x),
                    "median_effect": np.median(x), "p_value": test.pvalue,
                })
        if all(len(arrays[group]) >= 3 for group in GROUPS):
            test = kruskal(*(arrays[group] for group in GROUPS))
            rows.append({
                "fibroblast_state": state, "feature": feature,
                "analysis": "Kruskal-Wallis omnibus",
                "stratum": "Control vs UC noninflamed vs UC inflamed",
                "n_patients": sum(map(len, arrays.values())),
                "median_effect": np.nan, "p_value": test.pvalue,
            })
        for first, second in PAIRWISE:
            x, y = arrays[first], arrays[second]
            if not len(x) or not len(y):
                continue
            test = mannwhitneyu(x, y, alternative="two-sided")
            rows.append({
                "fibroblast_state": state, "feature": feature,
                "analysis": "Pairwise Mann-Whitney",
                "stratum": f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}",
                "n_patients": len(x) + len(y),
                "median_effect": np.median(x) - np.median(y), "p_value": test.pvalue,
            })
    out = pd.DataFrame(rows)
    out["fdr_within_state_analysis_stratum"] = out.groupby(
        ["fibroblast_state", "analysis", "stratum"], dropna=False
    ).p_value.transform(bh_adjust)
    return out


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


cells = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv", low_memory=False)
crc = pd.read_csv(ROOT / "neutrophil_fibroblast_crc_marker_map" / "crc_associated_marker_expression.csv")
d = cells.merge(crc, on="cell_ID", how="left", validate="one_to_one")
thresholds = pd.read_csv(SOURCE / "tables" / "00_fibroblast_state_thresholds.csv").set_index("marker")
fap_threshold = float(thresholds.loc["FAP", "primary_threshold_clr"])
a5_threshold = float(thresholds.loc["a5B1", "primary_threshold_clr"])

d["is_fibroblast"] = d.cell_type.eq("Fibroblast")
d["is_neutrophil"] = d.cell_type.str.startswith("Neutrophil", na=False)
fap_pos = d.FAPa_cell_clr.ge(fap_threshold)
a5_pos = d.a5B1_cell_clr.ge(a5_threshold)
d["fibroblast_state"] = "Non-fibroblast"
d.loc[d.is_fibroblast & fap_pos & a5_pos, "fibroblast_state"] = "FAP+/a5B1+"
d.loc[d.is_fibroblast & fap_pos & ~a5_pos, "fibroblast_state"] = "FAP+/a5B1-"
d.loc[d.is_fibroblast & ~fap_pos & a5_pos, "fibroblast_state"] = "FAP-/a5B1+"
d.loc[d.is_fibroblast & ~fap_pos & ~a5_pos, "fibroblast_state"] = "FAP-/a5B1-"

marker_source = {
    "OSM": "raw_OSM", "CXCR4": "raw_CXCR4", "PDL1": "raw_PDL1",
    "HLA_ABC": "raw_HLA_ABC", "a5B1": "raw_a5B1", "PADI4": "PADI4",
    "MX1": "MX.1", "CD66b": "CD66b", "CD16": "CD16", "CD11b": "CD11b",
    "Ki67": "raw_Ki67",
}

sample_counts = d.groupby(["PatientID", "Diagnosis2"]).agg(
    total_cells=("cell_ID", "size"),
    n_fibroblasts=("is_fibroblast", "sum"),
    n_neutrophils=("is_neutrophil", "sum"),
).reset_index()
sample_counts["eligible_equal_count_spatial"] = (
    sample_counts.n_fibroblasts.ge(N_FIBROBLASTS)
    & sample_counts.n_neutrophils.ge(N_NEUTROPHILS)
)
sample_counts.to_csv(TAB / "00_sample_cell_counts_and_eligibility.csv", index=False)

# Direct count-exposure normalization of the original observed state-edge counts.
original_spatial = pd.read_csv(SOURCE / "tables" / "03_spatial_enrichment_by_patient.csv")
direct = original_spatial[original_spatial.radius_um.eq(RADIUS)].copy()
direct["edges_per_million_neutrophil_state_fibroblast_pairs"] = (
    direct.observed_state_edges / (direct.n_neutrophils * direct.n_state_fibroblasts) * 1e6
)
direct["all_edges_per_million_neutrophil_fibroblast_pairs"] = (
    direct.total_edges / (direct.n_neutrophils * direct.n_fibroblasts) * 1e6
)
direct["log2_state_to_all_pair_rate_ratio"] = np.log2(
    (direct.observed_state_edges + .5) / direct.n_state_fibroblasts
    / ((direct.total_edges + .5) / direct.n_fibroblasts)
)
direct.to_csv(TAB / "01_direct_pair_exposure_normalized_spatial_metrics.csv", index=False)

direct_stats_rows = []
for state, q in direct.groupby("fibroblast_state"):
    for metric in (
        "edges_per_million_neutrophil_state_fibroblast_pairs",
        "log2_state_to_all_pair_rate_ratio",
    ):
        direct_stats_rows.extend(group_tests(q, metric, f"{state}; {metric}", zero_tests=metric.startswith("log2")))
direct_stats = pd.DataFrame(direct_stats_rows)
direct_stats["fdr_within_endpoint_analysis"] = direct_stats.groupby(
    ["endpoint", "analysis"]
).p_value.transform(bh_adjust)
direct_stats.to_csv(TAB / "02_direct_pair_exposure_normalized_statistics.csv", index=False)

# Equal-count bootstrap: 250 fibroblasts and 30 neutrophils per sample/iteration.
rng = np.random.default_rng(SEED)
spatial_boot_rows = []
for pid, z in d.groupby("PatientID", sort=True):
    fib = z[z.is_fibroblast].reset_index(drop=True)
    neut = z[z.is_neutrophil].reset_index(drop=True)
    if len(fib) < N_FIBROBLASTS or len(neut) < N_NEUTROPHILS:
        continue
    for iteration in range(1, N_BOOT + 1):
        fi = rng.choice(len(fib), N_FIBROBLASTS, replace=False)
        ni = rng.choice(len(neut), N_NEUTROPHILS, replace=False)
        fs = fib.iloc[fi].reset_index(drop=True)
        ns = neut.iloc[ni]
        fib_xy = fs[["x", "y"]].to_numpy()
        degrees = np.zeros(N_FIBROBLASTS, dtype=np.int64)
        for ids in cKDTree(fib_xy).query_ball_point(ns[["x", "y"]].to_numpy(), RADIUS):
            if ids:
                np.add.at(degrees, np.asarray(ids, dtype=int), 1)
        total_edges = int(degrees.sum())
        for state in STATES:
            mask = fs.fibroblast_state.eq(state).to_numpy()
            n_state = int(mask.sum())
            if n_state == 0:
                continue
            state_edges = int(degrees[mask].sum())
            expected_state_edges = total_edges * n_state / N_FIBROBLASTS
            pair_rate = state_edges / (N_NEUTROPHILS * n_state) * 1e6
            all_pair_rate = total_edges / (N_NEUTROPHILS * N_FIBROBLASTS) * 1e6
            rate_ratio = ((state_edges + .5) / n_state) / ((total_edges + .5) / N_FIBROBLASTS)
            spatial_boot_rows.append({
                "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
                "iteration": iteration, "fibroblast_state": state,
                "sampled_fibroblasts": N_FIBROBLASTS,
                "sampled_neutrophils": N_NEUTROPHILS,
                "sampled_state_fibroblasts": n_state,
                "state_edges": state_edges, "expected_state_edges": expected_state_edges,
                "total_edges": total_edges,
                "edges_per_million_possible_pairs": pair_rate,
                "all_fibroblast_edges_per_million_possible_pairs": all_pair_rate,
                "log2_state_to_all_pair_rate_ratio": np.log2(rate_ratio),
            })

spatial_boot = pd.DataFrame(spatial_boot_rows)
spatial_boot.to_csv(TAB / "03_equal_count_spatial_bootstrap_iterations.csv", index=False)
spatial_patient = spatial_boot.groupby(["PatientID", "Diagnosis2", "fibroblast_state"]).agg(
    n_valid_iterations=("iteration", "nunique"),
    median_sampled_state_fibroblasts=("sampled_state_fibroblasts", "median"),
    median_edges_per_million_pairs=("edges_per_million_possible_pairs", "median"),
    summed_observed_state_edges=("state_edges", "sum"),
    summed_expected_state_edges=("expected_state_edges", "sum"),
    median_log2_state_to_all_pair_rate_ratio=("log2_state_to_all_pair_rate_ratio", "median"),
    ci025_log2_rate_ratio=("log2_state_to_all_pair_rate_ratio", lambda x: x.quantile(.025)),
    ci975_log2_rate_ratio=("log2_state_to_all_pair_rate_ratio", lambda x: x.quantile(.975)),
).reset_index()
spatial_patient["aggregate_log2_observed_expected"] = np.log2(
    (spatial_patient.summed_observed_state_edges + .5)
    / (spatial_patient.summed_expected_state_edges + .5)
)
spatial_patient["passes_minimum_valid_bootstraps"] = spatial_patient.n_valid_iterations.ge(MIN_VALID_BOOTSTRAPS)
spatial_patient.to_csv(TAB / "04_equal_count_spatial_summary_by_patient.csv", index=False)

spatial_valid = spatial_patient[spatial_patient.passes_minimum_valid_bootstraps].copy()
spatial_stats_rows = []
for state, q in spatial_valid.groupby("fibroblast_state"):
    spatial_stats_rows.extend(group_tests(
        q, "aggregate_log2_observed_expected",
        f"{state}; equal 250 fibroblast/30 neutrophil aggregate log2 O/E",
        zero_tests=True,
    ))
spatial_stats = pd.DataFrame(spatial_stats_rows)
spatial_stats["fdr_within_endpoint_analysis"] = spatial_stats.groupby(
    ["endpoint", "analysis"]
).p_value.transform(bh_adjust)
spatial_stats.to_csv(TAB / "05_equal_count_spatial_statistics.csv", index=False)

# Equal near/far counts for marker effects. Each patient contributes a bootstrap
# median based on exactly five proximal and ten distant neutrophils per iteration.
distances = pd.read_csv(SOURCE / "tables" / "04_neutrophil_fibroblast_state_distances.csv")
distance_wide = distances.pivot(index="cell_ID", columns="fibroblast_state", values="nearest_state_distance_um")
neut = d[d.is_neutrophil].copy().merge(
    distance_wide, left_on="cell_ID", right_index=True, how="left", validate="one_to_one"
)
for marker, source in marker_source.items():
    neut[f"log2_{marker}"] = np.log2(neut[source].astype(float).clip(lower=0) + .1)
for marker in PROGRAM_MARKERS:
    col = f"log2_{marker}"
    neut[f"z_{marker}"] = neut.groupby("PatientID")[col].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0.0
    )
neut[PROGRAM_NAME] = neut[[f"z_{marker}" for marker in PROGRAM_MARKERS]].mean(axis=1)

marker_boot_rows = []
rng_marker = np.random.default_rng(SEED + 1)
features = list(MARKERS) + [PROGRAM_NAME]
for (pid, diagnosis), q in neut.groupby(["PatientID", "Diagnosis2"]):
    for state in STATES:
        prox = q[q[state].le(RADIUS)]
        distant = q[q[state].gt(100)]
        if len(prox) < N_PROXIMAL or len(distant) < N_DISTANT:
            continue
        prox_indices = rng_marker.integers(0, len(prox), size=(N_BOOT, N_PROXIMAL))
        dist_indices = rng_marker.integers(0, len(distant), size=(N_BOOT, N_DISTANT))
        for feature in features:
            p = (prox[PROGRAM_NAME] if feature == PROGRAM_NAME else prox[f"log2_{feature}"]).to_numpy()
            f = (distant[PROGRAM_NAME] if feature == PROGRAM_NAME else distant[f"log2_{feature}"]).to_numpy()
            effects = np.median(p[prox_indices], axis=1) - np.median(f[dist_indices], axis=1)
            marker_boot_rows.append({
                "PatientID": pid, "Diagnosis2": diagnosis, "fibroblast_state": state,
                "feature": feature, "n_available_proximal": len(prox),
                "n_available_distant": len(distant), "sampled_proximal_per_iteration": N_PROXIMAL,
                "sampled_distant_per_iteration": N_DISTANT, "n_bootstraps": N_BOOT,
                "bootstrap_median_effect": np.median(effects),
                "bootstrap_ci025": np.quantile(effects, .025),
                "bootstrap_ci975": np.quantile(effects, .975),
            })

marker_boot = pd.DataFrame(marker_boot_rows)
marker_boot.to_csv(TAB / "06_equal_count_marker_effects_by_patient.csv", index=False)
marker_stats = feature_statistics(marker_boot)
marker_stats.to_csv(TAB / "07_equal_count_marker_statistics.csv", index=False)

# Plot normalized results.
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})
group_order = [GROUP_LABEL[g] for g in GROUPS]
palette = {GROUP_LABEL[g]: GROUP_COLOR[g] for g in GROUPS}

fig, axes = plt.subplots(2, 2, figsize=(16, 13))

q = direct[direct.fibroblast_state.eq("FAP+/a5B1+")].copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="Group", y="log2_state_to_all_pair_rate_ratio", hue="Group",
            order=group_order, palette=palette, legend=False, showfliers=False, ax=axes[0, 0])
sns.stripplot(data=q, x="Group", y="log2_state_to_all_pair_rate_ratio", order=group_order,
              color="black", size=5, ax=axes[0, 0])
axes[0, 0].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0, 0].set_title("A  Direct pair-exposure normalization")
axes[0, 0].set_xlabel("")
axes[0, 0].set_ylabel("FAP+/a5B1+ log2 state/all edge rate")
axes[0, 0].tick_params(axis="x", rotation=15)

q = spatial_valid[spatial_valid.fibroblast_state.eq("FAP+/a5B1+")].copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="Group", y="aggregate_log2_observed_expected", hue="Group",
            order=group_order, palette=palette, legend=False, showfliers=False, ax=axes[0, 1])
sns.stripplot(data=q, x="Group", y="aggregate_log2_observed_expected", order=group_order,
              color="black", size=5, ax=axes[0, 1])
axes[0, 1].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0, 1].set_title("B  Equal 250-fibroblast/30-neutrophil bootstrap")
axes[0, 1].set_xlabel("")
axes[0, 1].set_ylabel("Aggregate bootstrap log2 observed / expected")
axes[0, 1].tick_params(axis="x", rotation=15)

sp_heat = spatial_valid.groupby(["Diagnosis2", "fibroblast_state"]).aggregate_log2_observed_expected.median().unstack()
sp_heat = sp_heat.reindex(index=GROUPS, columns=STATES)
sns.heatmap(sp_heat, annot=True, fmt=".2f", cmap="vlag", center=0,
            cbar_kws={"label": "Median normalized log2 rate"}, ax=axes[1, 0])
axes[1, 0].set_yticklabels(group_order, rotation=0)
axes[1, 0].set_title("C  Equal-count spatial effects by state")
axes[1, 0].set_xlabel("")
axes[1, 0].set_ylabel("")

heat = marker_boot[marker_boot.fibroblast_state.eq("FAP+/a5B1-")].groupby(
    ["Diagnosis2", "feature"]
).bootstrap_median_effect.median().unstack().reindex(index=GROUPS, columns=features)
sns.heatmap(heat, annot=False, cmap="vlag", center=0,
            cbar_kws={"label": "Equal-count median near-far effect"}, ax=axes[1, 1])
axes[1, 1].set_yticklabels(group_order, rotation=0)
axes[1, 1].set_title("D  FAP+/a5B1- neutrophil phenotype")
axes[1, 1].set_xlabel("")
axes[1, 1].set_ylabel("")
axes[1, 1].tick_params(axis="x", rotation=45, labelsize=9)

fig.suptitle("Cell-count-normalized FAP/a5B1 fibroblast-neutrophil analysis", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_Cell_Count_Normalized_FAP_A5B1_Neutrophil")


def fmt(x):
    return f"{x:.3g}"


dp_endpoint = "FAP+/a5B1+; equal 250 fibroblast/30 neutrophil aggregate log2 O/E"
dp_stats = spatial_stats[spatial_stats.endpoint.eq(dp_endpoint)]
dp_omnibus = dp_stats[dp_stats.analysis.eq("Kruskal-Wallis omnibus")].iloc[0]
dp_pairwise = dp_stats[dp_stats.analysis.eq("Pairwise Mann-Whitney")].set_index("comparison")
dp_within = dp_stats[dp_stats.analysis.eq("One-sample Wilcoxon versus zero")].set_index("group_1")
dp_medians = spatial_valid[spatial_valid.fibroblast_state.eq("FAP+/a5B1+")].groupby(
    "Diagnosis2"
).aggregate_log2_observed_expected.median()

pair_lines = []
for first, second in PAIRWISE:
    label = f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}"
    row = dp_pairwise.loc[label]
    pair_lines.append(
        f"- {label}: difference {row.effect:.3f}, P={fmt(row.p_value)}, "
        f"FDR={fmt(row.fdr_within_endpoint_analysis)}."
    )
within_lines = []
for group in GROUPS:
    if group not in dp_within.index:
        continue
    row = dp_within.loc[group]
    within_lines.append(
        f"- {GROUP_LABEL[group]}: median {row.effect:.3f}, P={fmt(row.p_value)}, "
        f"FDR={fmt(row.fdr_within_endpoint_analysis)}."
    )

fap_only_sig = marker_stats[
    marker_stats.fibroblast_state.eq("FAP+/a5B1-")
    & marker_stats.analysis.eq("One-sample Wilcoxon versus zero")
    & marker_stats.fdr_within_state_analysis_stratum.lt(.05)
].sort_values(["stratum", "median_effect"], ascending=[True, False])
marker_lines = []
for group in GROUPS:
    label = GROUP_LABEL[group]
    q = fap_only_sig[fap_only_sig.stratum.eq(label)]
    hits = [f"{r.feature} ({r.median_effect:+.2f})" for r in q.itertuples()]
    marker_lines.append(f"- {label}: {', '.join(hits) if hits else 'none after FDR correction'}.")

state_medians = spatial_valid.groupby(
    ["fibroblast_state", "Diagnosis2"]
).aggregate_log2_observed_expected.median()
dn_endpoint = "FAP-/a5B1-; equal 250 fibroblast/30 neutrophil aggregate log2 O/E"
dn_contrast = spatial_stats[
    spatial_stats.endpoint.eq(dn_endpoint)
    & spatial_stats.analysis.eq("Pairwise Mann-Whitney")
    & spatial_stats.comparison.eq("UC inflamed vs UC noninflamed")
].iloc[0]
a5_only_endpoint = "FAP-/a5B1+; equal 250 fibroblast/30 neutrophil aggregate log2 O/E"
a5_only_within = spatial_stats[
    spatial_stats.endpoint.eq(a5_only_endpoint)
    & spatial_stats.analysis.eq("One-sample Wilcoxon versus zero")
].set_index("group_1")

excluded = sample_counts[~sample_counts.eligible_equal_count_spatial]
excluded_text = ", ".join(
    f"{row.PatientID} ({int(row.n_neutrophils)} neutrophils, {int(row.n_fibroblasts)} fibroblasts)"
    for row in excluded.itertuples()
) or "none"

report = f"""# Cell-count-normalized FAP/a5B1 fibroblast-neutrophil analysis

## Normalization strategy

The original permutation enrichment already gave each patient equal inferential weight and controlled fibroblast-state abundance. This sensitivity analysis adds two explicit cell-count controls: (1) observed edges divided by the number of possible neutrophil-state-fibroblast pairs, and (2) {N_BOOT} repeated subsamples containing exactly {N_FIBROBLASTS} fibroblasts and {N_NEUTROPHILS} neutrophils per sample. Marker comparisons used exactly {N_PROXIMAL} proximal and {N_DISTANT} distant neutrophils per bootstrap iteration. Fibroblast-state abundance remained a fraction of total fibroblasts per sample.

The equal-count spatial analysis excluded: {excluded_text}. Samples were retained only when at least {MIN_VALID_BOOTSTRAPS} iterations contained the fibroblast state being evaluated.

## Equal-count FAP+/a5B1+ spatial result

Median aggregate bootstrap log2 observed/expected edge ratios were {dp_medians.get('Control', np.nan):.3f} in controls, {dp_medians.get('UC_Noninflamed', np.nan):.3f} in UC noninflamed, and {dp_medians.get('UC_Inflamed', np.nan):.3f} in UC inflamed tissue. The three-group omnibus P value was {fmt(dp_omnibus.p_value)}. Pairwise results were:

{chr(10).join(pair_lines)}

Within-group tests against no preferential association were:

{chr(10).join(within_lines)}

## Other fibroblast states after equal-count normalization

The FAP-/a5B1+ state also remained depleted in UC noninflamed tissue (median {state_medians.loc[('FAP-/a5B1+', 'UC_Noninflamed')]:.3f}, FDR={fmt(a5_only_within.loc['UC_Noninflamed', 'fdr_within_endpoint_analysis'])}) and UC inflamed tissue (median {state_medians.loc[('FAP-/a5B1+', 'UC_Inflamed')]:.3f}, FDR={fmt(a5_only_within.loc['UC_Inflamed', 'fdr_within_endpoint_analysis'])}). Conversely, FAP-/a5B1- fibroblasts were enriched near neutrophils in UC noninflamed tissue (median {state_medians.loc[('FAP-/a5B1-', 'UC_Noninflamed')]:.3f}) but not inflamed tissue ({state_medians.loc[('FAP-/a5B1-', 'UC_Inflamed')]:.3f}); the inflamed-minus-noninflamed difference was {dn_contrast.effect:.3f} (P={fmt(dn_contrast.p_value)}, FDR={fmt(dn_contrast.fdr_within_endpoint_analysis)}).

## Equal-count neutrophil phenotype near FAP+/a5B1- fibroblasts

After fixing near/far neutrophil counts in every iteration, significant marker effects were:

{chr(10).join(marker_lines)}

## Interpretation

These analyses remove differences in the numbers of fibroblasts and neutrophils contributing to each sample-level estimate. Results should still be interpreted as spatial and phenotypic associations rather than evidence of fibroblast-driven neutrophil reprogramming or causal signaling.
"""
(OUT / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

manifest = {
    "source_analysis": str(SOURCE), "radius_um": RADIUS,
    "equal_fibroblasts_per_iteration": N_FIBROBLASTS,
    "equal_neutrophils_per_iteration": N_NEUTROPHILS,
    "equal_proximal_neutrophils_per_marker_iteration": N_PROXIMAL,
    "equal_distant_neutrophils_per_marker_iteration": N_DISTANT,
    "bootstrap_iterations": N_BOOT, "minimum_valid_iterations": MIN_VALID_BOOTSTRAPS,
    "seed": SEED, "biological_replicate": "patient",
    "excluded_equal_count_samples": excluded.PatientID.tolist(),
}
(OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print("Completed cell-count-normalized FAP/a5B1 neutrophil analysis.")
print(f"Saved outputs to {OUT}")
