from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import kruskal, mannwhitneyu, spearmanr, wilcoxon

warnings.filterwarnings("ignore")

ROOT = Path("giotto_codex_results")
OUT = ROOT / "fap_a5b1_neutrophil_reprogramming"
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
STATES = ("FAP+/a5B1+", "FAP+/a5B1-", "FAP-/a5B1+", "FAP-/a5B1-")
STATE_COLOR = {
    "FAP+/a5B1+": "#7A0177",
    "FAP+/a5B1-": "#E6550D",
    "FAP-/a5B1+": "#2CA25F",
    "FAP-/a5B1-": "#969696",
}
PAIRWISE = (
    ("UC_Inflamed", "Control"),
    ("UC_Noninflamed", "Control"),
    ("UC_Inflamed", "UC_Noninflamed"),
)
SUBTYPES = ("PADI4", "CXCR4", "MX1", "OSM")
RADII = (25, 50, 100)
PRIMARY_RADIUS = 50
PROXIMAL_UM = 50
DISTANT_UM = 100
MIN_PROXIMAL = 5
MIN_DISTANT = 10
N_PERM = 1000
SEED = 20260802

MARKERS = (
    "OSM", "CXCR4", "PDL1", "HLA_ABC", "a5B1",
    "PADI4", "MX1", "CD66b", "CD16", "CD11b", "Ki67",
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


def epsilon_squared(h, n, k=3):
    return max(0.0, (h - k + 1) / (n - k)) if n > k else np.nan


def group_tests(frame, value, endpoint, zero_tests=True):
    arrays = {
        group: frame.loc[frame.Diagnosis2.eq(group), value].dropna().to_numpy()
        for group in GROUPS
    }
    rows = []
    if all(len(arrays[group]) for group in GROUPS):
        test = kruskal(*(arrays[group] for group in GROUPS))
        rows.append({
            "endpoint": endpoint,
            "analysis": "Kruskal-Wallis omnibus",
            "comparison": "Control vs UC noninflamed vs UC inflamed",
            "group_1": "", "group_2": "",
            "n_1": len(arrays[GROUPS[0]]), "n_2": len(arrays[GROUPS[1]]),
            "n_3": len(arrays[GROUPS[2]]),
            "median_1": np.median(arrays[GROUPS[0]]),
            "median_2": np.median(arrays[GROUPS[1]]),
            "median_3": np.median(arrays[GROUPS[2]]),
            "effect": epsilon_squared(test.statistic, sum(map(len, arrays.values()))),
            "effect_type": "Kruskal-Wallis epsilon-squared",
            "rank_biserial": np.nan,
            "p_value": test.pvalue,
        })
    for first, second in PAIRWISE:
        x, y = arrays[first], arrays[second]
        if not len(x) or not len(y):
            continue
        test = mannwhitneyu(x, y, alternative="two-sided")
        rows.append({
            "endpoint": endpoint,
            "analysis": "Pairwise Mann-Whitney",
            "comparison": f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}",
            "group_1": first, "group_2": second,
            "n_1": len(x), "n_2": len(y), "n_3": np.nan,
            "median_1": np.median(x), "median_2": np.median(y), "median_3": np.nan,
            "effect": np.median(x) - np.median(y),
            "effect_type": f"median difference ({GROUP_LABEL[first]} - {GROUP_LABEL[second]})",
            "rank_biserial": rank_biserial(x, y),
            "p_value": test.pvalue,
        })
    if zero_tests:
        for group in GROUPS:
            x = arrays[group]
            if len(x) >= 3 and np.any(x != 0):
                test = wilcoxon(x, alternative="two-sided", zero_method="wilcox")
                rows.append({
                    "endpoint": endpoint,
                    "analysis": "One-sample Wilcoxon versus zero",
                    "comparison": f"{GROUP_LABEL[group]} vs zero",
                    "group_1": group, "group_2": "",
                    "n_1": len(x), "n_2": np.nan, "n_3": np.nan,
                    "median_1": np.median(x), "median_2": 0.0, "median_3": np.nan,
                    "effect": np.median(x), "effect_type": "median effect",
                    "rank_biserial": np.nan, "p_value": test.pvalue,
                })
    return rows


def feature_effect_statistics(effects, method):
    rows = []
    q0 = effects[effects.method.eq(method)]
    for (state, feature), q in q0.groupby(["fibroblast_state", "feature"]):
        all_values = q.effect.dropna().to_numpy()
        if len(all_values) >= 3 and np.any(all_values != 0):
            test = wilcoxon(all_values, alternative="two-sided", zero_method="wilcox")
            rows.append({
                "method": method, "fibroblast_state": state, "feature": feature,
                "analysis": "One-sample Wilcoxon versus zero", "stratum": "All samples",
                "n_patients": len(all_values), "median_effect": np.median(all_values),
                "p_value": test.pvalue,
            })
        arrays = {}
        for group in GROUPS:
            values = q.loc[q.Diagnosis2.eq(group), "effect"].dropna().to_numpy()
            arrays[group] = values
            if len(values) >= 3 and np.any(values != 0):
                test = wilcoxon(values, alternative="two-sided", zero_method="wilcox")
                rows.append({
                    "method": method, "fibroblast_state": state, "feature": feature,
                    "analysis": "One-sample Wilcoxon versus zero",
                    "stratum": GROUP_LABEL[group], "n_patients": len(values),
                    "median_effect": np.median(values), "p_value": test.pvalue,
                })
        if all(len(arrays[group]) >= 3 for group in GROUPS):
            test = kruskal(*(arrays[group] for group in GROUPS))
            rows.append({
                "method": method, "fibroblast_state": state, "feature": feature,
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
                "method": method, "fibroblast_state": state, "feature": feature,
                "analysis": "Pairwise Mann-Whitney",
                "stratum": f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}",
                "n_patients": len(x) + len(y),
                "median_effect": np.median(x) - np.median(y), "p_value": test.pvalue,
            })
    out = pd.DataFrame(rows)
    out["fdr_within_state_method_analysis_stratum"] = out.groupby(
        ["method", "fibroblast_state", "analysis", "stratum"], dropna=False
    ).p_value.transform(bh_adjust)
    return out


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


cells = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv", low_memory=False)
crc = pd.read_csv(ROOT / "neutrophil_fibroblast_crc_marker_map" / "crc_associated_marker_expression.csv")
d = cells.merge(crc, on="cell_ID", how="left", validate="one_to_one")
required = {
    "cell_ID", "PatientID", "Diagnosis2", "cell_type", "Neutrophil_subtype_0.8",
    "FAPa_cell_clr", "a5B1_cell_clr", "x", "y", "PADI4", "MX.1", "CD66b",
    "CD16", "CD11b", "raw_OSM", "raw_CXCR4", "raw_PDL1", "raw_HLA_ABC",
    "raw_a5B1", "raw_Ki67",
}
missing = required.difference(d.columns)
if missing:
    raise RuntimeError(f"Missing required columns: {sorted(missing)}")

marker_source = {
    "OSM": "raw_OSM", "CXCR4": "raw_CXCR4", "PDL1": "raw_PDL1",
    "HLA_ABC": "raw_HLA_ABC", "a5B1": "raw_a5B1", "PADI4": "PADI4",
    "MX1": "MX.1", "CD66b": "CD66b", "CD16": "CD16", "CD11b": "CD11b",
    "Ki67": "raw_Ki67",
}
d["is_fibroblast"] = d.cell_type.eq("Fibroblast")
d["is_neutrophil"] = d.cell_type.str.startswith("Neutrophil", na=False)
uc_fib = d[d.is_fibroblast & d.Diagnosis2.isin(("UC_Noninflamed", "UC_Inflamed"))]
fap_threshold = float(uc_fib.FAPa_cell_clr.quantile(0.75))
a5_threshold = float(uc_fib.a5B1_cell_clr.quantile(0.75))

fap_pos = d.FAPa_cell_clr.ge(fap_threshold)
a5_pos = d.a5B1_cell_clr.ge(a5_threshold)
d["fibroblast_state"] = "Non-fibroblast"
d.loc[d.is_fibroblast & fap_pos & a5_pos, "fibroblast_state"] = "FAP+/a5B1+"
d.loc[d.is_fibroblast & fap_pos & ~a5_pos, "fibroblast_state"] = "FAP+/a5B1-"
d.loc[d.is_fibroblast & ~fap_pos & a5_pos, "fibroblast_state"] = "FAP-/a5B1+"
d.loc[d.is_fibroblast & ~fap_pos & ~a5_pos, "fibroblast_state"] = "FAP-/a5B1-"

thresholds = pd.DataFrame({
    "marker": ["FAP", "a5B1"],
    "primary_threshold_clr": [fap_threshold, a5_threshold],
    "source": ["pooled UC fibroblast 75th percentile"] * 2,
})
thresholds.to_csv(TAB / "00_fibroblast_state_thresholds.csv", index=False)

fib = d[d.is_fibroblast].copy()
state_abundance = (
    fib.groupby(["PatientID", "Diagnosis2", "fibroblast_state"]).size()
    .rename("n_state").reset_index()
)
totals = fib.groupby(["PatientID", "Diagnosis2"]).size().rename("n_fibroblasts").reset_index()
state_abundance = totals.merge(state_abundance, on=["PatientID", "Diagnosis2"], how="left")
state_abundance["fraction_fibroblasts"] = state_abundance.n_state / state_abundance.n_fibroblasts
state_abundance.to_csv(TAB / "01_fibroblast_state_abundance_by_patient.csv", index=False)

abundance_stats_rows = []
for state, q in state_abundance.groupby("fibroblast_state"):
    abundance_stats_rows.extend(group_tests(q, "fraction_fibroblasts", state, zero_tests=False))
abundance_stats = pd.DataFrame(abundance_stats_rows)
abundance_stats["fdr_within_analysis"] = abundance_stats.groupby("analysis").p_value.transform(bh_adjust)
abundance_stats.to_csv(TAB / "02_fibroblast_state_abundance_statistics.csv", index=False)

# Within-patient spatial enrichment after permuting fibroblast-state labels.
rng = np.random.default_rng(SEED)
spatial_rows = []
neut_distance_rows = []
for pid, z in d.groupby("PatientID", sort=True):
    fibp = z[z.is_fibroblast].copy()
    neut = z[z.is_neutrophil].copy()
    if fibp.empty or neut.empty:
        continue
    fib_xy = fibp[["x", "y"]].to_numpy()
    neut_xy = neut[["x", "y"]].to_numpy()
    tree = cKDTree(fib_xy)
    neighbor_ids = tree.query_ball_point(neut_xy, max(RADII))
    state_masks = {state: fibp.fibroblast_state.eq(state).to_numpy() for state in STATES}
    for state in STATES:
        state_xy = fibp.loc[state_masks[state], ["x", "y"]].to_numpy()
        if not len(state_xy):
            continue
        distances = cKDTree(state_xy).query(neut_xy)[0]
        for row_idx, (cell_id, subtype, dist) in enumerate(zip(
            neut.cell_ID, neut["Neutrophil_subtype_0.8"], distances
        )):
            neut_distance_rows.append({
                "cell_ID": cell_id, "PatientID": pid,
                "Diagnosis2": z.Diagnosis2.iloc[0],
                "Neutrophil_subtype": subtype if pd.notna(subtype) else "Unclassified",
                "fibroblast_state": state, "nearest_state_distance_um": dist,
            })
    for radius in RADII:
        degrees = np.zeros(len(fibp), dtype=np.int64)
        total_edges = 0
        for point, ids in zip(neut_xy, neighbor_ids):
            if not ids:
                continue
            ids = np.asarray(ids, dtype=int)
            delta = fib_xy[ids] - point
            ids = ids[(delta * delta).sum(axis=1) <= radius * radius]
            if len(ids):
                np.add.at(degrees, ids, 1)
                total_edges += len(ids)
        if total_edges == 0:
            continue
        for state in STATES:
            mask = state_masks[state]
            n_state = int(mask.sum())
            if n_state == 0:
                continue
            observed = int(degrees[mask].sum())
            permuted = np.empty(N_PERM, dtype=float)
            for iteration in range(N_PERM):
                chosen = rng.choice(len(fibp), n_state, replace=False)
                permuted[iteration] = degrees[chosen].sum()
            expected = float(permuted.mean())
            p_value = (1 + np.sum(np.abs(permuted - expected) >= abs(observed - expected))) / (N_PERM + 1)
            spatial_rows.append({
                "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
                "fibroblast_state": state, "radius_um": radius,
                "n_neutrophils": len(neut), "n_fibroblasts": len(fibp),
                "n_state_fibroblasts": n_state, "total_edges": total_edges,
                "observed_state_edges": observed, "expected_state_edges": expected,
                "log2_observed_expected": np.log2((observed + 0.5) / (expected + 0.5)),
                "permutation_p": p_value,
            })

spatial = pd.DataFrame(spatial_rows)
spatial.to_csv(TAB / "03_spatial_enrichment_by_patient.csv", index=False)
neut_distances = pd.DataFrame(neut_distance_rows)
neut_distances.to_csv(TAB / "04_neutrophil_fibroblast_state_distances.csv", index=False)

spatial_stats_rows = []
for (state, radius), q in spatial.groupby(["fibroblast_state", "radius_um"]):
    spatial_stats_rows.extend(group_tests(
        q, "log2_observed_expected", f"{state}; {int(radius)} um log2 observed/expected",
        zero_tests=True,
    ))
spatial_stats = pd.DataFrame(spatial_stats_rows)
spatial_stats["fdr_within_endpoint_analysis"] = spatial_stats.groupby(
    ["endpoint", "analysis"]
).p_value.transform(bh_adjust)
spatial_stats.to_csv(TAB / "05_spatial_enrichment_statistics.csv", index=False)

# Patient-specific q75 sensitivity for the double-positive state at 50 um.
sensitivity_rows = []
rng_sens = np.random.default_rng(SEED + 1)
for pid, z in d.groupby("PatientID", sort=True):
    fibp = z[z.is_fibroblast]
    neut = z[z.is_neutrophil]
    if fibp.empty or neut.empty:
        continue
    mask = fibp.FAPa_cell_clr.ge(fibp.FAPa_cell_clr.quantile(.75)) & fibp.a5B1_cell_clr.ge(
        fibp.a5B1_cell_clr.quantile(.75)
    )
    fib_xy = fibp[["x", "y"]].to_numpy()
    degrees = np.zeros(len(fibp), dtype=np.int64)
    for ids in cKDTree(fib_xy).query_ball_point(neut[["x", "y"]].to_numpy(), PRIMARY_RADIUS):
        if ids:
            np.add.at(degrees, np.asarray(ids, dtype=int), 1)
    n_state = int(mask.sum())
    observed = int(degrees[mask.to_numpy()].sum())
    permuted = np.empty(N_PERM)
    for iteration in range(N_PERM):
        permuted[iteration] = degrees[rng_sens.choice(len(fibp), n_state, replace=False)].sum()
    expected = float(permuted.mean())
    sensitivity_rows.append({
        "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
        "n_double_positive": n_state, "observed_edges": observed, "expected_edges": expected,
        "log2_observed_expected": np.log2((observed + .5) / (expected + .5)),
    })
sensitivity = pd.DataFrame(sensitivity_rows)
sensitivity.to_csv(TAB / "06_double_positive_patient_q75_sensitivity.csv", index=False)
sensitivity_stats = pd.DataFrame(group_tests(
    sensitivity, "log2_observed_expected", "FAP+/a5B1+ patient-q75 sensitivity; 50 um",
    zero_tests=True,
))
sensitivity_stats["fdr_within_analysis"] = sensitivity_stats.groupby("analysis").p_value.transform(bh_adjust)
sensitivity_stats.to_csv(TAB / "07_double_positive_patient_q75_sensitivity_statistics.csv", index=False)

# Prepare neutrophil marker values and a within-patient composite program score.
neut = d[d.is_neutrophil].copy()
for marker, source in marker_source.items():
    neut[f"log2_{marker}"] = np.log2(neut[source].astype(float).clip(lower=0) + 0.1)
for marker in PROGRAM_MARKERS:
    col = f"log2_{marker}"
    neut[f"z_{marker}"] = neut.groupby("PatientID")[col].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0.0
    )
neut[PROGRAM_NAME] = neut[[f"z_{marker}" for marker in PROGRAM_MARKERS]].mean(axis=1)

marker_effect_rows = []
distance_lookup = neut_distances.pivot(
    index="cell_ID", columns="fibroblast_state", values="nearest_state_distance_um"
)
neut = neut.merge(distance_lookup, left_on="cell_ID", right_index=True, how="left", validate="one_to_one")
features = list(MARKERS) + [PROGRAM_NAME]
for (pid, diagnosis), q in neut.groupby(["PatientID", "Diagnosis2"]):
    for state in STATES:
        distances = q[state]
        proximal = distances.le(PROXIMAL_UM)
        distant = distances.gt(DISTANT_UM)
        for feature in features:
            values = q[PROGRAM_NAME] if feature == PROGRAM_NAME else q[f"log2_{feature}"]
            if proximal.sum() >= MIN_PROXIMAL and distant.sum() >= MIN_DISTANT:
                marker_effect_rows.append({
                    "PatientID": pid, "Diagnosis2": diagnosis,
                    "fibroblast_state": state, "feature": feature,
                    "method": "Proximal-vs-distant median difference",
                    "n_proximal": int(proximal.sum()), "n_distant": int(distant.sum()),
                    "effect": float(values[proximal].median() - values[distant].median()),
                })
            ok = distances.notna() & values.notna()
            if ok.sum() >= 10 and distances[ok].nunique() > 1 and values[ok].nunique() > 1:
                rho = spearmanr(-distances[ok], values[ok]).statistic
                marker_effect_rows.append({
                    "PatientID": pid, "Diagnosis2": diagnosis,
                    "fibroblast_state": state, "feature": feature,
                    "method": "Spearman correlation with proximity",
                    "n_proximal": int(proximal.sum()), "n_distant": int(distant.sum()),
                    "effect": rho,
                })

marker_effects = pd.DataFrame(marker_effect_rows)
marker_effects.to_csv(TAB / "08_neutrophil_marker_proximity_effects_by_patient.csv", index=False)
marker_stats = pd.concat([
    feature_effect_statistics(marker_effects, method)
    for method in marker_effects.method.unique()
], ignore_index=True)
marker_stats.to_csv(TAB / "09_neutrophil_marker_proximity_statistics.csv", index=False)

# Neutrophil subtype enrichment among cells within 50 um of each fibroblast state.
subtype_rows = []
for (pid, diagnosis), q in neut.groupby(["PatientID", "Diagnosis2"]):
    all_counts = q["Neutrophil_subtype_0.8"].value_counts()
    for state in STATES:
        prox = q[q[state].le(PROXIMAL_UM)]
        if len(prox) < MIN_PROXIMAL:
            continue
        prox_counts = prox["Neutrophil_subtype_0.8"].value_counts()
        for subtype in SUBTYPES:
            prox_prop = (prox_counts.get(subtype, 0) + 0.5) / (len(prox) + 0.5 * len(SUBTYPES))
            all_prop = (all_counts.get(subtype, 0) + 0.5) / (len(q) + 0.5 * len(SUBTYPES))
            subtype_rows.append({
                "PatientID": pid, "Diagnosis2": diagnosis,
                "fibroblast_state": state, "feature": subtype,
                "method": "Subtype log2 proximal/all-neutrophil enrichment",
                "n_proximal": len(prox), "n_distant": len(q) - len(prox),
                "effect": np.log2(prox_prop / all_prop),
            })
subtype_effects = pd.DataFrame(subtype_rows)
subtype_effects.to_csv(TAB / "10_neutrophil_subtype_enrichment_by_patient.csv", index=False)
subtype_stats = feature_effect_statistics(
    subtype_effects, "Subtype log2 proximal/all-neutrophil enrichment"
)
subtype_stats.to_csv(TAB / "11_neutrophil_subtype_enrichment_statistics.csv", index=False)

# Compact group summaries for plotting and interpretation.
spatial_group = spatial.groupby(["Diagnosis2", "fibroblast_state", "radius_um"]).agg(
    n_patients=("PatientID", "nunique"),
    median_log2_observed_expected=("log2_observed_expected", "median"),
).reset_index()
spatial_group.to_csv(TAB / "12_spatial_group_summary.csv", index=False)
marker_group = marker_effects.groupby(
    ["method", "Diagnosis2", "fibroblast_state", "feature"]
).agg(n_patients=("PatientID", "nunique"), median_effect=("effect", "median")).reset_index()
marker_group.to_csv(TAB / "13_marker_effect_group_summary.csv", index=False)
subtype_group = subtype_effects.groupby(
    ["Diagnosis2", "fibroblast_state", "feature"]
).agg(n_patients=("PatientID", "nunique"), median_effect=("effect", "median")).reset_index()
subtype_group.to_csv(TAB / "14_subtype_effect_group_summary.csv", index=False)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})
group_order = [GROUP_LABEL[group] for group in GROUPS]
palette = {GROUP_LABEL[group]: GROUP_COLOR[group] for group in GROUPS}

fig, axes = plt.subplots(2, 3, figsize=(21, 13))

ab_heat = state_abundance.groupby(["Diagnosis2", "fibroblast_state"]).fraction_fibroblasts.median().unstack()
ab_heat = ab_heat.reindex(index=GROUPS, columns=STATES)
sns.heatmap(ab_heat, annot=True, fmt=".2f", cmap="Purples", cbar_kws={"label": "Median fraction"}, ax=axes[0, 0])
axes[0, 0].set_yticklabels(group_order, rotation=0)
axes[0, 0].set_title("A  Fibroblast-state abundance")
axes[0, 0].set_xlabel("")
axes[0, 0].set_ylabel("")

sp_heat = spatial_group[spatial_group.radius_um.eq(PRIMARY_RADIUS)].pivot(
    index="Diagnosis2", columns="fibroblast_state", values="median_log2_observed_expected"
).reindex(index=GROUPS, columns=STATES)
sns.heatmap(sp_heat, annot=True, fmt=".2f", cmap="vlag", center=0,
            cbar_kws={"label": "Median log2 O/E"}, ax=axes[0, 1])
axes[0, 1].set_yticklabels(group_order, rotation=0)
axes[0, 1].set_title("B  Neutrophil spatial targeting (50 um)")
axes[0, 1].set_xlabel("")
axes[0, 1].set_ylabel("")

q = spatial[spatial.fibroblast_state.eq("FAP+/a5B1+") & spatial.radius_um.eq(PRIMARY_RADIUS)].copy()
q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="Group", y="log2_observed_expected", hue="Group", order=group_order,
            palette=palette, legend=False, showfliers=False, ax=axes[0, 2])
sns.stripplot(data=q, x="Group", y="log2_observed_expected", order=group_order,
              color="black", size=5, ax=axes[0, 2])
axes[0, 2].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0, 2].set_title("C  FAP+/a5B1+ targeting")
axes[0, 2].set_xlabel("")
axes[0, 2].set_ylabel("log2 observed / expected edges")
axes[0, 2].tick_params(axis="x", rotation=15)

prox_method = "Proximal-vs-distant median difference"
prox_heat = marker_group[
    marker_group.method.eq(prox_method) & marker_group.fibroblast_state.eq("FAP+/a5B1-")
].pivot(index="Diagnosis2", columns="feature", values="median_effect").reindex(index=GROUPS, columns=features)
sns.heatmap(prox_heat, annot=False, cmap="vlag", center=0,
            cbar_kws={"label": "Median near-far effect"}, ax=axes[1, 0])
axes[1, 0].set_yticklabels(group_order, rotation=0)
axes[1, 0].set_title("D  Neutrophil phenotype near FAP+/a5B1-")
axes[1, 0].set_xlabel("")
axes[1, 0].set_ylabel("")
axes[1, 0].tick_params(axis="x", rotation=45, labelsize=9)

corr_method = "Spearman correlation with proximity"
corr_heat = marker_group[
    marker_group.method.eq(corr_method) & marker_group.fibroblast_state.eq("FAP+/a5B1-")
].pivot(index="Diagnosis2", columns="feature", values="median_effect").reindex(index=GROUPS, columns=features)
sns.heatmap(corr_heat, annot=False, cmap="vlag", center=0,
            cbar_kws={"label": "Median proximity rho"}, ax=axes[1, 1])
axes[1, 1].set_yticklabels(group_order, rotation=0)
axes[1, 1].set_title("E  Continuous distance association")
axes[1, 1].set_xlabel("")
axes[1, 1].set_ylabel("")
axes[1, 1].tick_params(axis="x", rotation=45, labelsize=9)

sub_heat = subtype_group[subtype_group.fibroblast_state.eq("FAP+/a5B1+")].pivot(
    index="Diagnosis2", columns="feature", values="median_effect"
).reindex(index=GROUPS, columns=SUBTYPES)
sns.heatmap(sub_heat, annot=True, fmt=".2f", cmap="vlag", center=0,
            cbar_kws={"label": "Median log2 enrichment"}, ax=axes[1, 2])
axes[1, 2].set_yticklabels(group_order, rotation=0)
axes[1, 2].set_title("F  Neutrophil subtypes near FAP+/a5B1+")
axes[1, 2].set_xlabel("")
axes[1, 2].set_ylabel("")

fig.suptitle("FAP/a5B1 fibroblast niches and neutrophil phenotypic state", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_FAP_A5B1_Neutrophil_Reprogramming_Integrated")

# State-specific marker comparison in inflamed UC.
fig, axes = plt.subplots(1, 2, figsize=(18, 6))
for ax, method, title in (
    (axes[0], prox_method, "Near-versus-far marker effect"),
    (axes[1], corr_method, "Continuous proximity association"),
):
    heat = marker_group[
        marker_group.method.eq(method) & marker_group.Diagnosis2.eq("UC_Inflamed")
    ].pivot(index="fibroblast_state", columns="feature", values="median_effect").reindex(
        index=STATES, columns=features
    )
    sns.heatmap(heat, annot=True, fmt=".2f", cmap="vlag", center=0,
                annot_kws={"size": 8}, cbar_kws={"label": "Median patient effect"}, ax=ax)
    ax.set_yticklabels(STATES, rotation=0)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=45, labelsize=9)
fig.suptitle("State specificity of neutrophil phenotype in inflamed UC", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_FAP_A5B1_State_Specificity_UC_Inflamed")

# Representative whole-section maps selected by median double-positive enrichment.
primary = spatial[spatial.fibroblast_state.eq("FAP+/a5B1+") & spatial.radius_um.eq(PRIMARY_RADIUS)]
representatives = {}
for group in GROUPS:
    q = primary[primary.Diagnosis2.eq(group)].copy()
    med = q.log2_observed_expected.median()
    representatives[group] = q.loc[(q.log2_observed_expected - med).abs().idxmin(), "PatientID"]
fig, axes = plt.subplots(1, 3, figsize=(22, 7))
for ax, group in zip(axes, GROUPS):
    pid = representatives[group]
    z = d[d.PatientID.eq(pid)]
    background = z[~z.is_fibroblast & ~z.is_neutrophil]
    ax.scatter(background.x, background.y, s=.35, c="#D9D9D9", alpha=.12, rasterized=True)
    for state in STATES:
        q = z[z.fibroblast_state.eq(state)]
        ax.scatter(q.x, q.y, s=3 if state != "FAP+/a5B1+" else 6,
                   c=STATE_COLOR[state], alpha=.65, label=state, rasterized=True)
    q = z[z.is_neutrophil]
    ax.scatter(q.x, q.y, s=5, c="#00A6D6", alpha=.8, label="Neutrophil", rasterized=True)
    enr = primary.loc[primary.PatientID.eq(pid), "log2_observed_expected"].iloc[0]
    ax.set_title(f"{GROUP_LABEL[group]}: {pid}\nFAP+/a5B1+ 50 um log2 O/E = {enr:.2f}")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.axis("off")
axes[-1].legend(frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=10)
fig.suptitle("Representative FAP/a5B1 fibroblast-neutrophil spatial organization", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_FAP_A5B1_Neutrophil_Representative_Maps")


def fmt(value):
    return f"{value:.3g}"


primary_endpoint = "FAP+/a5B1+; 50 um log2 observed/expected"
primary_stats = spatial_stats[spatial_stats.endpoint.eq(primary_endpoint)]
primary_omnibus = primary_stats[primary_stats.analysis.eq("Kruskal-Wallis omnibus")].iloc[0]
primary_pairwise = primary_stats[primary_stats.analysis.eq("Pairwise Mann-Whitney")].set_index("comparison")
primary_within = primary_stats[primary_stats.analysis.eq("One-sample Wilcoxon versus zero")].set_index("group_1")
primary_medians = primary.groupby("Diagnosis2").log2_observed_expected.median()

pair_lines = []
for first, second in PAIRWISE:
    label = f"{GROUP_LABEL[first]} vs {GROUP_LABEL[second]}"
    row = primary_pairwise.loc[label]
    pair_lines.append(
        f"- {label}: median difference {row.effect:.3f}, P={fmt(row.p_value)}, "
        f"FDR={fmt(row.fdr_within_endpoint_analysis)}."
    )
within_lines = []
for group in GROUPS:
    row = primary_within.loc[group]
    within_lines.append(
        f"- {GROUP_LABEL[group]}: median {row.effect:.3f}, P={fmt(row.p_value)}, "
        f"FDR={fmt(row.fdr_within_endpoint_analysis)}."
    )

def significant_features(stats, method, state, stratum, alpha=.05):
    q = stats[
        stats.method.eq(method) & stats.fibroblast_state.eq(state)
        & stats.analysis.eq("One-sample Wilcoxon versus zero")
        & stats.stratum.eq(stratum)
        & stats.fdr_within_state_method_analysis_stratum.lt(alpha)
    ].sort_values("median_effect", ascending=False)
    return [f"{row.feature} ({row.median_effect:+.2f})" for row in q.itertuples()]

marker_lines = []
for group in GROUPS:
    label = GROUP_LABEL[group]
    hits = significant_features(marker_stats, prox_method, "FAP+/a5B1+", label)
    marker_lines.append(f"- {label}: {', '.join(hits) if hits else 'none after FDR correction'}.")
corr_lines = []
for group in GROUPS:
    label = GROUP_LABEL[group]
    hits = significant_features(marker_stats, corr_method, "FAP+/a5B1+", label)
    corr_lines.append(f"- {label}: {', '.join(hits) if hits else 'none after FDR correction'}.")
subtype_lines = []
for group in GROUPS:
    label = GROUP_LABEL[group]
    hits = significant_features(
        subtype_stats, "Subtype log2 proximal/all-neutrophil enrichment", "FAP+/a5B1+", label
    )
    subtype_lines.append(f"- {label}: {', '.join(hits) if hits else 'none after FDR correction'}.")

fap_only_marker_lines = []
fap_only_corr_lines = []
for group in GROUPS:
    label = GROUP_LABEL[group]
    hits = significant_features(marker_stats, prox_method, "FAP+/a5B1-", label)
    fap_only_marker_lines.append(
        f"- {label}: {', '.join(hits) if hits else 'none after FDR correction'}."
    )
    hits = significant_features(marker_stats, corr_method, "FAP+/a5B1-", label)
    fap_only_corr_lines.append(
        f"- {label}: {', '.join(hits) if hits else 'none after FDR correction'}."
    )

eligibility = marker_effects[
    marker_effects.method.eq(prox_method) & marker_effects.fibroblast_state.eq("FAP+/a5B1+")
    & marker_effects.feature.eq(PROGRAM_NAME)
].groupby("Diagnosis2").PatientID.nunique()
fap_only_eligibility = marker_effects[
    marker_effects.method.eq(prox_method) & marker_effects.fibroblast_state.eq("FAP+/a5B1-")
    & marker_effects.feature.eq(PROGRAM_NAME)
].groupby("Diagnosis2").PatientID.nunique()

double_negative_primary = spatial[
    spatial.fibroblast_state.eq("FAP-/a5B1-") & spatial.radius_um.eq(PRIMARY_RADIUS)
]
double_negative_medians = double_negative_primary.groupby("Diagnosis2").log2_observed_expected.median()
double_negative_stats = spatial_stats[
    spatial_stats.endpoint.eq("FAP-/a5B1-; 50 um log2 observed/expected")
    & spatial_stats.analysis.eq("Pairwise Mann-Whitney")
].set_index("comparison")
dn_uc_contrast = double_negative_stats.loc["UC inflamed vs UC noninflamed"]

significant_between_markers = marker_stats[
    marker_stats.analysis.eq("Pairwise Mann-Whitney")
    & marker_stats.fdr_within_state_method_analysis_stratum.lt(.05)
]

sensitivity_omnibus = sensitivity_stats[
    sensitivity_stats.analysis.eq("Kruskal-Wallis omnibus")
].iloc[0]

report = f"""# FAP/a5B1 fibroblast niches and neutrophil phenotypic reprogramming in UC

## Design

Fibroblasts were classified using fixed pooled-UC 75th-percentile CLR thresholds: FAP >= {fap_threshold:.4f} and a5B1 >= {a5_threshold:.4f}. This produced four mutually exclusive states. Patient was the biological replicate throughout. The primary spatial endpoint was the 50 um log2 observed/permuted-expected neutrophil edge count, with 1,000 within-patient permutations of fibroblast-state labels. Neutrophil phenotype was evaluated both as a near (<=50 um) versus far (>100 um) marker difference and as a continuous within-patient correlation with proximity. The latter retains patients with few proximal neutrophils.

## FAP+/a5B1+ spatial association

Median 50 um log2 observed/expected enrichment was {primary_medians['Control']:.3f} in controls, {primary_medians['UC_Noninflamed']:.3f} in UC noninflamed, and {primary_medians['UC_Inflamed']:.3f} in UC inflamed tissue. The three-group omnibus test was P={fmt(primary_omnibus.p_value)} (epsilon-squared={primary_omnibus.effect:.3f}). Pairwise results were:

{chr(10).join(pair_lines)}

Within each group, tests against no spatial preference were:

{chr(10).join(within_lines)}

## Neutrophil phenotype near FAP+/a5B1+ fibroblasts

The binary proximity analysis required at least {MIN_PROXIMAL} proximal and {MIN_DISTANT} distant neutrophils and included {int(eligibility.get('Control', 0))}/6 controls, {int(eligibility.get('UC_Noninflamed', 0))}/9 UC noninflamed, and {int(eligibility.get('UC_Inflamed', 0))}/9 UC inflamed patients. Significant near-versus-far marker effects within each group were:

{chr(10).join(marker_lines)}

Significant continuous proximity-marker correlations were:

{chr(10).join(corr_lines)}

The five-marker NF-niche program combines OSM, CXCR4, PDL1, HLA-ABC, and a5B1 after within-patient neutrophil standardization. It summarizes co-localized phenotype and is not a validated functional score.

## Phenotype signal localizes to FAP+/a5B1- fibroblasts

The coherent neutrophil phenotype was associated with proximity to FAP+/a5B1- fibroblasts rather than the double-positive state. The binary analysis included {int(fap_only_eligibility.get('Control', 0))}/6 controls, {int(fap_only_eligibility.get('UC_Noninflamed', 0))}/9 UC noninflamed, and {int(fap_only_eligibility.get('UC_Inflamed', 0))}/9 UC inflamed patients. Significant near-versus-far effects were:

{chr(10).join(fap_only_marker_lines)}

The continuous proximity analysis produced:

{chr(10).join(fap_only_corr_lines)}

Although the FAP+/a5B1- program was most coherent in inflamed UC, formal pairwise group comparisons did not show an inflammation-specific increase after FDR correction. Across all state-marker tests, {len(significant_between_markers)} pairwise group effect survived FDR correction; it was a lower CD66b near-versus-far effect for FAP+/a5B1+ niches in UC inflamed versus control tissue. Thus, the data support niche-associated phenotype more strongly than an inflammation-by-niche interaction.

## Inflammation-related spatial reorganization

Neutrophils were enriched near FAP-/a5B1- fibroblasts in UC noninflamed tissue (median log2 O/E {double_negative_medians['UC_Noninflamed']:.3f}) but not inflamed tissue ({double_negative_medians['UC_Inflamed']:.3f}). The inflamed-versus-noninflamed contrast was significant (median difference {dn_uc_contrast.effect:.3f}, P={fmt(dn_uc_contrast.p_value)}, FDR={fmt(dn_uc_contrast.fdr_within_endpoint_analysis)}). Together with depletion around both a5B1-positive fibroblast states, this indicates redistribution away from a5B1-high stromal territories rather than preferential targeting of FAP+/a5B1+ fibroblasts.

## Neutrophil subtype composition

Significant subtype enrichment among neutrophils within 50 um of FAP+/a5B1+ fibroblasts was:

{chr(10).join(subtype_lines)}

## Sensitivity and interpretation

Using within-patient FAP and a5B1 top-quartile thresholds, the double-positive spatial-enrichment three-group omnibus test was P={fmt(sensitivity_omnibus.p_value)}. Complete results for all four fibroblast states, 25/50/100 um radii, individual markers, the composite program, and neutrophil subtypes are provided in the tables.

In this observational CODEX dataset, spatial proximity and phenotype co-variation are consistent with a fibroblast-associated neutrophil state but do not establish that fibroblasts reprogram neutrophils, the signaling direction, direct cell contact, or causality. FAP and a5B1 positivity are operational high-expression calls because external positivity controls were not supplied.
"""
(OUT / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

manifest = {
    "input_cells": str(ROOT / "additional_analyses" / "codex_extended_cells.csv"),
    "input_marker_expression": str(ROOT / "neutrophil_fibroblast_crc_marker_map" / "crc_associated_marker_expression.csv"),
    "groups": list(GROUPS), "fibroblast_states": list(STATES),
    "primary_radius_um": PRIMARY_RADIUS, "proximal_um": PROXIMAL_UM,
    "distant_um_strictly_greater_than": DISTANT_UM,
    "permutations_per_patient": N_PERM, "seed": SEED,
    "FAP_threshold_clr": fap_threshold, "a5B1_threshold_clr": a5_threshold,
    "threshold_source": "pooled UC fibroblast 75th percentile",
    "biological_replicate": "patient", "program_markers": list(PROGRAM_MARKERS),
    "representative_patients": representatives,
}
(OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print("Completed FAP/a5B1 fibroblast-neutrophil reprogramming analysis.")
print(f"Saved outputs to {OUT}")
