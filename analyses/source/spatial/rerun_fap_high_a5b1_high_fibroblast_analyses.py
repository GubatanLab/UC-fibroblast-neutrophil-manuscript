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
OUT = ROOT / "fap_high_a5b1_high_rerun"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GROUPS = ("UC_Noninflamed", "UC_Inflamed")
GROUP_LABEL = {"UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}
GROUP_COLOR = {"UC_Noninflamed": "#D4AD35", "UC_Inflamed": "#D95F5F"}
TARGETS = ("FAP-high", "a5B1-high", "FAP-high/a5B1-high")
RADII = (25, 50, 100)
PRIMARY_RADIUS = 50
N_FIBROBLASTS = 250
N_NEUTROPHILS = 10
N_PROXIMAL = 3
N_DISTANT = 3
N_BOOT = 200
SEED = 20260806
FDR_ALPHA = 0.10
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


def uc_tests(frame, value, endpoint, zero_tests=True):
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
                    "endpoint": endpoint, "analysis": "Effect versus zero",
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
crc = pd.read_csv(ROOT / "neutrophil_fibroblast_crc_marker_map" / "crc_associated_marker_expression.csv")
d = cells.merge(crc, on="cell_ID", how="left", validate="one_to_one")
d = d[d.Diagnosis2.isin(GROUPS)].copy()
d["is_fibroblast"] = d.cell_type.eq("Fibroblast")
d["is_neutrophil"] = d.cell_type.str.startswith("Neutrophil", na=False)

# High is defined separately within every patient using fibroblast CLR q75.
fib = d[d.is_fibroblast].copy()
thresholds = fib.groupby(["PatientID", "Diagnosis2"]).agg(
    FAP_q75_clr=("FAPa_cell_clr", lambda x: x.quantile(.75)),
    a5B1_q75_clr=("a5B1_cell_clr", lambda x: x.quantile(.75)),
    n_fibroblasts=("cell_ID", "size"),
).reset_index()
thresholds.to_csv(TAB / "00_patient_high_thresholds.csv", index=False)
d = d.merge(thresholds[["PatientID", "FAP_q75_clr", "a5B1_q75_clr"]], on="PatientID", how="left", validate="many_to_one")
d["FAP-high"] = d.is_fibroblast & d.FAPa_cell_clr.ge(d.FAP_q75_clr)
d["a5B1-high"] = d.is_fibroblast & d.a5B1_cell_clr.ge(d.a5B1_q75_clr)
d["FAP-high/a5B1-high"] = d["FAP-high"] & d["a5B1-high"]

# Target abundance.
abundance_rows = []
for (pid, diagnosis), z in d.groupby(["PatientID", "Diagnosis2"]):
    f = z[z.is_fibroblast]
    for target in TARGETS:
        abundance_rows.append({
            "PatientID": pid, "Diagnosis2": diagnosis, "target": target,
            "n_fibroblasts": len(f), "n_target": int(f[target].sum()),
            "fraction_fibroblasts": float(f[target].mean()),
        })
abundance = pd.DataFrame(abundance_rows)
abundance.to_csv(TAB / "01_high_target_abundance_by_patient.csv", index=False)
abundance_stats_rows = []
for target, q in abundance.groupby("target"):
    abundance_stats_rows.extend(uc_tests(q, "fraction_fibroblasts", target, zero_tests=False))
abundance_stats = pd.DataFrame(abundance_stats_rows)
abundance_stats["fdr"] = bh_adjust(abundance_stats.p_value)
abundance_stats.to_csv(TAB / "02_high_target_abundance_statistics.csv", index=False)

# Raw multiscale edge enrichment and exact distances.
spatial_rows = []
distance_rows = []
for pid, z in d.groupby("PatientID", sort=True):
    f = z[z.is_fibroblast].reset_index(drop=True)
    n = z[z.is_neutrophil].reset_index(drop=True)
    if f.empty or n.empty:
        continue
    fib_xy = f[["x", "y"]].to_numpy()
    neut_xy = n[["x", "y"]].to_numpy()
    neighbors = cKDTree(fib_xy).query_ball_point(neut_xy, max(RADII))
    for target in TARGETS:
        mask = f[target].to_numpy()
        target_xy = fib_xy[mask]
        distances = cKDTree(target_xy).query(neut_xy)[0]
        distance_rows.extend({
            "cell_ID": cid, "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
            "target": target, "nearest_target_distance_um": dist,
        } for cid, dist in zip(n.cell_ID, distances))
    for radius in RADII:
        degrees = np.zeros(len(f), dtype=np.int64)
        for point, ids in zip(neut_xy, neighbors):
            if not ids:
                continue
            ids = np.asarray(ids, dtype=int)
            delta = fib_xy[ids] - point
            ids = ids[(delta * delta).sum(axis=1) <= radius * radius]
            if len(ids):
                np.add.at(degrees, ids, 1)
        total_edges = int(degrees.sum())
        for target in TARGETS:
            mask = f[target].to_numpy()
            n_target = int(mask.sum())
            observed = int(degrees[mask].sum())
            expected = total_edges * n_target / len(f)
            spatial_rows.append({
                "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
                "target": target, "radius_um": radius,
                "n_neutrophils": len(n), "n_fibroblasts": len(f),
                "n_target": n_target, "total_edges": total_edges,
                "observed_target_edges": observed, "expected_target_edges": expected,
                "log2_observed_expected": np.log2((observed + .5) / (expected + .5)),
            })
spatial = pd.DataFrame(spatial_rows)
distances = pd.DataFrame(distance_rows)
spatial.to_csv(TAB / "03_multiscale_spatial_enrichment_by_patient.csv", index=False)
distances.to_csv(TAB / "04_neutrophil_nearest_high_target_distances.csv", index=False)

spatial_stats_rows = []
for (target, radius), q in spatial.groupby(["target", "radius_um"]):
    spatial_stats_rows.extend(uc_tests(
        q, "log2_observed_expected", f"{target}; {int(radius)} um log2 observed/expected"
    ))
spatial_stats = pd.DataFrame(spatial_stats_rows)
spatial_stats["fdr_within_analysis"] = spatial_stats.groupby("analysis").p_value.transform(bh_adjust)
spatial_stats.to_csv(TAB / "05_multiscale_spatial_statistics.csv", index=False)

raw_distance = distances.groupby(["PatientID", "Diagnosis2", "target"]).agg(
    n_neutrophils=("cell_ID", "size"),
    median_nearest_target_um=("nearest_target_distance_um", "median"),
    q25_nearest_target_um=("nearest_target_distance_um", lambda x: x.quantile(.25)),
    q75_nearest_target_um=("nearest_target_distance_um", lambda x: x.quantile(.75)),
).reset_index()
for radius in RADII:
    frac = distances.assign(hit=distances.nearest_target_distance_um.le(radius)).groupby(
        ["PatientID", "Diagnosis2", "target"]
    ).hit.mean().rename(f"fraction_within_{radius}um").reset_index()
    raw_distance = raw_distance.merge(frac, on=["PatientID", "Diagnosis2", "target"])
raw_distance.to_csv(TAB / "06_raw_distance_summary_by_patient.csv", index=False)

raw_distance_stats_rows = []
for target, q in raw_distance.groupby("target"):
    raw_distance_stats_rows.extend(uc_tests(q, "median_nearest_target_um", f"{target}; raw median nearest distance", zero_tests=False))
    for radius in RADII:
        raw_distance_stats_rows.extend(uc_tests(q, f"fraction_within_{radius}um", f"{target}; fraction within {radius} um", zero_tests=False))
raw_distance_stats = pd.DataFrame(raw_distance_stats_rows)
raw_distance_stats["fdr"] = bh_adjust(raw_distance_stats.p_value)
raw_distance_stats.to_csv(TAB / "07_raw_distance_statistics.csv", index=False)

# Equal-count edge bootstrap plus abundance-adjusted proximity null.
rng = np.random.default_rng(SEED)
edge_boot_rows = []
proximity_boot_rows = []
eligibility_rows = []
for pid, z in d.groupby("PatientID", sort=True):
    f = z[z.is_fibroblast].reset_index(drop=True)
    n = z[z.is_neutrophil].reset_index(drop=True)
    eligible = len(f) >= N_FIBROBLASTS and len(n) >= N_NEUTROPHILS
    eligibility_rows.append({
        "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
        "n_fibroblasts": len(f), "n_neutrophils": len(n), "eligible": eligible,
    })
    if not eligible:
        continue
    fib_xy_all = f[["x", "y"]].to_numpy()
    neut_xy_all = n[["x", "y"]].to_numpy()
    observed_full_distances = {
        target: cKDTree(fib_xy_all[f[target].to_numpy()]).query(neut_xy_all)[0]
        for target in TARGETS
    }
    for iteration in range(1, N_BOOT + 1):
        fi = rng.choice(len(f), N_FIBROBLASTS, replace=False)
        ni = rng.choice(len(n), N_NEUTROPHILS, replace=False)
        fs = f.iloc[fi].reset_index(drop=True)
        ns_xy = neut_xy_all[ni]
        fs_xy = fs[["x", "y"]].to_numpy()
        degrees = np.zeros(N_FIBROBLASTS, dtype=np.int64)
        for ids in cKDTree(fs_xy).query_ball_point(ns_xy, PRIMARY_RADIUS):
            if ids:
                np.add.at(degrees, np.asarray(ids, dtype=int), 1)
        total_edges = int(degrees.sum())
        for target in TARGETS:
            mask = fs[target].to_numpy()
            n_target_sampled = int(mask.sum())
            if n_target_sampled:
                observed = int(degrees[mask].sum())
                expected = total_edges * n_target_sampled / N_FIBROBLASTS
                edge_boot_rows.append({
                    "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
                    "iteration": iteration, "target": target,
                    "sampled_fibroblasts": N_FIBROBLASTS,
                    "sampled_neutrophils": N_NEUTROPHILS,
                    "sampled_target_fibroblasts": n_target_sampled,
                    "observed_edges": observed, "expected_edges": expected,
                })
            full_mask = f[target].to_numpy()
            n_target_full = int(full_mask.sum())
            observed_dist = observed_full_distances[target][ni]
            random_idx = rng.choice(len(f), n_target_full, replace=False)
            expected_dist = cKDTree(fib_xy_all[random_idx]).query(ns_xy)[0]
            proximity_boot_rows.append({
                "PatientID": pid, "Diagnosis2": z.Diagnosis2.iloc[0],
                "iteration": iteration, "target": target,
                "observed_median_distance_um": np.median(observed_dist),
                "expected_median_random_distance_um": np.median(expected_dist),
                "log2_observed_random_distance": np.log2(
                    (np.median(observed_dist) + .5) / (np.median(expected_dist) + .5)
                ),
                "observed_fraction_within_50um": np.mean(observed_dist <= 50),
                "expected_fraction_within_50um": np.mean(expected_dist <= 50),
                "difference_fraction_within_50um": np.mean(observed_dist <= 50) - np.mean(expected_dist <= 50),
            })

eligibility = pd.DataFrame(eligibility_rows)
eligibility.to_csv(TAB / "08_equal_count_eligibility.csv", index=False)
edge_boot = pd.DataFrame(edge_boot_rows)
proximity_boot = pd.DataFrame(proximity_boot_rows)
edge_boot.to_csv(TAB / "09_equal_count_edge_bootstrap_iterations.csv", index=False)
proximity_boot.to_csv(TAB / "10_abundance_adjusted_proximity_bootstrap_iterations.csv", index=False)

edge_patient = edge_boot.groupby(["PatientID", "Diagnosis2", "target"]).agg(
    n_iterations=("iteration", "nunique"),
    summed_observed_edges=("observed_edges", "sum"),
    summed_expected_edges=("expected_edges", "sum"),
    median_sampled_target=("sampled_target_fibroblasts", "median"),
).reset_index()
edge_patient["aggregate_log2_observed_expected"] = np.log2(
    (edge_patient.summed_observed_edges + .5) / (edge_patient.summed_expected_edges + .5)
)
edge_patient.to_csv(TAB / "11_equal_count_spatial_by_patient.csv", index=False)
edge_stats_rows = []
for target, q in edge_patient.groupby("target"):
    edge_stats_rows.extend(uc_tests(q, "aggregate_log2_observed_expected", f"{target}; equal-count 50 um log2 O/E"))
edge_stats = pd.DataFrame(edge_stats_rows)
edge_stats["fdr_within_analysis"] = edge_stats.groupby("analysis").p_value.transform(bh_adjust)
edge_stats.to_csv(TAB / "12_equal_count_spatial_statistics.csv", index=False)

proximity_patient = proximity_boot.groupby(["PatientID", "Diagnosis2", "target"]).agg(
    n_iterations=("iteration", "nunique"),
    median_log2_observed_random_distance=("log2_observed_random_distance", "median"),
    median_difference_fraction_within_50um=("difference_fraction_within_50um", "median"),
).reset_index()
proximity_patient.to_csv(TAB / "13_abundance_adjusted_proximity_by_patient.csv", index=False)
proximity_stats_rows = []
for target, q in proximity_patient.groupby("target"):
    proximity_stats_rows.extend(uc_tests(q, "median_log2_observed_random_distance", f"{target}; normalized nearest distance"))
    proximity_stats_rows.extend(uc_tests(q, "median_difference_fraction_within_50um", f"{target}; normalized 50 um contact fraction"))
proximity_stats = pd.DataFrame(proximity_stats_rows)
proximity_stats["fdr_within_analysis"] = proximity_stats.groupby("analysis").p_value.transform(bh_adjust)
proximity_stats.to_csv(TAB / "14_abundance_adjusted_proximity_statistics.csv", index=False)

# Equal-count proximal-versus-distant neutrophil marker analysis.
marker_source = {
    "OSM": "raw_OSM", "CXCR4": "raw_CXCR4", "PDL1": "raw_PDL1",
    "HLA_ABC": "raw_HLA_ABC", "a5B1": "raw_a5B1", "PADI4": "PADI4",
    "MX1": "MX.1", "CD66b": "CD66b", "CD16": "CD16", "CD11b": "CD11b",
    "Ki67": "raw_Ki67",
}
neut = d[d.is_neutrophil].copy()
wide_dist = distances.pivot(index="cell_ID", columns="target", values="nearest_target_distance_um")
wide_dist = wide_dist.rename(columns={target: f"distance__{target}" for target in TARGETS})
neut = neut.merge(wide_dist, left_on="cell_ID", right_index=True, how="left", validate="one_to_one")
for marker, source in marker_source.items():
    neut[f"log2_{marker}"] = np.log2(neut[source].astype(float).clip(lower=0) + .1)
for marker in PROGRAM_MARKERS:
    col = f"log2_{marker}"
    neut[f"z_{marker}"] = neut.groupby("PatientID")[col].transform(
        lambda x: (x - x.mean()) / x.std(ddof=0) if x.std(ddof=0) > 0 else 0.0
    )
neut[PROGRAM_NAME] = neut[[f"z_{m}" for m in PROGRAM_MARKERS]].mean(axis=1)

rng_marker = np.random.default_rng(SEED + 1)
features = list(MARKERS) + [PROGRAM_NAME]
marker_rows = []
for (pid, diagnosis), q in neut.groupby(["PatientID", "Diagnosis2"]):
    for target in TARGETS:
        ranked = q.sort_values(f"distance__{target}")
        pool_size = max(N_PROXIMAL, int(np.ceil(len(ranked) / 3)))
        prox = ranked.head(pool_size)
        far = ranked.tail(pool_size)
        pi = rng_marker.integers(0, len(prox), size=(N_BOOT, N_PROXIMAL))
        di = rng_marker.integers(0, len(far), size=(N_BOOT, N_DISTANT))
        for feature in features:
            p = (prox[PROGRAM_NAME] if feature == PROGRAM_NAME else prox[f"log2_{feature}"]).to_numpy()
            f = (far[PROGRAM_NAME] if feature == PROGRAM_NAME else far[f"log2_{feature}"]).to_numpy()
            effects = np.median(p[pi], axis=1) - np.median(f[di], axis=1)
            marker_rows.append({
                "PatientID": pid, "Diagnosis2": diagnosis, "target": target,
                "feature": feature, "n_proximal": len(prox), "n_distant": len(far),
                "bootstrap_median_effect": np.median(effects),
                "bootstrap_ci025": np.quantile(effects, .025),
                "bootstrap_ci975": np.quantile(effects, .975),
            })
marker_effects = pd.DataFrame(marker_rows)
marker_effects.to_csv(TAB / "15_equal_count_marker_effects_by_patient.csv", index=False)

marker_stats_rows = []
for (target, feature), q in marker_effects.groupby(["target", "feature"]):
    marker_stats_rows.extend(uc_tests(q, "bootstrap_median_effect", f"{target}; {feature}; equal near-far effect"))
marker_stats = pd.DataFrame(marker_stats_rows)
marker_stats["target"] = marker_stats.endpoint.str.split(";").str[0]
marker_stats["feature"] = marker_stats.endpoint.str.split(";").str[1].str.strip()
marker_stats["fdr_within_target_analysis"] = marker_stats.groupby(["target", "analysis"]).p_value.transform(bh_adjust)
marker_stats.to_csv(TAB / "16_equal_count_marker_statistics.csv", index=False)

# Figure.
sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})
order = [GROUP_LABEL[g] for g in GROUPS]
palette = {GROUP_LABEL[g]: GROUP_COLOR[g] for g in GROUPS}
fig, axes = plt.subplots(2, 3, figsize=(21, 13))

q = abundance.copy(); q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="target", y="fraction_fibroblasts", hue="Group", hue_order=order,
            palette=palette, showfliers=False, ax=axes[0, 0])
axes[0, 0].set_title("A  Within-patient high-state abundance")
axes[0, 0].set_xlabel(""); axes[0, 0].set_ylabel("Fraction of fibroblasts")
axes[0, 0].tick_params(axis="x", rotation=25, labelsize=9)
axes[0, 0].legend(title="", frameon=False, fontsize=9)

multi = spatial.groupby(["target", "radius_um", "Diagnosis2"]).log2_observed_expected.median().unstack("Diagnosis2")
multi["difference"] = multi.UC_Inflamed - multi.UC_Noninflamed
multi = multi.difference.unstack("radius_um").reindex(index=TARGETS)
sns.heatmap(multi, annot=True, fmt=".2f", cmap="vlag", center=0,
            cbar_kws={"label": "Inflamed - noninflamed"}, ax=axes[0, 1])
axes[0, 1].set_title("B  Neutrophil proximity by radius")
axes[0, 1].set_xlabel("Radius (um)"); axes[0, 1].set_ylabel("")
axes[0, 1].tick_params(axis="y", rotation=0, labelsize=9)

q = edge_patient.copy(); q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="target", y="aggregate_log2_observed_expected", hue="Group",
            hue_order=order, palette=palette, showfliers=False, ax=axes[0, 2])
axes[0, 2].axhline(0, color="black", linestyle="--", linewidth=1)
axes[0, 2].set_title("C  Neutrophils within 50 um of fibroblasts")
axes[0, 2].set_xlabel(""); axes[0, 2].set_ylabel("log2 observed / expected proximity")
axes[0, 2].tick_params(axis="x", rotation=25, labelsize=9)
axes[0, 2].legend(title="", frameon=False, fontsize=9)

q = proximity_patient.copy(); q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(data=q, x="target", y="median_log2_observed_random_distance", hue="Group",
            hue_order=order, palette=palette, showfliers=False, ax=axes[1, 0])
axes[1, 0].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1, 0].set_title("D  Nearest neutrophil-to-fibroblast distance")
axes[1, 0].set_xlabel(""); axes[1, 0].set_ylabel("log2 observed / random distance")
axes[1, 0].tick_params(axis="x", rotation=25, labelsize=9)
axes[1, 0].legend(title="", frameon=False, fontsize=9)

contact = proximity_patient.groupby(["Diagnosis2", "target"]).median_difference_fraction_within_50um.median().unstack().reindex(index=GROUPS, columns=TARGETS)
sns.heatmap(contact, annot=True, fmt=".2f", cmap="vlag", center=0,
            cbar_kws={"label": "Observed - random fraction"}, ax=axes[1, 1])
axes[1, 1].set_yticklabels(order, rotation=0)
axes[1, 1].set_title("E  Neutrophils within 50 um of fibroblasts")
axes[1, 1].set_xlabel(""); axes[1, 1].set_ylabel("")

heat = marker_effects.groupby(["Diagnosis2", "target", "feature"]).bootstrap_median_effect.median().reset_index()
heat["row"] = heat.Diagnosis2.map(GROUP_LABEL) + " | " + heat.target
heat = heat.pivot(index="row", columns="feature", values="bootstrap_median_effect")
row_order = [f"{GROUP_LABEL[g]} | {t}" for g in GROUPS for t in TARGETS]
heat = heat.reindex(index=row_order, columns=features)
sns.heatmap(heat, cmap="vlag", center=0, cbar_kws={"label": "Closest - farthest third"}, ax=axes[1, 2])
axes[1, 2].set_title("F  Neutrophil state by fibroblast proximity")
axes[1, 2].set_xlabel(""); axes[1, 2].set_ylabel("")
axes[1, 2].tick_params(axis="x", rotation=45, labelsize=8)
axes[1, 2].tick_params(axis="y", labelsize=8)

fig.suptitle("Neutrophil proximity to FAP-high and a5B1-high fibroblasts in UC", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_FAP_High_a5B1_High_Rerun")


def fmt(x):
    return f"{x:.3g}"


edge_between = edge_stats[edge_stats.analysis.eq("UC inflamed versus UC noninflamed")].copy()
edge_between["target"] = edge_between.endpoint.str.split(";").str[0]
edge_between = edge_between.set_index("target")
edge_within = edge_stats[edge_stats.analysis.eq("Effect versus zero")].copy()
edge_within["target"] = edge_within.endpoint.str.split(";").str[0]
edge_medians = edge_patient.groupby(["target", "Diagnosis2"]).aggregate_log2_observed_expected.median()
edge_lines = []
within_lines = []
for target in TARGETS:
    row = edge_between.loc[target]
    edge_lines.append(
        f"- {target}: UC noninflamed {edge_medians.loc[(target, 'UC_Noninflamed')]:.3f}; "
        f"UC inflamed {edge_medians.loc[(target, 'UC_Inflamed')]:.3f}; difference "
        f"{row.effect_inflamed_minus_noninflamed:.3f}, P={fmt(row.p_value)}, "
        f"FDR={fmt(row.fdr_within_analysis)}."
    )
    for group in GROUPS:
        r = edge_within[edge_within.target.eq(target) & edge_within.group.eq(group)].iloc[0]
        within_lines.append(
            f"- {target}, {GROUP_LABEL[group]}: median {r.effect_inflamed_minus_noninflamed:.3f}, "
            f"P={fmt(r.p_value)}, FDR={fmt(r.fdr_within_analysis)}."
        )

prox_between = proximity_stats[
    proximity_stats.analysis.eq("UC inflamed versus UC noninflamed")
    & proximity_stats.endpoint.str.contains("normalized nearest distance", regex=False)
].copy()
prox_between["target"] = prox_between.endpoint.str.split(";").str[0]
prox_between = prox_between.set_index("target")
prox_within = proximity_stats[
    proximity_stats.analysis.eq("Effect versus zero")
    & proximity_stats.endpoint.str.contains("normalized nearest distance", regex=False)
].copy()
prox_within["target"] = prox_within.endpoint.str.split(";").str[0]
prox_medians = proximity_patient.groupby(["target", "Diagnosis2"]).median_log2_observed_random_distance.median()
prox_lines = []
for target in TARGETS:
    r = prox_between.loc[target]
    prox_lines.append(
        f"- {target}: UC noninflamed {prox_medians.loc[(target, 'UC_Noninflamed')]:.3f}; "
        f"UC inflamed {prox_medians.loc[(target, 'UC_Inflamed')]:.3f}; difference "
        f"{r.effect_inflamed_minus_noninflamed:.3f}, P={fmt(r.p_value)}, "
        f"FDR={fmt(r.fdr_within_analysis)}."
    )

prox_within_lines = []
for target in TARGETS:
    for group in GROUPS:
        r = prox_within[prox_within.target.eq(target) & prox_within.group.eq(group)].iloc[0]
        prox_within_lines.append(
            f"- {target}, {GROUP_LABEL[group]}: normalized distance "
            f"{r.effect_inflamed_minus_noninflamed:.3f}, P={fmt(r.p_value)}, "
            f"FDR={fmt(r.fdr_within_analysis)}."
        )

contact_within = proximity_stats[
    proximity_stats.analysis.eq("Effect versus zero")
    & proximity_stats.endpoint.str.contains("normalized 50 um contact fraction", regex=False)
].copy()
contact_within["target"] = contact_within.endpoint.str.split(";").str[0]
contact_lines = []
for target in TARGETS:
    for group in GROUPS:
        r = contact_within[contact_within.target.eq(target) & contact_within.group.eq(group)].iloc[0]
        contact_lines.append(
            f"- {target}, {GROUP_LABEL[group]}: observed-minus-random fraction "
            f"{r.effect_inflamed_minus_noninflamed:.3f}, P={fmt(r.p_value)}, "
            f"FDR={fmt(r.fdr_within_analysis)}."
        )

abundance_between = abundance_stats[
    abundance_stats.analysis.eq("UC inflamed versus UC noninflamed")
].set_index("endpoint")
abundance_lines = []
for target in TARGETS:
    r = abundance_between.loc[target]
    abundance_lines.append(
        f"- {target}: UC noninflamed {r.median_UC_noninflamed:.3f}; "
        f"UC inflamed {r.median_UC_inflamed:.3f}; P={fmt(r.p_value)}, FDR={fmt(r.fdr)}."
    )

marker_between = marker_stats[
    marker_stats.analysis.eq("UC inflamed versus UC noninflamed")
    & marker_stats.fdr_within_target_analysis.lt(FDR_ALPHA)
]
marker_within = marker_stats[
    marker_stats.analysis.eq("Effect versus zero")
    & marker_stats.fdr_within_target_analysis.lt(FDR_ALPHA)
]
marker_lines = []
for target in TARGETS:
    for group in GROUPS:
        q = marker_within[marker_within.target.eq(target) & marker_within.group.eq(group)].sort_values(
            "effect_inflamed_minus_noninflamed", ascending=False
        )
        hits = [
            f"{r.feature} ({r.effect_inflamed_minus_noninflamed:+.2f}; "
            f"P={fmt(r.p_value)}, FDR={fmt(r.fdr_within_target_analysis)})"
            for r in q.itertuples()
        ]
        marker_lines.append(
            f"- {target}, {GROUP_LABEL[group]}: "
            f"{', '.join(hits) if hits else f'none at FDR < {FDR_ALPHA:.2f}'}."
        )
between_lines = []
for target in TARGETS:
    q = marker_between[marker_between.target.eq(target)].sort_values("effect_inflamed_minus_noninflamed", ascending=False)
    hits = [
        f"{r.feature} ({r.effect_inflamed_minus_noninflamed:+.2f}; "
        f"P={fmt(r.p_value)}, FDR={fmt(r.fdr_within_target_analysis)})"
        for r in q.itertuples()
    ]
    between_lines.append(
        f"- {target}: "
        f"{', '.join(hits) if hits else f'no UC-state marker difference at FDR < {FDR_ALPHA:.2f}'}."
    )

excluded = eligibility[~eligibility.eligible]
excluded_text = ", ".join(f"{r.PatientID} ({int(r.n_neutrophils)} neutrophils)" for r in excluded.itertuples()) or "none"

report = f"""# Neutrophil proximity to FAP-high and a5B1-high fibroblasts in UC

## Definitions and normalization

The comparison was restricted to UC noninflamed and UC inflamed samples. FAP-high and a5B1-high were defined as the top fibroblast CLR quartile separately within each UC patient; double-high required both calls. Proximity analyses used {N_BOOT} iterations with exactly {N_FIBROBLASTS} fibroblasts and {N_NEUTROPHILS} neutrophils per sample, the UC-cohort minimum, so every UC patient was retained. Nearest-distance normalization compared each high-state fibroblast population with an equal-sized random fibroblast subset from the same patient. Marker comparisons sampled exactly {N_PROXIMAL} neutrophils from the closest third and {N_DISTANT} from the farthest third within each patient. Patient was the biological replicate. Patient exclusions: {excluded_text}. Raw two-sided P values are reported alongside Benjamini-Hochberg adjusted values; the prespecified discovery threshold is FDR < {FDR_ALPHA:.2f}.

## Neutrophil proximity within 50 um of high-state fibroblasts

{chr(10).join(edge_lines)}

Within-group proximity effects relative to random fibroblast labeling were:

{chr(10).join(within_lines)}

## High-state abundance

{chr(10).join(abundance_lines)}

## Nearest neutrophil-to-fibroblast distance

Positive values indicate lower proximity (neutrophils are farther from the high-state fibroblast population than from matched random fibroblasts). Negative inflamed-minus-noninflamed effects indicate greater neutrophil proximity in inflamed UC:

{chr(10).join(prox_lines)}

Within-group neutrophil-to-fibroblast distance effects were:

{chr(10).join(prox_within_lines)}

Normalized fractions of neutrophils within 50 um of fibroblasts were:

{chr(10).join(contact_lines)}

## Neutrophil phenotype by relative proximity to high fibroblast states

Significant closest-third-versus-farthest-third marker effects within each UC state were:

{chr(10).join(marker_lines)}

Formal inflamed-versus-noninflamed marker-effect comparisons were:

{chr(10).join(between_lines)}

Complete abundance, multiscale 25/50/100 um proximity, raw distance, equal-count proximity, nearest-distance, and marker results are provided in the tables with raw P values and adjusted FDR values. These proximity associations do not establish causal reprogramming or signaling direction.
"""
(OUT / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

manifest = {
    "high_definition": "within-patient fibroblast CLR 75th percentile",
    "included_groups": list(GROUPS),
    "comparison_scope": "UC inflamed versus UC noninflamed only",
    "targets": list(TARGETS), "radii_um": list(RADII),
    "equal_fibroblasts": N_FIBROBLASTS, "equal_neutrophils": N_NEUTROPHILS,
    "equal_proximal_neutrophils": N_PROXIMAL, "equal_distant_neutrophils": N_DISTANT,
    "marker_pool_definition": "closest third versus farthest third within patient",
    "iterations": N_BOOT, "excluded_equal_count_samples": excluded.PatientID.tolist(),
    "FDR_reporting_threshold": FDR_ALPHA,
    "seed": SEED, "biological_replicate": "patient",
}
(OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print("Completed FAP-high/a5B1-high rerun.")
print(f"Saved outputs to {OUT}")
