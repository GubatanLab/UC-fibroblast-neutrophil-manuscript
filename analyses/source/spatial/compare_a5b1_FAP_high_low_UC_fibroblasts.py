from pathlib import Path
import json
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

warnings.filterwarnings("ignore")

ROOT = Path("giotto_codex_results")
OUT = ROOT / "a5b1_FAP_high_low_UC_fibroblasts"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GROUPS = ("UC_Noninflamed", "UC_Inflamed")
GROUP_LABEL = {"UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}
GROUP_COLOR = {"UC_Noninflamed": "#D4AD35", "UC_Inflamed": "#D95F5F"}
A5_POSITIVE_THRESHOLD = 0.9315


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    x = p[ok]
    order = np.argsort(x)
    ranked = x[order] * len(x) / np.arange(1, len(x) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.minimum(ranked, 1)
    out[ok] = adjusted
    return out


def paired_summary(frame, definition):
    rows = []
    for (pid, diagnosis), q in frame.groupby(["PatientID", "Diagnosis2"]):
        low = q[q.FAP_group == "FAP-low"]
        high = q[q.FAP_group == "FAP-high"]
        if low.empty or high.empty:
            continue
        low_clr = low.a5B1_cell_clr.median()
        high_clr = high.a5B1_cell_clr.median()
        low_raw = low.a5B1_cell_raw.median()
        high_raw = high.a5B1_cell_raw.median()
        rows.append(
            {
                "PatientID": pid,
                "Diagnosis2": diagnosis,
                "definition": definition,
                "n_FAP_low": len(low),
                "n_FAP_high": len(high),
                "median_a5B1_clr_FAP_low": low_clr,
                "median_a5B1_clr_FAP_high": high_clr,
                "delta_a5B1_clr_high_minus_low": high_clr - low_clr,
                "median_a5B1_raw_FAP_low": low_raw,
                "median_a5B1_raw_FAP_high": high_raw,
                "log2_raw_a5B1_high_low": np.log2((high_raw + 0.1) / (low_raw + 0.1)),
                "fraction_a5B1_positive_FAP_low": np.mean(low.a5B1_cell_clr >= A5_POSITIVE_THRESHOLD),
                "fraction_a5B1_positive_FAP_high": np.mean(high.a5B1_cell_clr >= A5_POSITIVE_THRESHOLD),
            }
        )
    out = pd.DataFrame(rows)
    out["delta_fraction_a5B1_positive_high_minus_low"] = (
        out.fraction_a5B1_positive_FAP_high - out.fraction_a5B1_positive_FAP_low
    )
    return out


def paired_tests(summary, definition):
    rows = []
    metrics = (
        ("delta_a5B1_clr_high_minus_low", "α5β1 CLR median difference"),
        ("log2_raw_a5B1_high_low", "log2 raw α5β1 median ratio"),
        ("delta_fraction_a5B1_positive_high_minus_low", "α5β1-positive fraction difference"),
    )
    strata = (("All UC", summary),) + tuple(
        (GROUP_LABEL[group], summary[summary.Diagnosis2 == group]) for group in GROUPS
    )
    for label, q in strata:
        for metric, endpoint in metrics:
            values = q[metric].dropna().to_numpy()
            test = wilcoxon(values, alternative="two-sided", zero_method="wilcox")
            rows.append(
                {
                    "definition": definition,
                    "comparison": "FAP-high versus FAP-low",
                    "stratum": label,
                    "endpoint": endpoint,
                    "n_patients": len(values),
                    "median_effect": np.median(values),
                    "mean_effect": np.mean(values),
                    "positive_effect_patients": int(np.sum(values > 0)),
                    "p_value": test.pvalue,
                }
            )
    for metric, endpoint in metrics:
        noninflamed = summary[summary.Diagnosis2 == "UC_Noninflamed"][metric].dropna()
        inflamed = summary[summary.Diagnosis2 == "UC_Inflamed"][metric].dropna()
        test = mannwhitneyu(inflamed, noninflamed, alternative="two-sided")
        rows.append(
            {
                "definition": definition,
                "comparison": "UC inflamed versus UC noninflamed effect",
                "stratum": "Between UC states",
                "endpoint": endpoint,
                "n_patients": len(inflamed) + len(noninflamed),
                "median_effect": np.median(inflamed) - np.median(noninflamed),
                "mean_effect": np.mean(inflamed) - np.mean(noninflamed),
                "positive_effect_patients": np.nan,
                "p_value": test.pvalue,
            }
        )
    out = pd.DataFrame(rows)
    out["fdr_within_definition"] = bh_adjust(out.p_value)
    return out


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


d = pd.read_csv(ROOT / "additional_analyses" / "codex_extended_cells.csv")
required = {
    "cell_ID", "PatientID", "Diagnosis2", "cell_type",
    "a5B1_cell_raw", "a5B1_cell_clr", "FAPa_cell_raw", "FAPa_cell_clr",
}
missing = required.difference(d.columns)
if missing:
    raise RuntimeError(f"Missing exported columns: {sorted(missing)}")

fib = d[
    d.cell_type.eq("Fibroblast") & d.Diagnosis2.isin(GROUPS)
].copy()
thresholds = (
    fib.groupby("PatientID").FAPa_cell_clr
    .quantile([0.25, 0.5, 0.75])
    .unstack()
    .rename(columns={0.25: "FAP_q25", 0.5: "FAP_median", 0.75: "FAP_q75"})
)
fib = fib.merge(thresholds, left_on="PatientID", right_index=True)
thresholds.reset_index().to_csv(TAB / "00_patient_FAP_thresholds.csv", index=False)

# Primary: within-patient median split.
primary_cells = fib.copy()
primary_cells["FAP_group"] = np.where(
    primary_cells.FAPa_cell_clr > primary_cells.FAP_median, "FAP-high", "FAP-low"
)
primary = paired_summary(primary_cells, "within-patient FAP median split")
primary.to_csv(TAB / "01_patient_a5B1_summary_median_split.csv", index=False)
primary_stats = paired_tests(primary, "within-patient FAP median split")
primary_stats.to_csv(TAB / "02_patient_level_statistics_median_split.csv", index=False)

# Sensitivity: retain only bottom and top patient-specific FAP quartiles.
quartile_cells = fib[
    (fib.FAPa_cell_clr <= fib.FAP_q25) | (fib.FAPa_cell_clr >= fib.FAP_q75)
].copy()
quartile_cells["FAP_group"] = np.where(
    quartile_cells.FAPa_cell_clr >= quartile_cells.FAP_q75, "FAP-high", "FAP-low"
)
quartile = paired_summary(quartile_cells, "within-patient bottom versus top FAP quartile")
quartile.to_csv(TAB / "03_patient_a5B1_summary_extreme_quartiles.csv", index=False)
quartile_stats = paired_tests(quartile, "within-patient bottom versus top FAP quartile")
quartile_stats.to_csv(TAB / "04_patient_level_statistics_extreme_quartiles.csv", index=False)

# Patient-level continuous FAP–α5β1 association.
corr_rows = []
for (pid, diagnosis), q in fib.groupby(["PatientID", "Diagnosis2"]):
    rho_raw, p_raw = spearmanr(q.FAPa_cell_raw, q.a5B1_cell_raw)
    rho_clr, p_clr = spearmanr(q.FAPa_cell_clr, q.a5B1_cell_clr)
    corr_rows.append(
        {
            "PatientID": pid,
            "Diagnosis2": diagnosis,
            "n_fibroblasts": len(q),
            "spearman_rho_FAP_a5B1_raw": rho_raw,
            "spearman_rho_FAP_a5B1_clr": rho_clr,
            "cell_level_raw_spearman_p_descriptive_only": p_raw,
            "cell_level_clr_spearman_p_descriptive_only": p_clr,
        }
    )
corr = pd.DataFrame(corr_rows)
corr.to_csv(TAB / "05_patient_FAP_a5B1_correlations.csv", index=False)
corr_stats = []
for label, q in (("All UC", corr),) + tuple(
    (GROUP_LABEL[group], corr[corr.Diagnosis2 == group]) for group in GROUPS
):
    values = q.spearman_rho_FAP_a5B1_raw.to_numpy()
    test = wilcoxon(values, alternative="two-sided", zero_method="wilcox")
    corr_stats.append(
        {
            "stratum": label,
            "n_patients": len(values),
            "median_spearman_rho": np.median(values),
            "positive_rho_patients": int(np.sum(values > 0)),
            "p_value": test.pvalue,
        }
    )
corr_stats = pd.DataFrame(corr_stats)
corr_stats["fdr"] = bh_adjust(corr_stats.p_value)
corr_stats.to_csv(TAB / "06_patient_correlation_statistics.csv", index=False)

sns.set_theme(style="whitegrid", context="talk")
plt.rcParams.update({"font.family": "Arial", "axes.titleweight": "bold", "axes.linewidth": 1.2})

fig, axes = plt.subplots(1, 3, figsize=(20, 6))
long = primary.melt(
    id_vars=["PatientID", "Diagnosis2"],
    value_vars=["median_a5B1_raw_FAP_low", "median_a5B1_raw_FAP_high"],
    var_name="FAP_group", value_name="patient_median_a5B1_raw",
)
long["FAP_group"] = long.FAP_group.map(
    {"median_a5B1_raw_FAP_low": "FAP-low", "median_a5B1_raw_FAP_high": "FAP-high"}
)
long["patient_median_log2_a5B1_raw"] = np.log2(long.patient_median_a5B1_raw + 0.1)
for _, q in long.groupby("PatientID"):
    axes[0].plot(
        [0, 1], q.set_index("FAP_group").loc[
            ["FAP-low", "FAP-high"], "patient_median_log2_a5B1_raw"
        ],
        color=GROUP_COLOR[q.Diagnosis2.iloc[0]], alpha=0.55, marker="o",
    )
axes[0].set_xticks([0, 1], ["FAP-low", "FAP-high"])
axes[0].set_ylabel("Patient median log2(raw α5β1 + 0.1)")
axes[0].set_title("A  Paired fibroblast comparison")

plot_primary = primary.copy()
plot_primary["Group"] = plot_primary.Diagnosis2.map(GROUP_LABEL)
palette = {GROUP_LABEL[k]: GROUP_COLOR[k] for k in GROUPS}
sns.boxplot(
    data=plot_primary, x="Group", y="log2_raw_a5B1_high_low",
    order=[GROUP_LABEL[g] for g in GROUPS], palette=palette, showfliers=False, ax=axes[1],
)
sns.stripplot(
    data=plot_primary, x="Group", y="log2_raw_a5B1_high_low",
    order=[GROUP_LABEL[g] for g in GROUPS], color="black", size=6, ax=axes[1],
)
axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
axes[1].set_xlabel("")
axes[1].set_ylabel("log2 raw α5β1 ratio: FAP-high / FAP-low")
axes[1].set_title("B  Effect by UC state")

plot_corr = corr.copy()
plot_corr["Group"] = plot_corr.Diagnosis2.map(GROUP_LABEL)
sns.boxplot(
    data=plot_corr, x="Group", y="spearman_rho_FAP_a5B1_raw",
    order=[GROUP_LABEL[g] for g in GROUPS], palette=palette, showfliers=False, ax=axes[2],
)
sns.stripplot(
    data=plot_corr, x="Group", y="spearman_rho_FAP_a5B1_raw",
    order=[GROUP_LABEL[g] for g in GROUPS], color="black", size=6, ax=axes[2],
)
axes[2].axhline(0, color="black", linestyle="--", linewidth=1)
axes[2].set_xlabel("")
axes[2].set_ylabel("Within-patient Spearman ρ")
axes[2].set_title("C  Continuous raw-intensity association")
fig.suptitle("α5β1 expression in FAP-high versus FAP-low UC fibroblasts", fontweight="bold")
fig.tight_layout()
save(fig, "Figure_FAP_A5B1_01_primary_comparison")

main = primary_stats[
    (primary_stats.stratum == "All UC")
    & (primary_stats.endpoint == "log2 raw α5β1 median ratio")
].iloc[0]
noninflamed = primary_stats[
    (primary_stats.stratum == "UC noninflamed")
    & (primary_stats.endpoint == "log2 raw α5β1 median ratio")
].iloc[0]
inflamed = primary_stats[
    (primary_stats.stratum == "UC inflamed")
    & (primary_stats.endpoint == "log2 raw α5β1 median ratio")
].iloc[0]
between = primary_stats[
    (primary_stats.stratum == "Between UC states")
    & (primary_stats.endpoint == "log2 raw α5β1 median ratio")
].iloc[0]
quartile_main = quartile_stats[
    (quartile_stats.stratum == "All UC")
    & (quartile_stats.endpoint == "log2 raw α5β1 median ratio")
].iloc[0]
clr_noninflamed = primary_stats[
    (primary_stats.stratum == "UC noninflamed")
    & (primary_stats.endpoint == "α5β1 CLR median difference")
].iloc[0]
clr_inflamed = primary_stats[
    (primary_stats.stratum == "UC inflamed")
    & (primary_stats.endpoint == "α5β1 CLR median difference")
].iloc[0]
clr_between = primary_stats[
    (primary_stats.stratum == "Between UC states")
    & (primary_stats.endpoint == "α5β1 CLR median difference")
].iloc[0]
corr_main = corr_stats[corr_stats.stratum == "All UC"].iloc[0]

manifest = {
    "input": str(ROOT / "additional_analyses" / "codex_extended_cells.csv"),
    "primary_FAP_definition": "within-patient median split",
    "primary_a5B1_endpoint": "within-patient log2 ratio of median raw cell intensity",
    "normalization_sensitivity": "CLR-normalized a5B1 median difference",
    "sensitivity_FAP_definition": "within-patient bottom versus top quartile",
    "biological_replicate": "patient",
    "n_UC_patients": int(primary.PatientID.nunique()),
    "n_UC_fibroblasts": int(len(fib)),
}
(OUT / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

report = f"""# α5β1 expression in FAP-high versus FAP-low UC fibroblasts

## Methods

The analysis included {len(fib):,} fibroblasts from {primary.PatientID.nunique()} UC patients. FAP-high and FAP-low were defined separately within each patient using the median cell-level CLR-normalized FAP intensity, which controls patient-specific staining and batch shifts. The primary α5β1 endpoint was the within-patient log2 ratio of median raw cell intensity. CLR-normalized α5β1 was evaluated as a sensitivity endpoint. Paired Wilcoxon tests used patients—not cells—as replicates.

## Results

Across all UC patients, FAP-high fibroblasts had a median α5β1 raw-intensity log2 ratio of {main.median_effect:.3f}, equivalent to {100 * (2 ** main.median_effect - 1):.1f}% higher median expression than FAP-low fibroblasts ({int(main.positive_effect_patients)}/{int(main.n_patients)} patients positive; paired Wilcoxon P={main.p_value:.4g}, FDR={main.fdr_within_definition:.4g}).

The median log2 ratio was {noninflamed.median_effect:.3f} in UC noninflamed tissue ({100 * (2 ** noninflamed.median_effect - 1):.1f}% higher; {int(noninflamed.positive_effect_patients)}/{int(noninflamed.n_patients)} patients; P={noninflamed.p_value:.4g}, FDR={noninflamed.fdr_within_definition:.4g}) and {inflamed.median_effect:.3f} in UC inflamed tissue ({100 * (2 ** inflamed.median_effect - 1):.1f}% higher; {int(inflamed.positive_effect_patients)}/{int(inflamed.n_patients)} patients; P={inflamed.p_value:.4g}, FDR={inflamed.fdr_within_definition:.4g}). The raw-intensity effect did not differ significantly between UC states (Mann–Whitney P={between.p_value:.4g}, FDR={between.fdr_within_definition:.4g}).

## Robustness

Restricting the analysis to the bottom and top patient-specific FAP quartiles produced a median raw α5β1 log2 ratio of {quartile_main.median_effect:.3f}, equivalent to {100 * (2 ** quartile_main.median_effect - 1):.1f}% higher expression ({int(quartile_main.positive_effect_patients)}/{int(quartile_main.n_patients)} patients positive; paired Wilcoxon P={quartile_main.p_value:.4g}, FDR={quartile_main.fdr_within_definition:.4g}).

As a continuous measure, the median within-patient FAP–α5β1 Spearman correlation was {corr_main.median_spearman_rho:.3f} ({int(corr_main.positive_rho_patients)}/{int(corr_main.n_patients)} patients positive; Wilcoxon P={corr_main.p_value:.4g}, FDR={corr_main.fdr:.4g}).

CLR normalization produced a state-dependent sensitivity result: FAP-high fibroblasts had higher relative α5β1 in UC noninflamed tissue (median difference {clr_noninflamed.median_effect:.3f}; FDR={clr_noninflamed.fdr_within_definition:.4g}) but not in inflamed tissue (median difference {clr_inflamed.median_effect:.3f}; FDR={clr_inflamed.fdr_within_definition:.4g}); the CLR effects differed between states (FDR={clr_between.fdr_within_definition:.4g}). Because CLR expresses α5β1 relative to each cell's full marker profile, this does not negate the modest increase in absolute α5β1 intensity.

## Interpretation

This analysis evaluates protein co-expression within phenotypically annotated fibroblasts. It does not establish that FAP regulates α5β1 or that the two proteins form a functional signaling unit.
"""
(OUT / "RESULTS_SUMMARY.md").write_text(report, encoding="utf-8")

print("Completed FAP-high versus FAP-low fibroblast comparison.")
print(f"Saved outputs to {OUT}")
