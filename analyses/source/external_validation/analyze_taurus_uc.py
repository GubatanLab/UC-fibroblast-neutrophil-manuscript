from pathlib import Path
import json
import warnings

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

ROOT = Path(r"input_data/mouse\Figure 6 TAURUS External Validation")
DATA = ROOT / "data"
TABLES = ROOT / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

FIB_MODULES = {
    "FAP inflammatory": ["FAP", "PDPN", "CXCL1", "CXCL5", "CXCL6", "CXCL8", "CSF3", "IL6", "ICAM1"],
    "alpha5beta1 adhesion": ["ITGA5", "ITGB1", "FN1", "PTK2", "SRC", "PXN", "VCL", "TLN1", "ACTN1", "RHOA", "ROCK1", "MYL9"],
    "neutrophil recruitment": ["CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CSF3"],
    "OSM response": ["OSMR", "LIFR", "IL6ST", "STAT3", "SOCS3", "CEBPD", "JUNB", "FOSL2", "CXCL1", "CXCL8"],
}
MYELOID_MODULES = {
    "myeloid feedback proxy": ["NAMPT", "OSM", "IL1B", "TNF", "S100A9"],
}
MIN_FIB = 20
MIN_MY = 10
MIN_ACTIVATED = 5
RNG = np.random.default_rng(20260816)


def extract_normalized(path, genes):
    """Extract only requested genes from backed raw counts, then log-normalize per cell."""
    obj = ad.read_h5ad(path, backed="r")
    genes = [g for g in genes if g in obj.var_names]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sub = obj[:, genes].to_memory()
    obs = obj.obs.copy()
    obj.file.close()
    x = sub.X.astype(np.float64)
    totals = obs["total_counts"].to_numpy(dtype=float)
    scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
    if sparse.issparse(x):
        x = sparse.diags(scale) @ x
        x.data = np.log1p(x.data)
        expr = pd.DataFrame.sparse.from_spmatrix(x, index=obs.index, columns=genes)
    else:
        expr = pd.DataFrame(np.log1p(x * scale[:, None]), index=obs.index, columns=genes)
    return obs, expr


def score_modules(obs, expr, modules, cell_mask, min_cells, prefix):
    use_obs = obs.loc[cell_mask].copy()
    use_expr = expr.loc[cell_mask]
    keys = ["sample_id", "Patient", "Site", "Treatment", "Inflammation", "Inflammation_score", "Remission_status"]
    base = use_obs[keys].drop_duplicates("sample_id").set_index("sample_id")
    counts = use_obs.groupby("sample_id", observed=True).size().rename("n_cells")
    out = base.join(counts)
    for label, genes in modules.items():
        present = [g for g in genes if g in use_expr.columns]
        cell_score = use_expr[present].mean(axis=1).astype(float)
        sample_score = cell_score.groupby(use_obs["sample_id"], observed=True).mean()
        out[label] = sample_score
        out.loc[out["n_cells"] < min_cells, label] = np.nan
    out["compartment"] = prefix
    return out.reset_index()


def paired_patient_deltas(sample_features, feature):
    z = sample_features.loc[
        sample_features["Remission_status"].isin(["Remission", "Non_Remission"]),
        ["Patient", "Site", "Treatment", "Remission_status", feature],
    ].dropna()
    wide = z.pivot_table(index=["Patient", "Site", "Remission_status"], columns="Treatment", values=feature, aggfunc="mean").reset_index()
    wide = wide.dropna(subset=["Pre", "Post"])
    wide["delta"] = wide["Post"] - wide["Pre"]
    patient = wide.groupby(["Patient", "Remission_status"], observed=True).agg(
        pre=("Pre", "mean"), post=("Post", "mean"), delta=("delta", "mean"), n_sites=("Site", "nunique")
    ).reset_index()
    patient["feature"] = feature
    return wide, patient


def bootstrap_median_ci(x, n_boot=10000):
    x = np.asarray(x, dtype=float)
    if len(x) == 0:
        return np.nan, np.nan
    idx = RNG.integers(0, len(x), size=(n_boot, len(x)))
    vals = np.median(x[idx], axis=1)
    return tuple(np.quantile(vals, [0.025, 0.975]))


def summarize_deltas(patient):
    rows = []
    for feature, dat in patient.groupby("feature", observed=True):
        for outcome, z in dat.groupby("Remission_status", observed=True):
            x = z["delta"].to_numpy(dtype=float)
            p = wilcoxon(x, alternative="two-sided", zero_method="wilcox").pvalue if len(x) >= 3 and np.any(x != 0) else np.nan
            lo, hi = bootstrap_median_ci(x)
            rows.append({"feature": feature, "comparison": outcome, "n_patients": len(x), "median_delta": np.median(x), "ci_low": lo, "ci_high": hi, "p_value": p})
        r = dat.loc[dat["Remission_status"].eq("Remission"), "delta"].to_numpy(dtype=float)
        nr = dat.loc[dat["Remission_status"].eq("Non_Remission"), "delta"].to_numpy(dtype=float)
        p = mannwhitneyu(r, nr, alternative="two-sided").pvalue if len(r) and len(nr) else np.nan
        rows.append({"feature": feature, "comparison": "delta: remission vs nonremission", "n_patients": len(r) + len(nr), "median_delta": np.median(r) - np.median(nr), "ci_low": np.nan, "ci_high": np.nan, "p_value": p})
    ans = pd.DataFrame(rows)
    ans["FDR"] = np.nan
    for comp in ans["comparison"].unique():
        ix = ans["comparison"].eq(comp) & ans["p_value"].notna()
        if ix.any():
            ans.loc[ix, "FDR"] = multipletests(ans.loc[ix, "p_value"], method="fdr_bh")[1]
    return ans


def baseline_patient(sample_features, feature):
    z = sample_features.loc[
        sample_features["Treatment"].eq("Pre") & sample_features["Remission_status"].isin(["Remission", "Non_Remission"]),
        ["Patient", "Remission_status", "Inflammation", "Inflammation_score", feature],
    ].dropna(subset=[feature])
    return z.groupby(["Patient", "Remission_status"], observed=True).agg(
        score=(feature, "mean"), inflammation_score=("Inflammation_score", "mean"), n_biopsies=(feature, "size")
    ).reset_index()


def inflammation_patient_deltas(sample_features, feature):
    z = sample_features.loc[
        sample_features["Treatment"].eq("Pre") & sample_features["Inflammation"].isin(["Inflamed", "Non_Inflamed"]),
        ["Patient", "Inflammation", feature],
    ].dropna()
    patient = z.groupby(["Patient", "Inflammation"], observed=True)[feature].mean().unstack()
    patient = patient.dropna(subset=["Inflamed", "Non_Inflamed"]).reset_index()
    patient["delta"] = patient["Inflamed"] - patient["Non_Inflamed"]
    patient["feature"] = feature
    return patient


def partial_rank_corr(x, y, covariates):
    dat = pd.DataFrame({"x": x, "y": y}).join(pd.DataFrame(covariates)).dropna()
    if len(dat) < 6:
        return np.nan, np.nan, len(dat)
    ranks = dat.rank(method="average")
    design = np.column_stack([np.ones(len(ranks)), ranks.drop(columns=["x", "y"]).to_numpy()])
    rx = ranks["x"].to_numpy() - design @ np.linalg.lstsq(design, ranks["x"].to_numpy(), rcond=None)[0]
    ry = ranks["y"].to_numpy() - design @ np.linalg.lstsq(design, ranks["y"].to_numpy(), rcond=None)[0]
    rho, p = spearmanr(rx, ry)
    return rho, p, len(dat)


fib_genes = sorted(set(sum(FIB_MODULES.values(), [])))
my_genes = sorted(set(sum(MYELOID_MODULES.values(), [])))
fib_obs, fib_expr = extract_normalized(DATA / "fibperi_final.h5ad", fib_genes)
my_obs, my_expr = extract_normalized(DATA / "myeloid_final.h5ad", my_genes)

uc_fib = fib_obs["Disease"].eq("UC")
fib_lineage = uc_fib & ~fib_obs["minor"].eq("Pericyte")
fib_sample = score_modules(fib_obs, fib_expr, FIB_MODULES, fib_lineage, MIN_FIB, "fibroblast")

# Within-biopsy state contrast: does the independently annotated activated-FAP state carry alpha5beta1?
activated_cells = uc_fib & fib_obs["final_analysis"].eq("THY1pos FAPpos PDPNpos fibroblast")
other_fib_cells = fib_lineage & ~fib_obs["final_analysis"].eq("THY1pos FAPpos PDPNpos fibroblast")
activated_state_scores = score_modules(fib_obs, fib_expr, {"alpha5beta1 adhesion": FIB_MODULES["alpha5beta1 adhesion"]}, activated_cells, MIN_ACTIVATED, "activated FAP fibroblast")
other_state_scores = score_modules(fib_obs, fib_expr, {"alpha5beta1 adhesion": FIB_MODULES["alpha5beta1 adhesion"]}, other_fib_cells, MIN_FIB, "other fibroblasts")
state_contrast = activated_state_scores[["sample_id", "Patient", "alpha5beta1 adhesion"]].merge(
    other_state_scores[["sample_id", "alpha5beta1 adhesion"]], on="sample_id", suffixes=("_activated", "_other")
).dropna()
state_contrast["difference"] = state_contrast["alpha5beta1 adhesion_activated"] - state_contrast["alpha5beta1 adhesion_other"]
state_patient = state_contrast.groupby("Patient", observed=True).agg(
    activated=("alpha5beta1 adhesion_activated", "mean"), other=("alpha5beta1 adhesion_other", "mean"), difference=("difference", "mean"), n_biopsies=("sample_id", "nunique")
).reset_index()
state_p = wilcoxon(state_patient["difference"], alternative="greater").pvalue if len(state_patient) >= 3 else np.nan

# Activated-fibroblast fraction follows the authors' denominator: all fibroblasts/pericytes.
fib_uc_obs = fib_obs.loc[uc_fib]
denom = fib_uc_obs.groupby("sample_id", observed=True).size().rename("n_stromal")
num = fib_uc_obs["final_analysis"].eq("THY1pos FAPpos PDPNpos fibroblast").groupby(fib_uc_obs["sample_id"], observed=True).sum().rename("n_activated_fap")
abund_meta = fib_uc_obs[["sample_id", "Patient", "Site", "Treatment", "Inflammation", "Inflammation_score", "Remission_status"]].drop_duplicates("sample_id").set_index("sample_id")
activated = abund_meta.join([denom, num]).fillna({"n_activated_fap": 0})
activated["activated FAP fibroblast fraction"] = activated["n_activated_fap"] / activated["n_stromal"]
activated.loc[activated["n_stromal"] < MIN_FIB, "activated FAP fibroblast fraction"] = np.nan
activated = activated.reset_index()

damage_mono = my_obs["final_analysis"].isin(["S100A8 A9hi mono", "S100A8 A9hi TNFhi IL6pos mono"])
my_mask = my_obs["Disease"].eq("UC") & damage_mono
my_sample = score_modules(my_obs, my_expr, MYELOID_MODULES, my_mask, MIN_MY, "S100A8/A9-high monocytes")

fib_sample.to_csv(TABLES / "taurus_uc_fibroblast_sample_scores.tsv", sep="\t", index=False)
activated.to_csv(TABLES / "taurus_uc_activated_fap_abundance.tsv", sep="\t", index=False)
my_sample.to_csv(TABLES / "taurus_uc_myeloid_sample_scores.tsv", sep="\t", index=False)
state_contrast.to_csv(TABLES / "taurus_uc_alpha5beta1_state_contrast_by_biopsy.tsv", sep="\t", index=False)
state_patient.to_csv(TABLES / "taurus_uc_alpha5beta1_state_contrast_by_patient.tsv", sep="\t", index=False)
pd.DataFrame([{"n_patients": len(state_patient), "median_activated_minus_other": state_patient["difference"].median(), "wilcoxon_greater_p": state_p}]).to_csv(
    TABLES / "taurus_uc_alpha5beta1_state_contrast_statistics.tsv", sep="\t", index=False
)

feature_frames = {
    "activated FAP fibroblast fraction": activated,
    **{k: fib_sample for k in FIB_MODULES},
    "myeloid feedback proxy": my_sample,
}
site_deltas = []
patient_deltas = []
baseline = []
inflammation_deltas = []
for feature, frame in feature_frames.items():
    site, patient = paired_patient_deltas(frame, feature)
    site["feature"] = feature
    site_deltas.append(site)
    patient_deltas.append(patient)
    b = baseline_patient(frame, feature)
    b["feature"] = feature
    baseline.append(b)
    inflammation_deltas.append(inflammation_patient_deltas(frame, feature))
site_deltas = pd.concat(site_deltas, ignore_index=True)
patient_deltas = pd.concat(patient_deltas, ignore_index=True)
baseline = pd.concat(baseline, ignore_index=True)
inflammation_deltas = pd.concat(inflammation_deltas, ignore_index=True)
summary = summarize_deltas(patient_deltas)
site_deltas.to_csv(TABLES / "taurus_uc_site_matched_deltas.tsv", sep="\t", index=False)
patient_deltas.to_csv(TABLES / "taurus_uc_patient_paired_deltas.tsv", sep="\t", index=False)
baseline.to_csv(TABLES / "taurus_uc_patient_baseline_scores.tsv", sep="\t", index=False)
summary.to_csv(TABLES / "taurus_uc_longitudinal_statistics.tsv", sep="\t", index=False)

infl_rows = []
for feature, z in inflammation_deltas.groupby("feature", observed=True):
    x = z["delta"].to_numpy(dtype=float)
    p = wilcoxon(x, alternative="greater").pvalue if len(x) >= 3 and np.any(x != 0) else np.nan
    lo, hi = bootstrap_median_ci(x)
    infl_rows.append({"feature": feature, "n_patients": len(x), "median_inflamed_minus_noninflamed": np.median(x), "ci_low": lo, "ci_high": hi, "p_value": p})
infl_summary = pd.DataFrame(infl_rows)
ix = infl_summary["p_value"].notna()
infl_summary.loc[ix, "FDR"] = multipletests(infl_summary.loc[ix, "p_value"], method="fdr_bh")[1]
inflammation_deltas.to_csv(TABLES / "taurus_uc_baseline_inflammation_patient_deltas.tsv", sep="\t", index=False)
infl_summary.to_csv(TABLES / "taurus_uc_baseline_inflammation_statistics.tsv", sep="\t", index=False)

# Baseline patient-level cross-compartment concordance using only biopsy IDs represented in both lineages.
fib_b = fib_sample.loc[fib_sample["Treatment"].eq("Pre"), ["sample_id", "Patient", "Remission_status", "Inflammation_score", "neutrophil recruitment"]]
my_b = my_sample.loc[my_sample["Treatment"].eq("Pre"), ["sample_id", "myeloid feedback proxy"]]
corr_sample = fib_b.merge(my_b, on="sample_id").dropna()
corr = corr_sample.groupby(["Patient", "Remission_status"], observed=True).agg(
    fib_recruitment=("neutrophil recruitment", "mean"), myeloid_feedback=("myeloid feedback proxy", "mean"),
    inflammation_score=("Inflammation_score", "mean"), n_biopsies=("sample_id", "nunique")
).reset_index()
rho, p_corr = spearmanr(corr["fib_recruitment"], corr["myeloid_feedback"]) if len(corr) >= 4 else (np.nan, np.nan)
corr["nonremission"] = corr["Remission_status"].eq("Non_Remission").astype(int)
partial_rho, partial_p, partial_n = partial_rank_corr(
    corr["fib_recruitment"], corr["myeloid_feedback"], corr[["inflammation_score", "nonremission"]]
)
corr_sample.to_csv(TABLES / "taurus_uc_baseline_cross_compartment_by_biopsy.tsv", sep="\t", index=False)
corr.to_csv(TABLES / "taurus_uc_baseline_cross_compartment.tsv", sep="\t", index=False)
pd.DataFrame([{"n_patients": len(corr), "spearman_rho": rho, "p_value": p_corr, "partial_spearman_rho": partial_rho, "partial_p_value": partial_p, "partial_n": partial_n}]).to_csv(TABLES / "taurus_uc_cross_compartment_statistics.tsv", sep="\t", index=False)

manifest = {
    "cohort": "TAURUS ulcerative colitis only",
    "n_uc_patients": int(fib_obs.loc[uc_fib, "Patient"].nunique()),
    "n_uc_biopsies": int(fib_obs.loc[uc_fib, "sample_id"].nunique()),
    "n_remission_patients": int(fib_obs.loc[uc_fib & fib_obs["Remission_status"].eq("Remission"), "Patient"].nunique()),
    "n_nonremission_patients": int(fib_obs.loc[uc_fib & fib_obs["Remission_status"].eq("Non_Remission"), "Patient"].nunique()),
    "normalization": "log1p(raw count / cell total counts * 10,000)",
    "fibroblast_min_cells_per_biopsy": MIN_FIB,
    "myeloid_min_cells_per_biopsy": MIN_MY,
    "activated_state_min_cells_per_biopsy": MIN_ACTIVATED,
    "neutrophils_present": False,
    "interpretation_limit": "Direct neutrophil-state replication is impossible because TAURUS did not capture neutrophils; myeloid feedback is a prespecified proxy in S100A8/A9-high monocytes.",
    "fibroblast_modules": FIB_MODULES,
    "myeloid_modules": MYELOID_MODULES,
}
(TABLES / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

print(json.dumps(manifest, indent=2))
print("\nLongitudinal statistics:\n", summary.to_string(index=False))
print("\nBaseline inflamed vs noninflamed:\n", infl_summary.to_string(index=False))
print("\nAlpha5beta1 state contrast:\n", {"n": len(state_patient), "median": state_patient["difference"].median(), "p": state_p})
print("\nCross-compartment:\n", {"n": len(corr), "rho": rho, "p": p_corr, "partial_rho": partial_rho, "partial_p": partial_p})
