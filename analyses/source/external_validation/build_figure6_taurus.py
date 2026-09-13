from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(r"input_data/mouse\Figure 6 TAURUS External Validation")
TABLES = ROOT / "tables"
OUT = ROOT / "final"
OUT.mkdir(parents=True, exist_ok=True)

PURPLE = "#6A51A3"
PURPLE_L = "#D9D2E9"
TEAL = "#147D73"
TEAL_L = "#CBE9E4"
ORANGE = "#D56A1C"
ORANGE_L = "#F5D5BD"
GOLD = "#B07A00"
NAVY = "#17365D"
GREY = "#667085"
LIGHT = "#EEF1F5"
BLACK = "#17202A"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def panel_label(ax, letter):
    ax.text(-0.12, 1.08, letter, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top", ha="left")


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3, width=0.7)


def fmt_p(p):
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return f"{p:.1e}"
    return f"{p:.3f}"


def rank_biserial(x):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x) & (x != 0)]
    if len(x) == 0:
        return np.nan
    r = rankdata(np.abs(x), method="average")
    return (r[x > 0].sum() - r[x < 0].sum()) / r.sum()


def rb_ci(x, seed):
    x = np.asarray(x, dtype=float)
    rng = np.random.default_rng(seed)
    vals = [rank_biserial(x[rng.integers(0, len(x), len(x))]) for _ in range(10000)]
    return np.nanquantile(vals, [0.025, 0.975])


infl = pd.read_csv(TABLES / "taurus_uc_baseline_inflammation_patient_deltas.tsv", sep="\t")
infl_stat = pd.read_csv(TABLES / "taurus_uc_baseline_inflammation_statistics.tsv", sep="\t")
state = pd.read_csv(TABLES / "taurus_uc_alpha5beta1_state_contrast_by_patient.tsv", sep="\t")
state_stat = pd.read_csv(TABLES / "taurus_uc_alpha5beta1_state_contrast_statistics.tsv", sep="\t").iloc[0]
paired = pd.read_csv(TABLES / "taurus_uc_patient_paired_deltas.tsv", sep="\t")
long_stat = pd.read_csv(TABLES / "taurus_uc_longitudinal_statistics.tsv", sep="\t")
corr = pd.read_csv(TABLES / "taurus_uc_baseline_cross_compartment.tsv", sep="\t")
corr_stat = pd.read_csv(TABLES / "taurus_uc_cross_compartment_statistics.tsv", sep="\t").iloc[0]

fig = plt.figure(figsize=(13.2, 8.7), facecolor="white")
outer = gridspec.GridSpec(2, 3, figure=fig, width_ratios=[0.95, 1.12, 1.08], height_ratios=[0.92, 1.08], hspace=0.48, wspace=0.42)
axA = fig.add_subplot(outer[0, 0])
axB = fig.add_subplot(outer[0, 1])
axC = fig.add_subplot(outer[0, 2])
axD = fig.add_subplot(outer[1, 0])
axE = fig.add_subplot(outer[1, 1])
axF = fig.add_subplot(outer[1, 2])

fig.suptitle("Independent UC validation links FAP⁺ α5β1-high fibroblasts to neutrophil-attracting inflammation and anti-TNF nonremission", x=0.5, y=0.982, fontsize=15, fontweight="bold", color=BLACK)
fig.text(0.5, 0.952, "TAURUS single-cell atlas  |  fixed manuscript-defined signatures  |  patient-level inference", ha="center", va="center", fontsize=9.5, color=GREY)

# A — design
axA.set_axis_off()
panel_label(axA, "A")
axA.set_title("External validation design", loc="left", fontweight="bold", pad=8)
boxes = [
    (0.06, 0.72, 0.88, 0.18, TEAL_L, "TAURUS UC cohort\n22 patients • 108 biopsies"),
    (0.06, 0.43, 0.88, 0.20, PURPLE_L, "Locked stromal tests\nFAP-inflammatory • α5β1 • recruitment"),
    (0.06, 0.12, 0.88, 0.22, ORANGE_L, "Patient-level contrasts\nsite-paired inflammation • pre→post outcome"),
]
for x, y, w, h, fc, txt in boxes:
    axA.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.015,rounding_size=0.025", facecolor=fc, edgecolor="white", linewidth=1.2))
    axA.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=8.4, color=BLACK, linespacing=1.35)
for y1, y2 in [(0.71, 0.64), (0.42, 0.35)]:
    axA.add_patch(FancyArrowPatch((0.5, y1), (0.5, y2), arrowstyle="-|>", mutation_scale=10, linewidth=1.1, color=GREY))
axA.text(0.06, 0.01, "Neutrophils were not captured in frozen-tissue 10x;\nS100A8/A9-high monocytes provide a labeled feedback proxy.", fontsize=7.1, color=GREY, va="bottom")

# B — baseline inflammation replication forest
panel_label(axB, "B")
axB.set_title("Pretreatment inflamed sites recapitulate\nthe stromal program", loc="left", fontweight="bold", fontsize=9.2, pad=8)
order = ["activated FAP fibroblast fraction", "FAP inflammatory", "alpha5beta1 adhesion", "neutrophil recruitment", "OSM response", "myeloid feedback proxy"]
labels = ["Activated FAP⁺ fraction", "FAP-inflammatory", "α5β1 adhesion", "Neutrophil recruitment", "OSM response", "Myeloid feedback proxy"]
effects = []
for i, feat in enumerate(order):
    x = infl.loc[infl["feature"].eq(feat), "delta"].to_numpy(float)
    stat = infl_stat.loc[infl_stat["feature"].eq(feat)].iloc[0]
    lo, hi = rb_ci(x, 100 + i)
    effects.append((rank_biserial(x), lo, hi, int(stat["n_patients"]), stat["FDR"]))
y = np.arange(len(order))[::-1]
for yy, feat, (eff, lo, hi, n, q) in zip(y, order, effects):
    color = GOLD if feat == "myeloid feedback proxy" else PURPLE
    axB.plot([lo, hi], [yy, yy], color=color, lw=2)
    axB.scatter(eff, yy, s=42, color=color, edgecolor="white", linewidth=0.6, zorder=3)
    axB.text(1.04, yy, f"n={n}  q={fmt_p(q)}", va="center", ha="left", fontsize=6.9, color=GREY, clip_on=False)
axB.axvline(0, color="#98A2B3", lw=0.9, ls="--")
axB.set_xlim(-1.05, 1.05)
axB.set_yticks(y, labels)
axB.set_xlabel("Matched-pairs rank-biserial effect\n(inflamed − noninflamed)")
clean(axB)

# C — alpha5beta1 enrichment in activated state
panel_label(axC, "C")
axC.set_title("α5β1 activity localizes to activated\nFAP⁺ fibroblasts", loc="left", fontweight="bold", fontsize=9.2, pad=8)
for _, r in state.iterrows():
    axC.plot([0, 1], [r["other"], r["activated"]], color="#B7BDC7", lw=0.8, alpha=0.8, zorder=1)
axC.scatter(np.zeros(len(state)), state["other"], s=22, color="#A7A9AC", alpha=0.9, edgecolor="white", linewidth=0.4, zorder=2)
axC.scatter(np.ones(len(state)), state["activated"], s=28, color=PURPLE, alpha=0.95, edgecolor="white", linewidth=0.4, zorder=2)
med = [state["other"].median(), state["activated"].median()]
axC.plot([0, 1], med, color=BLACK, lw=2.5, zorder=3)
axC.set_xlim(-0.35, 1.35)
axC.set_xticks([0, 1], ["Other\nfibroblasts", "THY1⁺FAP⁺PDPN⁺\nfibroblasts"])
axC.set_ylabel("α5β1 adhesion module score")
axC.text(0.5, 0.97, f"n={int(state_stat['n_patients'])} patients  •  P={fmt_p(state_stat['wilcoxon_greater_p'])}", transform=axC.transAxes, ha="center", va="top", fontsize=7.6, color=PURPLE)
clean(axC)


def paired_panel(ax, feature, title, ylabel, percent=False):
    z = paired.loc[paired["feature"].eq(feature)].copy()
    colors = {"Remission": TEAL, "Non_Remission": ORANGE}
    xmap = {"Remission": (0, 1), "Non_Remission": (2.1, 3.1)}
    for outcome in ["Remission", "Non_Remission"]:
        q = z.loc[z["Remission_status"].eq(outcome)]
        xx = xmap[outcome]
        for _, r in q.iterrows():
            vals = np.array([r["pre"], r["post"]]) * (100 if percent else 1)
            ax.plot(xx, vals, color=colors[outcome], alpha=0.36, lw=1)
            ax.scatter(xx, vals, color=colors[outcome], s=18, alpha=0.7, edgecolor="white", linewidth=0.35, zorder=2)
        meds = [q["pre"].median(), q["post"].median()]
        if percent:
            meds = np.array(meds) * 100
        ax.plot(xx, meds, color=colors[outcome], lw=3.2, marker="o", markersize=4.5, zorder=4)
    s_r = long_stat.loc[(long_stat["feature"].eq(feature)) & (long_stat["comparison"].eq("Remission"))].iloc[0]
    s_n = long_stat.loc[(long_stat["feature"].eq(feature)) & (long_stat["comparison"].eq("Non_Remission"))].iloc[0]
    s_i = long_stat.loc[(long_stat["feature"].eq(feature)) & (long_stat["comparison"].str.startswith("delta:"))].iloc[0]
    ax.set_xticks([0, 1, 2.1, 3.1], ["Pre", "Post", "Pre", "Post"])
    ax.text(0.5, -0.17, f"Remission\nn={int(s_r['n_patients'])}", transform=ax.get_xaxis_transform(), ha="center", color=TEAL, fontsize=7.3)
    ax.text(2.6, -0.17, f"Nonremission\nn={int(s_n['n_patients'])}", transform=ax.get_xaxis_transform(), ha="center", color=ORANGE, fontsize=7.3)
    ax.set_title(title, loc="left", fontweight="bold", pad=8)
    ax.set_ylabel(ylabel)
    ax.text(0.02, 0.98, f"within R P={fmt_p(s_r['p_value'])}  |  within NR P={fmt_p(s_n['p_value'])}\nΔ outcome P={fmt_p(s_i['p_value'])}; q={fmt_p(s_i['FDR'])}", transform=ax.transAxes, ha="left", va="top", fontsize=6.8, color=GREY)
    clean(ax)


# D — outcome divergence in abundance
panel_label(axD, "D")
paired_panel(axD, "activated FAP fibroblast fraction", "Activated FAP⁺ fibroblasts diverge\nby treatment outcome", "Activated state (% of stroma)", percent=True)

# E — outcome divergence in recruitment program
panel_label(axE, "E")
paired_panel(axE, "neutrophil recruitment", "The recruitment program falls\ndirectionally in remission", "Fibroblast recruitment module score")

# F — cross-compartment concordance
panel_label(axF, "F")
axF.set_title("Stromal recruitment tracks an\ninflammatory feedback compartment", loc="left", fontweight="bold", fontsize=9.2, pad=8)
for outcome, color, marker in [("Remission", TEAL, "o"), ("Non_Remission", ORANGE, "s"), ("Not_avail", GREY, "^")]:
    z = corr.loc[corr["Remission_status"].eq(outcome)]
    if len(z):
        axF.scatter(z["fib_recruitment"], z["myeloid_feedback"], s=40, color=color, marker=marker, edgecolor="white", linewidth=0.6, alpha=0.9, label=outcome.replace("_", " "))
if len(corr) >= 3:
    coef = np.polyfit(corr["fib_recruitment"], corr["myeloid_feedback"], 1)
    xx = np.linspace(corr["fib_recruitment"].min(), corr["fib_recruitment"].max(), 100)
    axF.plot(xx, np.polyval(coef, xx), color=NAVY, lw=1.5, alpha=0.8)
axF.set_xlabel("Fibroblast neutrophil-recruitment score")
axF.set_ylabel("S100A8/A9-high monocyte\nfeedback-proxy score")
axF.text(0.03, 0.97, f"n={int(corr_stat['n_patients'])} patients\nSpearman ρ={corr_stat['spearman_rho']:.2f}, P={fmt_p(corr_stat['p_value'])}\npartial ρ={corr_stat['partial_spearman_rho']:.2f}, P={fmt_p(corr_stat['partial_p_value'])}", transform=axF.transAxes, ha="left", va="top", fontsize=7.2, color=NAVY, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=LIGHT, alpha=0.94))
axF.legend(frameon=False, fontsize=6.8, loc="lower right", handletextpad=0.4)
clean(axF)

fig.text(0.5, 0.015, "Conclusion: an independent longitudinal UC cohort validates an inflammation-linked, α5β1-high activated-FAP fibroblast state that persists with anti-TNF nonremission; direct neutrophil-state replication awaits a neutrophil-preserving cohort.", ha="center", va="bottom", fontsize=8.3, color=BLACK, bbox=dict(boxstyle="round,pad=0.4", fc="#F7F5FB", ec=PURPLE_L))
fig.subplots_adjust(top=0.91, bottom=0.105, left=0.075, right=0.955)

pdf = OUT / "Figure_6_TAURUS_external_validation.pdf"
png = OUT / "Figure_6_TAURUS_external_validation.png"
fig.savefig(pdf, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(png, dpi=350, bbox_inches="tight", facecolor="white")
plt.close(fig)

legend = """Figure 6. Independent validation of the activated FAP–α5β1 stromal program in the longitudinal TAURUS ulcerative colitis cohort. (A) Validation design. Public TAURUS lineage-specific objects were restricted to ulcerative colitis (22 patients, 108 biopsies). The manuscript-defined FAP-inflammatory, α5β1-adhesion, neutrophil-recruitment and OSM-response gene sets were fixed before outcome testing. Counts were normalized within cell as log1p(counts per 10,000). Biopsies with at least 20 fibroblasts (or 10 S100A8/A9-high monocytes) were retained. Multiple biopsy sites were averaged within patient before inference. (B) Pretreatment, within-patient comparison of inflamed and noninflamed sites. Points show matched-pairs rank-biserial effects; horizontal lines are 95% bootstrap confidence intervals. Positive values indicate higher activity in inflamed tissue. P values were from prespecified one-sided Wilcoxon signed-rank tests and were adjusted across the six displayed features by the Benjamini–Hochberg method; q values and patient counts are shown. The myeloid-feedback proxy comprised NAMPT, OSM, IL1B, TNF and S100A9 in S100A8/A9-high monocytes. (C) Patient-level, within-biopsy enrichment of the fixed α5β1-adhesion module in the authors’ THY1+FAP+PDPN+ activated fibroblast state versus all other nonpericyte fibroblasts. Biopsy contrasts were averaged within patient; the P value is from a one-sided Wilcoxon signed-rank test. (D,E) Site-matched longitudinal changes after adalimumab in activated-FAP fibroblast abundance (D) and the fibroblast neutrophil-recruitment module (E), stratified by clinical remission. Thin lines denote patients and thick lines denote group medians; biopsy-site changes were averaged within patient. Within-group P values are two-sided Wilcoxon signed-rank tests. Outcome-difference P values compare patient-level changes with a two-sided Mann–Whitney test; q values are Benjamini–Hochberg adjusted across the six prespecified longitudinal features. (F) Baseline patient-level association between fibroblast neutrophil-recruitment activity and the inflammatory myeloid-feedback proxy, using only biopsy identifiers represented in both compartments. Spearman correlation and partial rank correlation controlling for the TAURUS inflammation score and remission outcome are shown. TAURUS frozen-tissue 10x data did not capture neutrophils; therefore, panels B and F support the fibroblast recruitment and inflammatory-feedback arms but do not constitute direct replication of neutrophil transcriptional states. The original TAURUS spatial imaging localized S100A9+MPO+CD66B+ neutrophil aggregates near THY1+FAP+PDPN+ fibroblast-associated epithelial damage regions.

Source: Thomas et al., Nature Immunology (2024), https://doi.org/10.1038/s41590-024-01994-8. Data: Zenodo, https://doi.org/10.5281/zenodo.13768607.
"""
(OUT / "Figure_6_TAURUS_external_validation_legend.txt").write_text(legend, encoding="utf-8")

methods = """TAURUS UC external-validation analysis note

The analysis used the public fibperi_final.h5ad and myeloid_final.h5ad lineage objects from the TAURUS Zenodo record. The UC subset contained 22 patients and 108 unique biopsies; outcome metadata identified 6 remission and 14 nonremission patients, with outcome unavailable for 2 patients. Raw counts were normalized per cell as log1p(count / total_counts × 10,000). No gene selection or model fitting used remission labels.

Primary fixed fibroblast modules were copied from the manuscript analysis: FAP inflammatory (FAP, PDPN, CXCL1, CXCL5, CXCL6, CXCL8, CSF3, IL6, ICAM1); alpha5beta1 adhesion (ITGA5, ITGB1, FN1, PTK2, SRC, PXN, VCL, TLN1, ACTN1, RHOA, ROCK1, MYL9); neutrophil recruitment (CXCL1, CXCL2, CXCL3, CXCL5, CXCL6, CXCL8, CSF3); and OSM response (OSMR, LIFR, IL6ST, STAT3, SOCS3, CEBPD, JUNB, FOSL2, CXCL1, CXCL8). Module scores were mean normalized expression across member genes and then averaged across fibroblasts within each biopsy. The activated-FAP abundance was the fraction of THY1+FAP+PDPN+ fibroblasts among all fibroblast/pericyte cells, following the TAURUS annotation and denominator.

The TAURUS authors explicitly report that neutrophils were not captured. The secondary myeloid-feedback proxy (NAMPT, OSM, IL1B, TNF, S100A9) was therefore scored only in the annotated S100A8/A9-high monocyte states and is not represented as direct neutrophil validation.

Inference used patients, not cells, as independent units. Pretreatment inflamed and noninflamed biopsy scores were averaged within patient and compared with a one-sided Wilcoxon signed-rank test in the prespecified direction. Site-matched pre/post changes were first calculated within patient and site, then averaged within patient; within-outcome tests used two-sided Wilcoxon signed-rank tests and change-by-outcome comparisons used two-sided Mann–Whitney tests. Benjamini–Hochberg adjustment was applied within the six-feature baseline and longitudinal families. Correlation analyses used baseline patient means from biopsy identifiers observed in both lineages; the partial rank correlation controlled for the supplied inflammation score and remission outcome.
"""
(OUT / "TAURUS_external_validation_methods.txt").write_text(methods, encoding="utf-8")

print(pdf)
print(png)
