from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.stats import mannwhitneyu, kruskal
import seaborn as sns


ROOT = Path(r"input_data/codex\giotto_codex_results")
CELL_FILE = ROOT / "additional_analyses" / "codex_extended_cells.csv"
UMAP_FILE = ROOT / "multipanel" / "umap_cells.csv"
OUT = Path(r"input_data/mouse\CODEX_Figure_2")
OUT.mkdir(parents=True, exist_ok=True)

GROUPS = ["Control", "UC_Noninflamed", "UC_Inflamed"]
GROUP_LABEL = {"Control": "Control", "UC_Noninflamed": "UC noninflamed", "UC_Inflamed": "UC inflamed"}
GROUP_COLOR = {"Control": "#4C78A8", "UC_Noninflamed": "#C9A227", "UC_Inflamed": "#D95F5F"}

BROAD_ORDER = [
    "Epithelial Cell", "Enteroendocrine Cell", "Endothelial Cell", "Fibroblast",
    "Macrophage", "Dendritic Cell", "Neutrophil", "CD4 T", "CD8 T", "TReg",
    "B Cell", "Plasma B Cell"
]
BROAD_LABEL = {
    "Epithelial Cell": "Epithelial", "Enteroendocrine Cell": "Enteroendocrine",
    "Endothelial Cell": "Endothelial", "Fibroblast": "Fibroblast",
    "Macrophage": "Macrophage", "Dendritic Cell": "Dendritic",
    "Neutrophil": "Neutrophil", "CD4 T": "CD4 T", "CD8 T": "CD8 T",
    "TReg": "Treg", "B Cell": "B cell", "Plasma B Cell": "Plasma B"
}
BROAD_COLOR = dict(zip(BROAD_ORDER, [
    "#4E79A7", "#A0CBE8", "#59A14F", "#8CD17D", "#E15759", "#FF9D9A",
    "#F28E2B", "#B6992D", "#EDC948", "#B07AA1", "#76B7B2", "#9C755F"
]))

NEUT_ORDER = ["Neutrophil", "Neutrophil_PADI4", "Neutrophil_MX1", "Neutrophil_CXCR4", "Neutrophil_OSM"]
NEUT_LABEL = {
    "Neutrophil": "Unpolarized", "Neutrophil_PADI4": "PADI4+",
    "Neutrophil_MX1": "MX1+", "Neutrophil_CXCR4": "CXCR4+", "Neutrophil_OSM": "OSM+"
}
NEUT_COLOR = {
    "Neutrophil": "#8F8F8F", "Neutrophil_PADI4": "#7A5195", "Neutrophil_MX1": "#2A9D8F",
    "Neutrophil_CXCR4": "#4C78A8", "Neutrophil_OSM": "#E45756"
}

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 10,
    "axes.labelsize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.8, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"
})
sns.set_style("ticks")


def panel_label(ax, letter, x=-0.10, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=16, fontweight="bold",
            va="top", ha="left", clip_on=False)


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.clip(adjusted, 0, 1)
    return out


def save_all(fig, stem):
    for ext in ["png", "pdf", "svg"]:
        kwargs = {"dpi": 400} if ext == "png" else {}
        fig.savefig(OUT / f"{stem}.{ext}", bbox_inches="tight", facecolor="white", **kwargs)


def umap_panel(ax, data, color_col, order, colors, labels, title, max_points=None,
               label_centers=False, label_offsets=None):
    z = data
    if max_points and len(z) > max_points:
        z = z.sample(max_points, random_state=17)
    for category in order:
        q = z[z[color_col] == category]
        ax.scatter(q.UMAP_1, q.UMAP_2, s=1.15, c=colors[category], alpha=.58,
                   linewidths=0, rasterized=True, label=labels[category])
    if label_centers:
        centers = data.groupby(color_col)[["UMAP_1", "UMAP_2"]].median()
        for category in order:
            if category in centers.index:
                x, y = centers.loc[category]
                dx, dy = (label_offsets or {}).get(category, (0, 0))
                x += dx; y += dy
                ax.text(x, y, labels[category], ha="center", va="center", fontsize=6.3,
                        fontweight="bold", color="#202020",
                        bbox=dict(boxstyle="round,pad=.16", fc="white", ec="none", alpha=.78))
    ax.set_title(title, loc="left", fontweight="bold")
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


cells = pd.read_csv(CELL_FILE, usecols=[
    "cell_ID", "PatientID", "Diagnosis2", "cell_type", "PADI4", "CXCR4", "MX.1",
    "OSM", "CD66b", "CD15", "CD16", "CD11b"
], low_memory=False)
umap = pd.read_csv(UMAP_FILE, usecols=["cell_ID", "UMAP_1", "UMAP_2"])
dat = umap.merge(cells, on="cell_ID", how="inner", validate="one_to_one")
dat["broad_type"] = np.where(dat.cell_type.str.startswith("Neutrophil", na=False), "Neutrophil", dat.cell_type)

# ---------------- Figure 2—figure supplement 1: complete cell atlas ----------------
fig = plt.figure(figsize=(18, 15.8), facecolor="white")
outer = GridSpec(3, 2, figure=fig, height_ratios=[1.15, 1.05, .95], hspace=.34, wspace=.26,
                 left=.055, right=.98, top=.94, bottom=.055)

ax = fig.add_subplot(outer[0, 0]); panel_label(ax, "A", -.08, 1.04)
umap_panel(ax, dat, "broad_type", BROAD_ORDER, BROAD_COLOR, BROAD_LABEL,
           "CODEX atlas: complete cell annotations", max_points=90000, label_centers=True,
           label_offsets={"B Cell": (-.9, -.15), "CD4 T": (.85, -.15), "CD8 T": (.55, .25),
                          "TReg": (-.35, -.25), "Plasma B Cell": (-.25, .25)})

facet = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0, 1], wspace=.08)
facet_axes = []
for i, group in enumerate(GROUPS):
    a = fig.add_subplot(facet[0, i]); facet_axes.append(a)
    if i == 0: panel_label(a, "B", -.16, 1.04)
    umap_panel(a, dat[dat.Diagnosis2 == group], "broad_type", BROAD_ORDER, BROAD_COLOR, BROAD_LABEL,
               GROUP_LABEL[group], max_points=30000, label_centers=False)
    a.set_xlabel("UMAP 1")
    if i: a.set_ylabel("")
handles = [mpl.lines.Line2D([], [], marker='o', ls='', color=BROAD_COLOR[k], label=BROAD_LABEL[k], markersize=5) for k in BROAD_ORDER]
facet_axes[1].legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.12), ncol=4,
                     frameon=False, fontsize=6.5, handletextpad=.25, columnspacing=.75)

comp = (dat.groupby(["PatientID", "Diagnosis2", "broad_type"]).size().rename("n").reset_index())
tot = dat.groupby("PatientID").size().rename("total")
comp = comp.join(tot, on="PatientID")
comp["fraction"] = comp.n / comp.total
complete = pd.MultiIndex.from_product([dat.PatientID.unique(), BROAD_ORDER], names=["PatientID", "broad_type"])
comp = comp.set_index(["PatientID", "broad_type"]).reindex(complete, fill_value=0).reset_index()
patient_group = dat.drop_duplicates("PatientID").set_index("PatientID").Diagnosis2
comp["Diagnosis2"] = comp.PatientID.map(patient_group)
comp["total"] = comp.PatientID.map(tot); comp["fraction"] = comp.n / comp.total

ax = fig.add_subplot(outer[1, 0]); panel_label(ax, "C", -.08, 1.04)
heat = comp.groupby(["broad_type", "Diagnosis2"]).fraction.median().unstack().reindex(index=BROAD_ORDER, columns=GROUPS) * 100
sns.heatmap(heat, cmap="Blues", annot=True, fmt=".1f", linewidths=.5, linecolor="white",
            cbar_kws={"label": "Median cells (%)", "shrink": .78}, ax=ax)
ax.set_title("Median patient cell composition", loc="left", fontweight="bold")
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_xticklabels([GROUP_LABEL[g] for g in GROUPS], rotation=0)
ax.set_yticklabels([BROAD_LABEL[x] for x in BROAD_ORDER], rotation=0)

ax = fig.add_subplot(outer[1, 1]); panel_label(ax, "D", -.08, 1.04)
rows = []
for celltype in BROAD_ORDER:
    q = comp[comp.broad_type == celltype]
    vals = [q.loc[q.Diagnosis2 == g, "fraction"].values for g in GROUPS]
    _, p_kw = kruskal(*vals)
    infl = np.median(vals[2]); ctrl = np.median(vals[0])
    rows.append((celltype, np.log2((infl + 1e-4) / (ctrl + 1e-4)), p_kw))
stats = pd.DataFrame(rows, columns=["celltype", "log2fc", "p"])
stats["fdr"] = bh_adjust(stats.p)
stats = stats.set_index("celltype").loc[BROAD_ORDER].reset_index()
y = np.arange(len(stats))
sig = stats.fdr < .05
ax.axvline(0, color="#777777", lw=.9)
ax.scatter(stats.log2fc[~sig], y[~sig], s=35, color="#A9A9A9", edgecolor="white", lw=.5)
ax.scatter(stats.log2fc[sig], y[sig], s=48, color="#D95F5F", edgecolor="white", lw=.5)
for yi, row in stats[sig].iterrows():
    ax.text(row.log2fc + (.06 if row.log2fc >= 0 else -.06), yi, f"FDR {row.fdr:.3f}",
            ha="left" if row.log2fc >= 0 else "right", va="center", fontsize=6.6)
ax.set_yticks(y, [BROAD_LABEL[x] for x in stats.celltype]); ax.invert_yaxis()
ax.set_xlabel("log2 median fraction ratio: UC inflamed / control")
ax.set_title("Condition-associated remodeling of all annotations", loc="left", fontweight="bold")
clean(ax)

ax = fig.add_subplot(outer[2, :]); panel_label(ax, "E", -.04, 1.04)
patient_order = []
for group in GROUPS:
    patient_order += sorted(dat.loc[dat.Diagnosis2 == group, "PatientID"].unique())
patient_heat = comp.pivot(index="broad_type", columns="PatientID", values="fraction").reindex(index=BROAD_ORDER, columns=patient_order) * 100
sns.heatmap(patient_heat, cmap="mako", linewidths=.25, linecolor="white",
            cbar_kws={"label": "Cells (%)", "shrink": .65, "pad": .015}, ax=ax)
ax.set_title("Patient-resolved cell composition", loc="left", fontweight="bold")
ax.set_xlabel("Patients ordered by condition"); ax.set_ylabel("")
ax.set_yticklabels([BROAD_LABEL[x] for x in BROAD_ORDER], rotation=0)
for boundary in [sum(patient_group.eq(GROUPS[0])), sum(patient_group.isin(GROUPS[:2]))]:
    ax.axvline(boundary, color="white", lw=2.2)
starts = np.cumsum([0] + [patient_group.eq(g).sum() for g in GROUPS[:-1]])
sizes = [patient_group.eq(g).sum() for g in GROUPS]
for start, size, group in zip(starts, sizes, GROUPS):
    ax.text(start + size / 2, -1.05, GROUP_LABEL[group], ha="center", va="bottom",
            color=GROUP_COLOR[group], fontsize=8, fontweight="bold", clip_on=False)

fig.suptitle("Figure 2—figure supplement 1. Complete CODEX cell atlas across control and UC conditions",
             x=.055, y=.982, ha="left", fontsize=16, fontweight="bold")
fig.text(.055, .958, "245,020 cells from 24 patients; all abundance summaries and tests use patients as replicates.",
         fontsize=9, ha="left")
save_all(fig, "Figure_2_supplement_1_CODEX_full_annotations")
plt.close(fig)

# ---------------- Figure 2—figure supplement 2: neutrophil states ----------------
neu = dat[dat.cell_type.isin(NEUT_ORDER)].copy()
fig = plt.figure(figsize=(18, 16.7), facecolor="white")
outer = GridSpec(3, 2, figure=fig, height_ratios=[1.15, 1.08, 1.02], hspace=.36, wspace=.26,
                 left=.055, right=.98, top=.94, bottom=.055)

ax = fig.add_subplot(outer[0, 0]); panel_label(ax, "A", -.08, 1.04)
umap_panel(ax, neu, "cell_type", NEUT_ORDER, NEUT_COLOR, NEUT_LABEL,
           "Neutrophil states in the global CODEX embedding", label_centers=True,
           label_offsets={"Neutrophil": (-1.0, -.2), "Neutrophil_PADI4": (.8, -.1),
                          "Neutrophil_MX1": (-.75, -.55), "Neutrophil_CXCR4": (-.3, .45),
                          "Neutrophil_OSM": (-.85, .25)})

facet = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0, 1], wspace=.08)
facet_axes = []
for i, group in enumerate(GROUPS):
    a = fig.add_subplot(facet[0, i]); facet_axes.append(a)
    if i == 0: panel_label(a, "B", -.16, 1.04)
    umap_panel(a, neu[neu.Diagnosis2 == group], "cell_type", NEUT_ORDER, NEUT_COLOR, NEUT_LABEL,
               GROUP_LABEL[group], label_centers=False)
    if i: a.set_ylabel("")
handles = [mpl.lines.Line2D([], [], marker='o', ls='', color=NEUT_COLOR[k], label=NEUT_LABEL[k], markersize=5) for k in NEUT_ORDER]
facet_axes[1].legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.12), ncol=3,
                     frameon=False, fontsize=6.8, handletextpad=.25, columnspacing=.8)

patients = sorted(dat.PatientID.unique())
idx = pd.MultiIndex.from_product([patients, NEUT_ORDER], names=["PatientID", "cell_type"])
neuc = neu.groupby(["PatientID", "cell_type"]).size().rename("n").reindex(idx, fill_value=0).reset_index()
neuc["Diagnosis2"] = neuc.PatientID.map(patient_group)
neuc["total_cells"] = neuc.PatientID.map(tot)
neut_tot = neu.groupby("PatientID").size().rename("neut_total")
neuc["neut_total"] = neuc.PatientID.map(neut_tot).fillna(0)
neuc["fraction_all"] = neuc.n / neuc.total_cells
neuc["fraction_neut"] = np.where(neuc.neut_total > 0, neuc.n / neuc.neut_total, 0)

sub = GridSpecFromSubplotSpec(1, 5, subplot_spec=outer[1, :], wspace=.34)
pvals = []
for subtype in NEUT_ORDER:
    q = neuc[neuc.cell_type == subtype]
    vals = [q.loc[q.Diagnosis2 == g, "fraction_all"].values for g in GROUPS]
    _, p = kruskal(*vals); pvals.append(p)
fdrs = bh_adjust(pvals)
for i, (subtype, fdr) in enumerate(zip(NEUT_ORDER, fdrs)):
    ax = fig.add_subplot(sub[0, i])
    if i == 0: panel_label(ax, "C", -.30, 1.04)
    q = neuc[neuc.cell_type == subtype].copy(); q["Group"] = q.Diagnosis2.map(GROUP_LABEL)
    order = [GROUP_LABEL[g] for g in GROUPS]
    sns.boxplot(data=q, x="Group", y="fraction_all", order=order,
                color=NEUT_COLOR[subtype], width=.58, showfliers=False, linewidth=1, ax=ax)
    sns.stripplot(data=q, x="Group", y="fraction_all", order=order,
                  color="#222222", size=3.6, jitter=.13, ax=ax)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1, decimals=1))
    ax.set_title(NEUT_LABEL[subtype], fontweight="bold")
    ax.set_xlabel(""); ax.tick_params(axis="x", rotation=28)
    ax.set_ylabel("Cells (% of all cells)" if i == 0 else "")
    ax.text(.5, .97, f"KW FDR = {fdr:.3g}", transform=ax.transAxes, ha="center", va="top", fontsize=6.8)
    clean(ax)

ax = fig.add_subplot(outer[2, 0]); panel_label(ax, "D", -.08, 1.04)
stack = neuc.groupby(["Diagnosis2", "cell_type"]).n.sum().unstack().reindex(index=GROUPS, columns=NEUT_ORDER)
stack = stack.div(stack.sum(axis=1), axis=0) * 100
bottom = np.zeros(len(GROUPS))
for subtype in NEUT_ORDER:
    vals = stack[subtype].values
    ax.bar(np.arange(len(GROUPS)), vals, bottom=bottom, color=NEUT_COLOR[subtype], width=.66,
           label=NEUT_LABEL[subtype])
    bottom += vals
ax.set_xticks(range(3), [GROUP_LABEL[g] for g in GROUPS]); ax.set_ylim(0, 100)
ax.set_ylabel("Neutrophil-state composition (%)")
ax.set_title("Composition of the neutrophil compartment", loc="left", fontweight="bold")
ax.legend(frameon=False, ncol=2, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0))
clean(ax)

sub = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[2, 1], width_ratios=[1.08, 1.0], wspace=.35)
ax = fig.add_subplot(sub[0, 0]); panel_label(ax, "E", -.16, 1.04)
marker_cols = ["PADI4", "CXCR4", "MX.1", "OSM", "CD16", "CD11b", "CD15", "CD66b"]
marker_labels = ["PADI4", "CXCR4", "MX1", "OSM", "CD16", "CD11b", "CD15", "CD66b"]
med = neu.groupby("cell_type")[marker_cols].median().reindex(NEUT_ORDER)
zmed = (med - med.mean(axis=0)) / med.std(axis=0, ddof=0).replace(0, 1)
sns.heatmap(zmed, cmap="vlag", center=0, vmin=-1.8, vmax=1.8, annot=True, fmt=".1f",
            linewidths=.5, linecolor="white", cbar_kws={"label": "Marker z-score", "shrink": .7}, ax=ax)
ax.set_title("Phenotypic validation", loc="left", fontweight="bold")
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_xticklabels(marker_labels, rotation=45, ha="right")
ax.set_yticklabels([NEUT_LABEL[x] for x in NEUT_ORDER], rotation=0)

ax = fig.add_subplot(sub[0, 1]); panel_label(ax, "F", -.16, 1.04)
neut_patient_order = []
for group in GROUPS:
    neut_patient_order += sorted(dat.loc[dat.Diagnosis2 == group, "PatientID"].unique())
ph = neuc.pivot(index="cell_type", columns="PatientID", values="fraction_neut").reindex(index=NEUT_ORDER, columns=neut_patient_order) * 100
sns.heatmap(ph, cmap="rocket_r", linewidths=.3, linecolor="white",
            cbar_kws={"label": "Within-neutrophil fraction (%)", "shrink": .7}, ax=ax)
ax.set_title("Patient-resolved states", loc="left", fontweight="bold")
ax.set_xlabel("Patients ordered by condition"); ax.set_ylabel("")
ax.set_yticklabels([NEUT_LABEL[x] for x in NEUT_ORDER], rotation=0)
for boundary in [patient_group.eq(GROUPS[0]).sum(), patient_group.isin(GROUPS[:2]).sum()]:
    ax.axvline(boundary, color="white", lw=2)

fig.suptitle("Figure 2—figure supplement 2. CODEX neutrophil states across control and UC conditions",
             x=.055, y=.982, ha="left", fontsize=16, fontweight="bold")
fig.text(.055, .958, f"{len(neu):,} neutrophils from 24 patients; abundance tests use patient-level fractions.",
         fontsize=9, ha="left")
save_all(fig, "Figure_2_supplement_2_CODEX_neutrophil_subtypes")
plt.close(fig)

manifest = {
    "source_cells": str(CELL_FILE), "source_umap": str(UMAP_FILE),
    "n_cells": int(len(dat)), "n_patients": int(dat.PatientID.nunique()),
    "n_neutrophils": int(len(neu)), "groups": GROUPS,
    "broad_annotations": BROAD_ORDER, "neutrophil_states": NEUT_ORDER,
    "outputs": [
        "Figure_2_supplement_1_CODEX_full_annotations.[png|pdf|svg]",
        "Figure_2_supplement_2_CODEX_neutrophil_subtypes.[png|pdf|svg]"
    ]
}
(OUT / "Figure_2_CODEX_annotation_supplements_manifest.json").write_text(json.dumps(manifest, indent=2))
