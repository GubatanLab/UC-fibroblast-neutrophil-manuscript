from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal
from sklearn.cluster import MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns

SEED = 0
K = 10
N_NEIGHBORHOODS = 10
GROUPS = ["Control", "UC_Noninflamed", "UC_Inflamed"]
COLORS = {"Control": "#4C78A8", "UC_Noninflamed": "#F2CF5B", "UC_Inflamed": "#E45756"}
OUT = Path("giotto_codex_results/nolan_neighborhoods")
FIG = Path("giotto_codex_results/figures")
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

cells = pd.read_csv(OUT / "nolan_input_cells.csv")
cells = cells[cells["Diagnosis2"].isin(GROUPS)].reset_index(drop=True)
cell_types = sorted(cells["cell_type"].astype(str).unique())
type_to_i = {v: i for i, v in enumerate(cell_types)}
type_idx = cells["cell_type"].map(type_to_i).to_numpy()
windows = np.zeros((len(cells), len(cell_types)), dtype=np.uint8)

# Nolan et al. notebook implementation: composition of the 10 nearest cells,
# computed independently inside each image/region; the index cell is included.
for region, idx in cells.groupby("region", sort=False).groups.items():
    idx = np.asarray(list(idx), dtype=int)
    xy = cells.loc[idx, ["x", "y"]].to_numpy()
    n_neighbors = min(K, len(idx))
    neighbors = NearestNeighbors(n_neighbors=n_neighbors).fit(xy).kneighbors(return_distance=False)
    neighbor_types = type_idx[idx[neighbors]]
    for local_row, global_row in enumerate(idx):
        windows[global_row] = np.bincount(neighbor_types[local_row], minlength=len(cell_types))

km = MiniBatchKMeans(n_clusters=N_NEIGHBORHOODS, random_state=SEED, batch_size=4096, n_init=20)
labels = km.fit_predict(windows)
centroids = km.cluster_centers_

# Stable, interpretable names based on the most overrepresented cell types.
global_prop = np.maximum((windows.sum(axis=0) / windows.sum()), 1e-9)
centroid_prop = centroids / np.maximum(centroids.sum(axis=1, keepdims=True), 1e-9)
log2_enrichment = np.log2((centroid_prop + 1e-4) / (global_prop + 1e-4))
names = []
for n in range(N_NEIGHBORHOODS):
    top = np.argsort(log2_enrichment[n])[::-1][:3]
    names.append(f"CN{n+1}: " + " / ".join(cell_types[i] for i in top))
label_to_name = dict(enumerate(names))
cells["neighborhood_id"] = labels + 1
cells["neighborhood"] = pd.Series(labels).map(label_to_name)

profile = pd.DataFrame(centroid_prop, index=names, columns=cell_types)
enrichment = pd.DataFrame(log2_enrichment, index=names, columns=cell_types)
profile.to_csv(OUT / "neighborhood_celltype_proportions.csv")
enrichment.to_csv(OUT / "neighborhood_celltype_log2_enrichment.csv")
pd.DataFrame({"neighborhood_id": range(1, N_NEIGHBORHOODS + 1), "neighborhood": names}).to_csv(
    OUT / "neighborhood_names.csv", index=False
)
cells.to_csv(OUT / "cells_with_nolan_neighborhoods.csv", index=False)

counts = cells.groupby(["PatientID", "Diagnosis2", "neighborhood"], observed=False).size().rename("cell_count").reset_index()
totals = cells.groupby(["PatientID", "Diagnosis2"]).size().rename("total_cells").reset_index()
patients = totals[["PatientID", "Diagnosis2"]].drop_duplicates()
grid = patients.assign(key=1).merge(pd.DataFrame({"neighborhood": names, "key": 1}), on="key").drop(columns="key")
freq = grid.merge(counts, how="left").merge(totals, on=["PatientID", "Diagnosis2"])
freq["cell_count"] = freq["cell_count"].fillna(0).astype(int)
freq["frequency"] = freq["cell_count"] / freq["total_cells"]
freq.to_csv(OUT / "neighborhood_frequency_by_patient.csv", index=False)

overall = []
pairwise = []
contrasts = [("UC_Noninflamed", "Control"), ("UC_Inflamed", "Control"), ("UC_Inflamed", "UC_Noninflamed")]
for neighborhood in names:
    sub = freq[freq.neighborhood == neighborhood]
    vals = [sub.loc[sub.Diagnosis2 == g, "frequency"].to_numpy() for g in GROUPS]
    overall.append({"neighborhood": neighborhood, "kruskal_p": kruskal(*vals).pvalue,
                    **{f"median_{g}": float(np.median(v)) for g, v in zip(GROUPS, vals)}})
    for a, b in contrasts:
        va = sub.loc[sub.Diagnosis2 == a, "frequency"].to_numpy()
        vb = sub.loc[sub.Diagnosis2 == b, "frequency"].to_numpy()
        pairwise.append({"neighborhood": neighborhood, "contrast": f"{a} vs {b}",
                         "group_1": a, "group_2": b, "n_group_1": len(va), "n_group_2": len(vb),
                         "median_group_1": np.median(va), "median_group_2": np.median(vb),
                         "median_difference": np.median(va) - np.median(vb),
                         "p_value": mannwhitneyu(va, vb, alternative="two-sided").pvalue})
overall = pd.DataFrame(overall)
overall["p_adj_BH"] = multipletests(overall.kruskal_p, method="fdr_bh")[1]
pairwise = pd.DataFrame(pairwise)
pairwise["p_adj_BH"] = pairwise.groupby("contrast")["p_value"].transform(
    lambda x: multipletests(x, method="fdr_bh")[1]
)
overall.to_csv(OUT / "statistics_neighborhood_overall_kruskal.csv", index=False)
pairwise.to_csv(OUT / "statistics_neighborhood_pairwise_wilcoxon.csv", index=False)

sns.set_theme(style="whitegrid", context="talk")
g = sns.clustermap(enrichment, cmap="vlag", center=0, vmin=-3, vmax=3,
                   row_cluster=False, figsize=(15, 9), cbar_kws={"label": "log2 enrichment"})
g.ax_heatmap.set_xlabel("Cell type")
g.ax_heatmap.set_ylabel("Cellular neighborhood")
g.fig.suptitle("Cellular-neighborhood composition (10 nearest cells)", y=1.02, fontweight="bold")
g.savefig(FIG / "Figure_08_nolan_neighborhood_composition.png", dpi=300, bbox_inches="tight")
g.savefig(FIG / "Figure_08_nolan_neighborhood_composition.pdf", bbox_inches="tight")
plt.close(g.fig)

fig, axes = plt.subplots(2, 5, figsize=(22, 10), sharey=False)
for ax, neighborhood in zip(axes.flat, names):
    sub = freq[freq.neighborhood == neighborhood]
    sns.boxplot(data=sub, x="Diagnosis2", y="frequency", order=GROUPS, palette=COLORS,
                width=.65, showfliers=False, ax=ax)
    sns.stripplot(data=sub, x="Diagnosis2", y="frequency", order=GROUPS, color="black", size=7, linewidth=.45, edgecolor="white", ax=ax)
    ax.set_title(neighborhood, fontsize=10, fontweight="bold")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=35, labelsize=9)
    ax.set_ylabel("Patient-level frequency" if ax in axes[:, 0] else "")
    sig = pairwise[(pairwise.neighborhood == neighborhood) & (pairwise.p_adj_BH < 0.05)].copy()
    if len(sig):
        ymax = sub["frequency"].max()
        ymin = sub["frequency"].min()
        step = max((ymax - ymin) * 0.12, ymax * 0.06, 0.008)
        base = ymax + step * 0.35
        xpos = {g: i for i, g in enumerate(GROUPS)}
        for j, row in enumerate(sig.itertuples(index=False)):
            x1, x2 = xpos[row.group_1], xpos[row.group_2]
            y = base + j * step
            h = step * 0.22
            ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="black", linewidth=1.8)
            stars = "***" if row.p_adj_BH < 0.001 else "**" if row.p_adj_BH < 0.01 else "*"
            ax.text((x1 + x2) / 2, y + h, stars, ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.set_ylim(top=base + len(sig) * step + step * 0.25)
fig.suptitle("Cellular neighborhoods across tissue groups", fontweight="bold")
fig.tight_layout(rect=[0, 0, 1, .96])
fig.savefig(FIG / "Figure_09_nolan_neighborhood_group_comparison.png", dpi=300, bbox_inches="tight")
fig.savefig(FIG / "Figure_09_nolan_neighborhood_group_comparison.pdf", bbox_inches="tight")
plt.close(fig)

print(f"Analyzed {len(cells):,} cells from {cells.PatientID.nunique()} patients")
print(pairwise.sort_values("p_adj_BH").head(15).to_string(index=False))
