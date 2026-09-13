from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

DATA=Path("giotto_codex_results/nolan_neighborhoods/cells_with_nolan_neighborhoods.csv")
FIG=Path("giotto_codex_results/figures")
BG="#07090D"; MUTED="#303640"
cells=pd.read_csv(DATA)
groups=[
    ("Control","P03",[3,6,7,8]),
    ("UC noninflamed","P17",[1,5,10]),
    ("UC inflamed","P07",[2,4,9]),
]
palette={1:"#00C8FF",2:"#FF7A00",3:"#35E65C",4:"#FF3154",5:"#A66BFF",
         6:"#D98C5F",7:"#FF5EC4",8:"#B9C0CA",9:"#F1E400",10:"#00E0C6"}
descriptions={
  1:"mixed stromal",2:"PADI4/CXCR4/MX1 neutrophil",3:"enteroendocrine/CD4/B",
  4:"epithelial/neutrophil/CD8",5:"macrophage/DC/endothelial",6:"CD4/enteroendocrine/B",
  7:"CD8/neutrophil/epithelial",8:"TReg/B/enteroendocrine",9:"plasma B/OSM-neutrophil/DC",
  10:"fibroblast-dominant stromal"
}

plt.style.use("dark_background")
fig,axes=plt.subplots(1,3,figsize=(20,8)); fig.patch.set_facecolor(BG)
for ax,(label,pid,enriched) in zip(axes,groups):
    d=cells[cells.PatientID==pid].copy(); ax.set_facecolor(BG)
    other=~d.neighborhood_id.isin(enriched)
    ax.scatter(d.loc[other,"x"],d.loc[other,"y"],s=.62,c=MUTED,alpha=.42,linewidths=0,rasterized=True)
    for cn in enriched:
        z=d[d.neighborhood_id==cn]
        ax.scatter(z.x,z.y,s=1.35,c=palette[cn],alpha=.98,linewidths=0,rasterized=True)
    ax.set_aspect("equal",adjustable="box"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(f"{label} ({pid})",fontweight="bold",fontsize=16,pad=12)
    handles=[Line2D([0],[0],marker="o",linestyle="",markersize=8,markerfacecolor=palette[n],
                    markeredgecolor="none",label=f"CN{n}: {descriptions[n]}") for n in enriched]
    ax.legend(handles=handles,loc="lower center",bbox_to_anchor=(.5,-.13),frameon=False,
              fontsize=9,ncol=1,handletextpad=.5)
fig.suptitle("Representative CODEX maps of group-increased cellular neighborhoods",
             fontweight="bold",fontsize=21,y=.98)
fig.tight_layout(rect=[0,.06,1,.94],w_pad=2)
fig.savefig(FIG/"Figure_20_group_increased_CODEX_neighborhoods.png",dpi=300,bbox_inches="tight",facecolor=BG)
fig.savefig(FIG/"Figure_20_group_increased_CODEX_neighborhoods.pdf",bbox_inches="tight",facecolor=BG)
plt.close(fig)
print("Saved Figure 20")
