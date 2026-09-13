from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
from scipy.spatial import cKDTree

DATA=Path("giotto_codex_results/nolan_neighborhoods/cells_with_nolan_neighborhoods.csv")
FIG=Path("giotto_codex_results/figures")
BG="#07090D"; MUTED="#272C35"
cells=pd.read_csv(DATA)
palette={1:"#00C8FF",2:"#FF7A00",3:"#35E65C",4:"#FF3154",5:"#A66BFF",
         6:"#D98C5F",7:"#FF5EC4",8:"#D7DEE8",9:"#F1E400",10:"#00E0C6"}
panels=[
    ("Control","P03",[3,6,7,8],[(7,8),(3,6)]),
    ("UC noninflamed","P17",[1,5,10],[(10,1),(10,5)]),
    ("UC inflamed","P07",[2,4,9],[(2,4),(2,9)]),
]

def best_crop(d, selected, half=430):
    s=d[d.neighborhood_id.isin(selected)]
    xy=s[["x","y"]].values
    if len(xy)==0: return d.x.median(),d.y.median(),half
    sample=xy[np.linspace(0,len(xy)-1,min(len(xy),700),dtype=int)]
    tree=cKDTree(xy); counts=np.array([len(tree.query_ball_point(p,half)) for p in sample])
    cx,cy=sample[counts.argmax()]
    return cx,cy,half

def closest_pair(d,a,b):
    da=d[d.neighborhood_id==a]; db=d[d.neighborhood_id==b]
    if len(da)==0 or len(db)==0: return None
    tree=cKDTree(db[["x","y"]].values); dist,j=tree.query(da[["x","y"]].values)
    i=int(np.argmin(dist)); return da.iloc[i][["x","y"]].values.astype(float),db.iloc[int(j[i])][["x","y"]].values.astype(float)

plt.style.use("dark_background")
fig,axes=plt.subplots(1,3,figsize=(20,7.5)); fig.patch.set_facecolor(BG)
for ax,(label,pid,selected,interfaces) in zip(axes,panels):
    full=cells[cells.PatientID==pid].copy(); cx,cy,h=best_crop(full,selected)
    d=full[(full.x.between(cx-h,cx+h))&(full.y.between(cy-h,cy+h))]
    ax.set_facecolor(BG)
    other=~d.neighborhood_id.isin(selected)
    ax.scatter(d.loc[other,"x"],d.loc[other,"y"],s=2.0,c=MUTED,alpha=.48,linewidths=0,rasterized=True)
    for cn in selected:
        z=d[d.neighborhood_id==cn]
        ax.scatter(z.x,z.y,s=4.0,c=palette[cn],alpha=.98,linewidths=0,rasterized=True)
    for k,(a,b) in enumerate(interfaces):
        pair=closest_pair(d,a,b)
        if pair is None: continue
        p1,p2=pair
        arrow=FancyArrowPatch(p1,p2,arrowstyle="<->",mutation_scale=19,linewidth=3.2,
                              color="white",path_effects=[],zorder=8)
        ax.add_patch(arrow)
        mid=(p1+p2)/2; offset=np.array([0,28 if k==0 else -32])
        ax.text(*(mid+offset),f"CN{a} ↔ CN{b}",ha="center",va="center",color="white",
                fontsize=10,fontweight="bold",bbox=dict(boxstyle="round,pad=.25",facecolor=BG,
                edgecolor="#B9C0CA",alpha=.90),zorder=9)
    ax.set_xlim(cx-h,cx+h); ax.set_ylim(cy+h,cy-h); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_color("#59616D"); sp.set_linewidth(1)
    ax.set_title(f"{label} ({pid})",fontweight="bold",fontsize=16,pad=10)
    handles=[Line2D([0],[0],marker="o",linestyle="",markersize=8,markerfacecolor=palette[n],
                    markeredgecolor="none",label=f"CN{n}") for n in selected]
    ax.legend(handles=handles,loc="lower right",frameon=False,ncol=2,fontsize=9)
fig.suptitle("Closeups of group-associated cellular-neighborhood interfaces",fontweight="bold",fontsize=20,y=.98)
fig.text(.5,.025,"Double-headed arrows mark representative observed spatial interfaces; they do not imply signaling direction.",
         ha="center",color="#D7DEE8",fontsize=11)
fig.tight_layout(rect=[0,.06,1,.94],w_pad=2)
fig.savefig(FIG/"Figure_21_CODEX_neighborhood_closeups_interactions.png",dpi=300,bbox_inches="tight",facecolor=BG)
fig.savefig(FIG/"Figure_21_CODEX_neighborhood_closeups_interactions.pdf",bbox_inches="tight",facecolor=BG)
plt.close(fig)
print("Saved Figure 21")
