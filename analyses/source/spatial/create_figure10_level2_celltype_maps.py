from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path("giotto_codex_results")
FIG=ROOT/"figures"
cells=pd.read_csv(ROOT/"nolan_neighborhoods/cells_with_nolan_neighborhoods.csv")

samples=[("Control","P03"),("UC noninflamed","P17"),("UC inflamed","P07")]
colors={
    "B Cell":"#00A6FF", "CD4 T":"#5B8CFF", "CD8 T":"#FF7A00",
    "Dendritic Cell":"#35E65C", "Endothelial Cell":"#8EF07C",
    "Enteroendocrine Cell":"#FF3154", "Epithelial Cell":"#A66BFF",
    "Fibroblast":"#C89BFF", "Macrophage":"#FF5EC4", "Neutrophil":"#FF9ACD",
    "Neutrophil_CXCR4":"#B9C0CA", "Neutrophil_MX1":"#F1E400",
    "Neutrophil_OSM":"#E8E27A", "Neutrophil_PADI4":"#00E0C6",
    "Plasma B Cell":"#73E6FF", "TReg":"#D98C5F"
}

def scalebar(ax,length_um=500):
    span=abs(ax.get_xlim()[1]-ax.get_xlim()[0]); frac=min(length_um/span,.34); x0,y=.065,.065
    ax.plot([x0,x0+frac],[y,y],transform=ax.transAxes,color="white",lw=4,
            solid_capstyle="butt",clip_on=False,zorder=20)
    ax.text(x0+frac/2,y+.022,f"{length_um} µm",transform=ax.transAxes,color="white",
            ha="center",va="bottom",fontsize=10,fontweight="bold",zorder=20)

fig,axes=plt.subplots(1,3,figsize=(18,7.4),facecolor="#07090D",gridspec_kw={"wspace":.045})
for ax,(group,pid) in zip(axes,samples):
    dat=cells[cells.PatientID==pid]
    ax.set_facecolor("#07090D")
    for celltype in colors:
        q=dat[dat.cell_type==celltype]
        if len(q):
            ax.scatter(q.x,q.y,s=1.25,c=colors[celltype],alpha=.96,linewidths=0,rasterized=True)
    xmin,xmax=dat.x.min(),dat.x.max(); ymin,ymax=dat.y.min(),dat.y.max()
    span=max(xmax-xmin,ymax-ymin)*1.04; xm,ym=(xmin+xmax)/2,(ymin+ymax)/2
    ax.set_xlim(xm-span/2,xm+span/2); ax.set_ylim(ym-span/2,ym+span/2)
    ax.set_aspect("equal"); ax.set_box_aspect(1); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{group} ({pid})",color="white",fontsize=16,fontweight="bold",pad=8)
    for spine in ax.spines.values(): spine.set_visible(True); spine.set_color("#4B525D"); spine.set_linewidth(1.2)
    scalebar(ax,500)

handles=[Line2D([0],[0],marker="o",linestyle="",markerfacecolor=colors[c],markeredgecolor="none",markersize=7,label=c)
         for c in colors if c in set(cells.cell_type)]
leg=fig.legend(handles=handles,loc="lower center",bbox_to_anchor=(.5,.012),ncol=4,frameon=False,
               labelcolor="white",title="Level 2 cell-type annotation",title_fontsize=11,fontsize=9,
               columnspacing=1.5,handletextpad=.5)
leg.get_title().set_color("white")
fig.suptitle("Representative CODEX tissue maps colored by level 2 cell type",
             color="white",fontsize=20,fontweight="bold",y=.985)
fig.tight_layout(rect=[.01,.15,.99,.94])
out=FIG/"Figure_10_representative_CODEX_neighborhood_maps"
fig.savefig(out.with_suffix(".png"),dpi=300,bbox_inches="tight",facecolor="#07090D")
fig.savefig(out.with_suffix(".pdf"),bbox_inches="tight",facecolor="#07090D")
plt.close(fig)
print("Saved Figure 10 level-2 cell-type CODEX maps")
