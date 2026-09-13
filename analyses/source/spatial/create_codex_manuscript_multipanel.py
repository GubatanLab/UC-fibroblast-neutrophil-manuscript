from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from scipy.spatial import cKDTree
import seaborn as sns

ROOT=Path("giotto_codex_results"); FIG=ROOT/"figures"; MP=ROOT/"multipanel"
umap=pd.read_csv(MP/"umap_cells.csv")
cells=pd.read_csv(ROOT/"nolan_neighborhoods/cells_with_nolan_neighborhoods.csv")
freq=pd.read_csv(ROOT/"nolan_neighborhoods/neighborhood_frequency_by_patient.csv")
stats=pd.read_csv(ROOT/"nolan_neighborhoods/statistics_neighborhood_pairwise_wilcoxon.csv")
enrich=pd.read_csv(ROOT/"nolan_neighborhoods/neighborhood_celltype_log2_enrichment.csv",index_col=0)
nf=pd.read_csv(ROOT/"neutrophil_fibroblast_subtypes/tables/neutrophil_subtype_fibroblast_enrichment_by_patient.csv")

groups=["Control","UC_Noninflamed","UC_Inflamed"]
glabels={"Control":"Control","UC_Noninflamed":"UC noninflamed","UC_Inflamed":"UC inflamed"}
gcols={"Control":"#4C78A8","UC_Noninflamed":"#E6BE4A","UC_Inflamed":"#D95F5F"}
ctypes=sorted(umap.Celltype_updated.dropna().unique()); cmap=plt.get_cmap("tab20",len(ctypes)); ctcols={c:cmap(i) for i,c in enumerate(ctypes)}
cncols={1:"#00C8FF",2:"#FF7A00",3:"#35E65C",4:"#FF3154",5:"#A66BFF",6:"#D98C5F",7:"#FF5EC4",8:"#9AA3AE",9:"#D9C900",10:"#00B9A8"}

sns.set_theme(style="whitegrid",context="talk")
plt.rcParams.update({"font.family":"Arial","font.weight":"medium","axes.titleweight":"bold",
                     "axes.linewidth":1.4,"lines.linewidth":2.3})
fig=plt.figure(figsize=(24,20),facecolor="white")
gs=GridSpec(4,12,figure=fig,height_ratios=[5.2,4.2,4.4,4.4],hspace=.40,wspace=.50,
            top=.965,bottom=.035,left=.045,right=.975)

def letter(ax,s,color="black"):
    ax.text(-.08,1.06,s,transform=ax.transAxes,fontsize=24,fontweight="bold",va="top",ha="left",color=color)

def add_scalebar(ax,length_um,color="white"):
    xmin,xmax=ax.get_xlim(); span=abs(xmax-xmin)
    frac=min(length_um/span,.35); x0,y=.07,.075
    ax.plot([x0,x0+frac],[y,y],transform=ax.transAxes,color=color,lw=4,
            solid_capstyle="butt",clip_on=False,zorder=30)
    ax.text(x0+frac/2,y+.025,f"{length_um} µm",transform=ax.transAxes,color=color,
            ha="center",va="bottom",fontsize=10,fontweight="bold",zorder=30)

# A. UMAP and patient-balanced cell-type composition
suba=GridSpecFromSubplotSpec(1,2,subplot_spec=gs[0,:8],width_ratios=[1,.72],wspace=.48)
ax=fig.add_subplot(suba[0,0])
ax.text(-.15,1.10,"A",transform=ax.transAxes,fontsize=24,fontweight="bold",va="top",ha="left")
for c in ctypes:
    d=umap[umap.Celltype_updated==c]
    ax.scatter(d.UMAP_1,d.UMAP_2,s=.42,c=[ctcols[c]],alpha=.72,linewidths=0,rasterized=True,label=c)
ax.set_title("Harmony UMAP of annotated CODEX cells",fontsize=18); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.grid(False)
ax.legend(loc="center left",bbox_to_anchor=(1.015,.50),frameon=False,fontsize=10,
          markerscale=5.5,ncol=1,labelspacing=.25,handletextpad=.34,borderaxespad=0)

ax=fig.add_subplot(suba[0,1])
pc=cells.groupby(["PatientID","Diagnosis2","cell_type"]).size().rename("n").reset_index()
pt=cells.groupby(["PatientID","Diagnosis2"]).size().rename("total").reset_index(); pc=pc.merge(pt); pc["fraction"]=pc.n/pc.total
grid=pd.MultiIndex.from_product([pc.PatientID.unique(),ctypes],names=["PatientID","cell_type"]).to_frame(index=False)
grid=grid.merge(cells[["PatientID","Diagnosis2"]].drop_duplicates(),on="PatientID").merge(pc[["PatientID","cell_type","fraction"]],how="left").fillna({"fraction":0})
comp=grid.groupby(["Diagnosis2","cell_type"]).fraction.mean().unstack(fill_value=0).loc[groups]
comp=comp.div(comp.sum(axis=1),axis=0); ctorder=comp.mean().sort_values(ascending=False).index.tolist(); comp=comp[ctorder]
bottom=np.zeros(len(groups)); x=np.arange(len(groups))
for c in ctorder:
    ax.bar(x,comp[c].values,bottom=bottom,color=ctcols[c],width=.72,edgecolor="white",linewidth=.25)
    bottom+=comp[c].values
ax.set_xticks(x,["Control","UC\nnoninflamed","UC\ninflamed"],rotation=0,ha="center"); ax.set_ylim(0,1)
ax.set_ylabel("Mean patient cell fraction"); ax.set_xlabel("")
ax.set_title("Cell-type composition",fontsize=18)

# B. CN1-CN10 neighborhood-composition heatmap
ax=fig.add_subplot(gs[0,8:]); letter(ax,"B")
enrich=enrich.loc[sorted(enrich.index,key=lambda x:int(x.split(':')[0][2:]))]
sns.heatmap(enrich,cmap="vlag",center=0,vmin=-3,vmax=3,ax=ax,cbar_kws={"label":"log2 enrichment"})
ax.set_title("Cellular-neighborhood composition",fontsize=18); ax.set_xlabel("Cell type"); ax.set_ylabel("")
ax.tick_params(axis="x",rotation=55,labelsize=8); ax.tick_params(axis="y",labelsize=8)
ax.set_yticklabels([f"CN{i}" for i in range(1,11)],rotation=0)

# C. focused abundance
subgs=GridSpecFromSubplotSpec(1,3,subplot_spec=gs[1,:],wspace=.30)
focus=[(7,"CD8 T / neutrophil / epithelial"),(10,"Fibroblast-dominant stromal"),(2,"PADI4 / CXCR4 / MX1 neutrophil")]
for j,(cn,title) in enumerate(focus):
    ax=fig.add_subplot(subgs[0,j]);
    if j==0: letter(ax,"C")
    d=freq[freq.neighborhood.str.startswith(f"CN{cn}:")].copy(); d["Group"]=d.Diagnosis2.map(glabels)
    orderg=[glabels[g] for g in groups]; pal={glabels[k]:v for k,v in gcols.items()}
    sns.boxplot(d,x="Group",y="frequency",order=orderg,palette=pal,width=.62,linewidth=1.7,showfliers=False,ax=ax)
    sns.stripplot(d,x="Group",y="frequency",order=orderg,color="black",size=7,linewidth=.45,edgecolor="white",ax=ax)
    ax.set_title(f"CN{cn}: {title}",fontsize=15); ax.set_xlabel(""); ax.set_ylabel("Patient-level frequency" if j==0 else ""); ax.tick_params(axis="x",rotation=25,labelsize=10)
    sig=stats[(stats.neighborhood.str.startswith(f"CN{cn}:"))&(stats.p_adj_BH<.05)].sort_values("p_adj_BH")
    ymax=d.frequency.max(); step=max(ymax*.10,.012); base=ymax+step*.2; xpos={g:i for i,g in enumerate(groups)}
    for k,r in enumerate(sig.itertuples()):
        x1,x2=xpos[r.group_1],xpos[r.group_2]; y=base+k*step; h=step*.18
        ax.plot([x1,x1,x2,x2],[y,y+h,y+h,y],c="black",lw=1.7)
        star="***" if r.p_adj_BH<.001 else "**" if r.p_adj_BH<.01 else "*"
        ax.text((x1+x2)/2,y+h,star,ha="center",va="bottom",fontweight="bold")
    if len(sig): ax.set_ylim(top=base+(len(sig)+.4)*step)

# D. whole-tissue maps
subgs=GridSpecFromSubplotSpec(1,3,subplot_spec=gs[2,:],wspace=.12)
mapsets=[("Control","P03",[3,6,7,8]),("UC noninflamed","P17",[1,5,10]),("UC inflamed","P07",[2,4,9])]
for j,(lab,pid,sel) in enumerate(mapsets):
    ax=fig.add_subplot(subgs[0,j],facecolor="#000000");
    if j==0: ax.text(.015,.98,"D",transform=ax.transAxes,color="white",fontsize=24,fontweight="bold",va="top",zorder=40)
    d=cells[cells.PatientID==pid]; other=~d.neighborhood_id.isin(sel)
    ax.scatter(d.loc[other,"x"],d.loc[other,"y"],s=.65,c="#343A46",alpha=.72,linewidths=0,rasterized=True)
    for n in sel:
        z=d[d.neighborhood_id==n]; ax.scatter(z.x,z.y,s=1.4,c=cncols[n],alpha=.98,linewidths=0,rasterized=True)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"{lab} ({pid})",color="black",fontsize=16)
    add_scalebar(ax,500)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color("#69717D"); sp.set_linewidth(1.7)
    handles=[Line2D([0],[0],marker="o",linestyle="",markerfacecolor=cncols[n],markeredgecolor="none",markersize=8,label=f"CN{n}") for n in sel]
    ax.legend(handles=handles,loc="lower center",ncol=len(sel),frameon=False,labelcolor="white",fontsize=9)

# E. closeups
def crop(d,sel,h=420):
    z=d[d.neighborhood_id.isin(sel)]; xy=z[["x","y"]].values; sample=xy[np.linspace(0,len(xy)-1,min(600,len(xy)),dtype=int)]; t=cKDTree(xy); p=sample[np.argmax([len(t.query_ball_point(q,h)) for q in sample])]; return p[0],p[1],h
def nearest(d,a,b):
    da=d[d.neighborhood_id==a]; db=d[d.neighborhood_id==b]
    if not len(da) or not len(db): return None
    t=cKDTree(db[["x","y"]]); dist,j=t.query(da[["x","y"]]); i=np.argmin(dist); return da.iloc[i][["x","y"]].values.astype(float),db.iloc[j[i]][["x","y"]].values.astype(float)
subgs=GridSpecFromSubplotSpec(1,3,subplot_spec=gs[3,:],wspace=.12)
close=[("Control","P03",[3,6,7,8],[(7,8)]),("UC noninflamed","P17",[1,5,10],[(10,1)]),("UC inflamed","P07",[2,4,9],[(2,4)])]
for j,(lab,pid,sel,pairs) in enumerate(close):
    ax=fig.add_subplot(subgs[0,j],facecolor="#000000");
    if j==0: ax.text(.02,.98,"E",transform=ax.transAxes,color="white",fontsize=22,fontweight="bold",va="top",zorder=40)
    full=cells[cells.PatientID==pid]; cx,cy,h=crop(full,sel); d=full[full.x.between(cx-h,cx+h)&full.y.between(cy-h,cy+h)]
    oth=~d.neighborhood_id.isin(sel); ax.scatter(d.loc[oth,"x"],d.loc[oth,"y"],s=1.3,c="#343A46",alpha=.75,linewidths=0,rasterized=True)
    for n in sel:
        z=d[d.neighborhood_id==n]; ax.scatter(z.x,z.y,s=3,c=cncols[n],linewidths=0,rasterized=True)
    for a,b in pairs:
        q=nearest(d,a,b)
        if q:
            p1,p2=q; ax.add_patch(FancyArrowPatch(p1,p2,arrowstyle="<->",mutation_scale=17,lw=3,color="white")); m=(p1+p2)/2
            ax.text(m[0],m[1]+25,f"CN{a} ↔ CN{b}",color="black",ha="center",fontsize=9,fontweight="bold",bbox=dict(facecolor="white",edgecolor="black",boxstyle="round,pad=.2"))
    ax.set_xlim(cx-h,cx+h); ax.set_ylim(cy+h,cy-h); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(lab,color="black",fontsize=13)
    add_scalebar(ax,200)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color("#69717D"); sp.set_linewidth(1.7)

fig.savefig(FIG/"Figure_22_CODEX_manuscript_multipanel.png",dpi=300,bbox_inches="tight",facecolor="white")
fig.savefig(FIG/"Figure_22_CODEX_manuscript_multipanel.pdf",bbox_inches="tight",facecolor="white")
plt.close(fig)
print("Saved Figure 22")
