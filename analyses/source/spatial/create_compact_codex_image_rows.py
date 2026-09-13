from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
from scipy.spatial import cKDTree

ROOT=Path("giotto_codex_results"); FIG=ROOT/"figures"
cells=pd.read_csv(ROOT/"nolan_neighborhoods/cells_with_nolan_neighborhoods.csv")
freq=pd.read_csv(ROOT/"nolan_neighborhoods/neighborhood_frequency_by_patient.csv")
stats=pd.read_csv(ROOT/"nolan_neighborhoods/statistics_neighborhood_pairwise_wilcoxon.csv")
cncols={1:"#00C8FF",2:"#FF7A00",3:"#35E65C",4:"#FF3154",5:"#A66BFF",6:"#D98C5F",7:"#FF5EC4",8:"#9AA3AE",9:"#D9C900",10:"#00B9A8"}
groups=["Control","UC_Noninflamed","UC_Inflamed"]
glabels={"Control":"Control","UC_Noninflamed":"UC noninflamed","UC_Inflamed":"UC inflamed"}
gcols={"Control":"#4C78A8","UC_Noninflamed":"#D4AD35","UC_Inflamed":"#D95F5F"}

# Compact, self-contained neighborhood-abundance panel C.
sns.set_theme(style="whitegrid",context="talk")
fig,axs=plt.subplots(1,10,figsize=(26,4.2),facecolor="white",gridspec_kw={"wspace":.58})
focus=list(range(1,11))
for j,cn in enumerate(focus):
    ax=axs[j]; q=freq[freq.neighborhood.str.startswith(f"CN{cn}:")].copy(); q["Group"]=q.Diagnosis2.map(glabels)
    order=[glabels[g] for g in groups]; pal={glabels[k]:v for k,v in gcols.items()}
    sns.boxplot(q,x="Group",y="frequency",order=order,palette=pal,width=.62,linewidth=1.5,showfliers=False,ax=ax)
    sns.stripplot(q,x="Group",y="frequency",order=order,color="black",size=5,linewidth=.35,edgecolor="white",ax=ax)
    ax.set_title(f"CN{cn}",fontsize=12,fontweight="bold"); ax.set_xlabel(""); ax.set_ylabel("")
    ax.set_xticklabels(["Control","Noninfl.","Inflamed"],rotation=55,ha="right",fontsize=7)
    ax.tick_params(axis="y",labelsize=7,pad=1)
    sig=stats[(stats.neighborhood.str.startswith(f"CN{cn}:"))&(stats.p_adj_BH<.05)].sort_values("p_adj_BH")
    ymax=q.frequency.max(); step=max(ymax*.10,.012); base=ymax+step*.2; xpos={g:i for i,g in enumerate(groups)}
    for k,r in enumerate(sig.itertuples()):
        x1,x2=xpos[r.group_1],xpos[r.group_2]; y=base+k*step; hh=step*.18
        ax.plot([x1,x1,x2,x2],[y,y+hh,y+hh,y],c="black",lw=1.4)
        star="***" if r.p_adj_BH<.001 else "**" if r.p_adj_BH<.01 else "*"; ax.text((x1+x2)/2,y+hh,star,ha="center",va="bottom",fontweight="bold")
    if len(sig): ax.set_ylim(top=base+(len(sig)+.45)*step)
fig.subplots_adjust(left=.065,right=.995,bottom=.22,top=.84)
fig.text(.018,.53,"Patient-level frequency",rotation=90,ha="center",va="center",fontsize=13)
fig.text(.010,.94,"C",fontsize=23,fontweight="bold",va="top")
fig.savefig(FIG/"Figure_22_panel_C_compact.png",dpi=300,bbox_inches="tight",facecolor="white"); plt.close(fig)

def scalebar(ax,length):
    span=abs(np.diff(ax.get_xlim())[0]); frac=min(length/span,.34); x,y=.065,.07
    ax.plot([x,x+frac],[y,y],transform=ax.transAxes,color="white",lw=4,solid_capstyle="butt",zorder=30)
    ax.text(x+frac/2,y+.025,f"{length} µm",transform=ax.transAxes,color="white",ha="center",va="bottom",fontsize=10,fontweight="bold",zorder=30)

def style(ax,title):
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title,fontsize=15,fontweight="bold",pad=7)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_color("#69717D"); sp.set_linewidth(1.7)

maps=[("Control","P03",[3,6,7,8]),("UC noninflamed","P17",[1,5,10]),("UC inflamed","P07",[2,4,9])]
fig,axs=plt.subplots(1,3,figsize=(14.8,4.75),facecolor="white",gridspec_kw={"wspace":.035})
for j,(group,pid,sel) in enumerate(maps):
    ax=axs[j]; ax.set_facecolor("black"); z=cells[cells.PatientID==pid]; other=~z.neighborhood_id.isin(sel)
    ax.scatter(z.loc[other,"x"],z.loc[other,"y"],s=.75,c="#343A46",alpha=.72,linewidths=0,rasterized=True)
    for n in sel:
        q=z[z.neighborhood_id==n]; ax.scatter(q.x,q.y,s=1.6,c=cncols[n],alpha=.98,linewidths=0,rasterized=True)
    xmin,xmax=z.x.min(),z.x.max(); ymin,ymax=z.y.min(),z.y.max(); span=max(xmax-xmin,ymax-ymin)
    xm,ym=(xmin+xmax)/2,(ymin+ymax)/2; ax.set_xlim(xm-span/2,xm+span/2); ax.set_ylim(ym-span/2,ym+span/2)
    ax.invert_yaxis(); style(ax,f"{group} ({pid})"); scalebar(ax,500)
    handles=[Line2D([0],[0],marker="o",linestyle="",markerfacecolor=cncols[n],markeredgecolor="none",markersize=7,label=f"CN{n}") for n in sel]
    ax.legend(handles=handles,loc="lower right",ncol=len(sel),frameon=False,labelcolor="white",fontsize=8,columnspacing=.8,handletextpad=.2)
axs[0].text(.015,.98,"D",transform=axs[0].transAxes,color="white",fontsize=23,fontweight="bold",va="top",zorder=40)
fig.savefig(FIG/"Figure_22_panel_D_compact_CODEX.png",dpi=300,bbox_inches="tight",facecolor="white"); plt.close(fig)

def crop_center(z,sel,h=420):
    q=z[z.neighborhood_id.isin(sel)]; xy=q[["x","y"]].to_numpy(); sample=xy[np.linspace(0,len(xy)-1,min(600,len(xy)),dtype=int)]
    tree=cKDTree(xy); p=sample[np.argmax([len(tree.query_ball_point(v,h)) for v in sample])]; return p[0],p[1],h

def nearest(z,a,b):
    x=z[z.neighborhood_id==a]; y=z[z.neighborhood_id==b]
    if not len(x) or not len(y): return None
    tree=cKDTree(y[["x","y"]]); dist,ix=tree.query(x[["x","y"]]); i=np.argmin(dist)
    return x.iloc[i][["x","y"]].to_numpy(float),y.iloc[ix[i]][["x","y"]].to_numpy(float)

close=[("Control","P03",[3,6,7,8],(7,8)),("UC noninflamed","P17",[1,5,10],(10,1)),("UC inflamed","P07",[2,4,9],(2,4))]
fig,axs=plt.subplots(1,3,figsize=(14.8,4.75),facecolor="white",gridspec_kw={"wspace":.035})
for j,(group,pid,sel,pair) in enumerate(close):
    ax=axs[j]; ax.set_facecolor("black"); full=cells[cells.PatientID==pid]; cx,cy,h=crop_center(full,sel)
    z=full[full.x.between(cx-h,cx+h)&full.y.between(cy-h,cy+h)]; other=~z.neighborhood_id.isin(sel)
    ax.scatter(z.loc[other,"x"],z.loc[other,"y"],s=1.4,c="#343A46",alpha=.74,linewidths=0,rasterized=True)
    for n in sel:
        q=z[z.neighborhood_id==n]; ax.scatter(q.x,q.y,s=3.2,c=cncols[n],linewidths=0,rasterized=True)
    link=nearest(z,*pair)
    if link:
        a,b=link; ax.add_patch(FancyArrowPatch(a,b,arrowstyle="<->",mutation_scale=18,lw=3,color="white")); m=(a+b)/2
        ax.text(m[0],m[1]+25,f"CN{pair[0]} ↔ CN{pair[1]}",color="white",ha="center",fontsize=9,fontweight="bold",bbox=dict(facecolor="black",edgecolor="white",boxstyle="round,pad=.2"))
    ax.set_xlim(cx-h,cx+h); ax.set_ylim(cy+h,cy-h); style(ax,group); scalebar(ax,200)
axs[0].text(.015,.98,"E",transform=axs[0].transAxes,color="white",fontsize=23,fontweight="bold",va="top",zorder=40)
fig.savefig(FIG/"Figure_22_panel_E_compact_CODEX.png",dpi=300,bbox_inches="tight",facecolor="white"); plt.close(fig)
print("Saved compact CODEX panels D and E")
