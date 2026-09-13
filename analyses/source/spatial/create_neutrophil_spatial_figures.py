from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore")
ROOT=Path("giotto_codex_results")
OUT=ROOT/"neutrophil_spatial_analyses"
FIG=OUT/"figures"; TAB=OUT/"tables"
FIG.mkdir(parents=True,exist_ok=True); TAB.mkdir(parents=True,exist_ok=True)

d=pd.read_csv(ROOT/"additional_analyses/codex_extended_cells.csv")
cn=pd.read_csv(ROOT/"nolan_neighborhoods/cells_with_nolan_neighborhoods.csv",usecols=["cell_ID","neighborhood_id"])
fs=pd.read_csv(ROOT/"additional_analyses/tables/08_fibroblast_state_assignments.csv",usecols=["cell_ID","fibroblast_state"])
d=d.merge(cn,on="cell_ID",how="left").merge(fs,on="cell_ID",how="left")
d["subtype"]=d["Neutrophil_subtype_0.8"].map(lambda x:f"Neutrophil_{x}" if pd.notna(x) else np.nan)
groups=["Control","UC_Noninflamed","UC_Inflamed"]
glabel={"Control":"Control","UC_Noninflamed":"UC noninflamed","UC_Inflamed":"UC inflamed"}
gcols={"Control":"#4C78A8","UC_Noninflamed":"#D4AD35","UC_Inflamed":"#D95F5F"}
subs=["Neutrophil_PADI4","Neutrophil_CXCR4","Neutrophil_MX1","Neutrophil_OSM"]
scols=dict(zip(subs,["#FF7A00","#FF3154","#D9C900","#00B9A8"]))
slabel={x:x.replace("Neutrophil_","") for x in subs}
sns.set_theme(style="whitegrid",context="talk")
plt.rcParams.update({"font.family":"Arial","axes.titleweight":"bold","axes.linewidth":1.3,"lines.linewidth":2.4})

def save(fig,name):
    fig.savefig(FIG/f"{name}.png",dpi=300,bbox_inches="tight",facecolor="white")
    fig.savefig(FIG/f"{name}.pdf",bbox_inches="tight",facecolor="white")
    plt.close(fig)

def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); q[o]=np.minimum.accumulate((p[o]*n/np.arange(1,n+1))[::-1])[::-1]
    return np.minimum(q,1)

def patient_trees(pid):
    z=d[d.PatientID==pid]
    trees={}
    for key,mask in {
        "Fibroblast":z.cell_type.eq("Fibroblast"),"Epithelial":z.cell_type.eq("Epithelial Cell"),
        "Endothelial":z.cell_type.eq("Endothelial Cell"),"Neutrophil":z.subtype.notna()}.items():
        xy=z.loc[mask,["x","y"]].to_numpy(); trees[key]=cKDTree(xy) if len(xy) else None
    return z,trees

# 1. Neutrophil subtype proximity to fibroblasts
prox=[]
cell_dist=[]
for pid in d.PatientID.unique():
    z,tr=patient_trees(pid)
    if tr["Fibroblast"] is None: continue
    for s in subs:
        n=z[z.subtype==s];
        if not len(n): continue
        dist=tr["Fibroblast"].query(n[["x","y"]])[0]
        prox.append([pid,n.Diagnosis2.iloc[0],s,len(n),np.median(dist),np.mean(dist<=50),np.mean(dist<=30)])
        cell_dist.extend(zip(n.cell_ID,[pid]*len(n),[s]*len(n),dist))
prox=pd.DataFrame(prox,columns=["PatientID","Diagnosis2","subtype","n_neutrophils","median_nearest_fibroblast_um","fraction_within_50um","fraction_within_30um"])
prox.to_csv(TAB/"01_neutrophil_fibroblast_proximity_by_patient.csv",index=False)
pd.DataFrame(cell_dist,columns=["cell_ID","PatientID","subtype","nearest_fibroblast_um"]).to_csv(TAB/"01_neutrophil_fibroblast_distances_by_cell.csv",index=False)
fig,axs=plt.subplots(1,2,figsize=(16,6))
for ax,y,title in zip(axs,["median_nearest_fibroblast_um","fraction_within_50um"],["Nearest fibroblast distance","Neutrophils within 50 µm of fibroblasts"]):
    q=prox.copy(); q["Group"]=q.Diagnosis2.map(glabel); q["Subtype"]=q.subtype.map(slabel)
    sns.boxplot(q,x="Subtype",y=y,hue="Group",hue_order=[glabel[g] for g in groups],palette={glabel[k]:v for k,v in gcols.items()},showfliers=False,ax=ax)
    sns.stripplot(q,x="Subtype",y=y,hue="Group",hue_order=[glabel[g] for g in groups],dodge=True,color="black",size=3.5,ax=ax,legend=False)
    ax.set_title(title); ax.set_xlabel(""); ax.set_ylabel("Distance (µm)" if "distance" in title.lower() else "Patient-level fraction")
axs[0].legend_.remove(); axs[1].legend(title="",frameon=False,bbox_to_anchor=(1.02,1),loc="upper left")
fig.suptitle("Neutrophil subtype proximity to fibroblasts",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_N01_neutrophil_fibroblast_proximity")

# 2. Fibroblast-state enrichment among nearby fibroblasts
fib=d[d.cell_type.eq("Fibroblast") & d.fibroblast_state.notna()]
states=sorted(fib.fibroblast_state.unique()); enr=[]
for pid in d.PatientID.unique():
    z=d[d.PatientID==pid]; f=fib[fib.PatientID==pid]
    if not len(f): continue
    tree=cKDTree(f[["x","y"]]); base=f.fibroblast_state.value_counts(normalize=True)
    for s in subs:
        n=z[z.subtype==s]
        if not len(n): continue
        dist,ix=tree.query(n[["x","y"]]); near=f.iloc[ix[dist<=50]].fibroblast_state
        for st in states:
            obs=(near==st).sum(); exp=max(len(near)*base.get(st,0),.5)
            enr.append([pid,n.Diagnosis2.iloc[0],s,st,len(near),np.log2((obs+.5)/(exp+.5))])
enr=pd.DataFrame(enr,columns=["PatientID","Diagnosis2","subtype","fibroblast_state","n_within_50um","log2_observed_expected"])
enr.to_csv(TAB/"02_neutrophil_fibroblast_state_enrichment_by_patient.csv",index=False)
fig,axs=plt.subplots(1,3,figsize=(20,5.5),sharey=True)
for ax,g in zip(axs,groups):
    m=enr[enr.Diagnosis2==g].groupby(["subtype","fibroblast_state"]).log2_observed_expected.mean().unstack().reindex(index=subs,columns=states)
    sns.heatmap(m,cmap="vlag",center=0,vmin=-2,vmax=2,ax=ax,cbar=ax is axs[-1],cbar_kws={"label":"Mean log2 observed/expected"})
    ax.set_title(glabel[g]); ax.set_xlabel("Fibroblast state"); ax.set_ylabel("Neutrophil subtype" if ax is axs[0] else "")
    ax.set_yticklabels([slabel.get(x.get_text(),x.get_text()) for x in ax.get_yticklabels()],rotation=0); ax.tick_params(axis="x",rotation=45,labelsize=8)
fig.suptitle("Neutrophil subtype enrichment around fibroblast states (≤50 µm)",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_N02_neutrophil_fibroblast_state_enrichment")

# 3. Distance-dependent fibroblast association curves
bins=np.array([0,20,40,60,100,150,250]); curves=[]
cd=pd.DataFrame(cell_dist,columns=["cell_ID","PatientID","subtype","distance"])
cd=cd.merge(d[["PatientID","Diagnosis2"]].drop_duplicates(),on="PatientID")
for (pid,g,s),q in cd.groupby(["PatientID","Diagnosis2","subtype"]):
    for r in bins[1:]: curves.append([pid,g,s,r,np.mean(q.distance<=r)])
curves=pd.DataFrame(curves,columns=["PatientID","Diagnosis2","subtype","radius_um","cumulative_fraction"]); curves.to_csv(TAB/"03_distance_dependent_fibroblast_association.csv",index=False)
fig,axs=plt.subplots(1,4,figsize=(20,4.8),sharex=True,sharey=True)
for ax,s in zip(axs,subs):
    sns.lineplot(curves[curves.subtype==s],x="radius_um",y="cumulative_fraction",hue="Diagnosis2",hue_order=groups,palette=gcols,estimator="mean",errorbar="se",marker="o",ax=ax)
    ax.set_title(slabel[s]); ax.set_xlabel("Radius from fibroblast (µm)"); ax.set_ylabel("Cumulative fraction" if ax is axs[0] else "")
    if ax.legend_ is not None: ax.legend_.remove()
handles=[Line2D([0],[0],color=gcols[g],marker="o",label=glabel[g]) for g in groups]
axs[-1].legend(handles=handles,title="",frameon=False,bbox_to_anchor=(1.02,1),loc="upper left")
fig.suptitle("Distance-dependent neutrophil–fibroblast association",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_N03_distance_dependent_fibroblast_association")

# 4. Neutrophil subtype composition within cellular neighborhoods
nd=d[d.subtype.notna() & d.neighborhood_id.notna()].copy(); nd["CN"]="CN"+nd.neighborhood_id.astype(int).astype(str)
cnt=nd.groupby(["PatientID","Diagnosis2","CN","subtype"]).size().rename("n").reset_index()
pcs=nd[["PatientID","Diagnosis2","CN"]].drop_duplicates()
full=pcs.merge(pd.DataFrame({"subtype":subs}),how="cross").merge(cnt,how="left").fillna({"n":0})
full["fraction"]=full.n/full.groupby(["PatientID","CN"]).n.transform("sum"); cnt=full
cnt.to_csv(TAB/"04_neutrophil_subtype_composition_by_patient_cn.csv",index=False)
mean=cnt.groupby(["Diagnosis2","CN","subtype"]).fraction.mean().unstack(fill_value=0)
fig,axs=plt.subplots(1,3,figsize=(19,5.5),sharey=True)
for ax,g in zip(axs,groups):
    m=mean.loc[g].reindex([f"CN{i}" for i in range(1,11)]).fillna(0); bottom=np.zeros(len(m))
    for s in subs: ax.bar(m.index,m[s],bottom=bottom,color=scols[s],label=slabel[s],edgecolor="white",linewidth=.3); bottom+=m[s].to_numpy()
    ax.set_title(glabel[g]); ax.set_xlabel("Cellular neighborhood"); ax.tick_params(axis="x",rotation=45); ax.set_ylabel("Mean patient subtype fraction" if ax is axs[0] else "")
axs[-1].legend(title="Neutrophil subtype",frameon=False,bbox_to_anchor=(1.02,1),loc="upper left")
fig.suptitle("Neutrophil subtype composition across cellular neighborhoods",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_N04_neutrophil_subtypes_by_neighborhood")

# 5. Same-subtype spatial aggregation
agg=[]
for (pid,g,s),q in d[d.subtype.notna()].groupby(["PatientID","Diagnosis2","subtype"]):
    xy=q[["x","y"]].to_numpy()
    if len(xy)<2: continue
    dd=cKDTree(xy).query(xy,k=2)[0][:,1]
    agg.append([pid,g,s,len(xy),np.median(dd),np.mean(dd<=30),np.mean(dd<=50)])
agg=pd.DataFrame(agg,columns=["PatientID","Diagnosis2","subtype","n_cells","median_same_subtype_nn_um","fraction_nn_within_30um","fraction_nn_within_50um"]); agg.to_csv(TAB/"05_neutrophil_subtype_aggregation_by_patient.csv",index=False)
fig,axs=plt.subplots(1,2,figsize=(16,6))
for ax,y,title in zip(axs,["median_same_subtype_nn_um","fraction_nn_within_30um"],["Same-subtype nearest-neighbor distance","Fraction in ≤30 µm same-subtype aggregates"]):
    q=agg.copy(); q["Group"]=q.Diagnosis2.map(glabel); q["Subtype"]=q.subtype.map(slabel)
    sns.boxplot(q,x="Subtype",y=y,hue="Group",hue_order=[glabel[g] for g in groups],palette={glabel[k]:v for k,v in gcols.items()},showfliers=False,ax=ax)
    sns.stripplot(q,x="Subtype",y=y,hue="Group",hue_order=[glabel[g] for g in groups],dodge=True,color="black",size=3.5,legend=False,ax=ax)
    ax.set_title(title); ax.set_xlabel(""); ax.set_ylabel("Distance (µm)" if "distance" in title.lower() else "Patient-level fraction")
axs[0].legend_.remove(); axs[1].legend(title="",frameon=False,bbox_to_anchor=(1.02,1),loc="upper left")
fig.suptitle("Neutrophil subtype spatial aggregation",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_N05_neutrophil_spatial_aggregation")

# 6. Nearest tissue compartment localization
loc=[]
for pid in d.PatientID.unique():
    z,tr=patient_trees(pid)
    for s in subs:
        n=z[z.subtype==s]
        if not len(n): continue
        ds=np.column_stack([tr[k].query(n[["x","y"]])[0] for k in ["Epithelial","Fibroblast","Endothelial"]])
        nearest=np.array(["Epithelial","Fibroblast","Endothelial"])[np.argmin(ds,axis=1)]
        for c in ["Epithelial","Fibroblast","Endothelial"]: loc.append([pid,n.Diagnosis2.iloc[0],s,c,np.mean(nearest==c)])
loc=pd.DataFrame(loc,columns=["PatientID","Diagnosis2","subtype","nearest_compartment","fraction"]); loc.to_csv(TAB/"06_neutrophil_compartment_localization_by_patient.csv",index=False)
fig,axs=plt.subplots(1,3,figsize=(18,5.5),sharey=True)
cc={"Epithelial":"#8DD3C7","Fibroblast":"#BC80BD","Endothelial":"#80B1D3"}
for ax,g in zip(axs,groups):
    m=loc[loc.Diagnosis2==g].groupby(["subtype","nearest_compartment"]).fraction.mean().unstack().reindex(subs).fillna(0); bottom=np.zeros(len(m))
    for c in cc: ax.bar([slabel[x] for x in m.index],m[c],bottom=bottom,color=cc[c],label=c,edgecolor="white",linewidth=.4); bottom+=m[c].to_numpy()
    ax.set_title(glabel[g]); ax.set_xlabel(""); ax.set_ylabel("Mean patient fraction" if ax is axs[0] else ""); ax.tick_params(axis="x",rotation=30)
axs[-1].legend(title="Nearest compartment",frameon=False,bbox_to_anchor=(1.02,1),loc="upper left")
fig.suptitle("Neutrophil subtype localization across tissue compartments",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_N06_neutrophil_compartment_localization")

# 7. Endothelial-to-epithelial spatial axis
axis=[]
for pid in d.PatientID.unique():
    z,tr=patient_trees(pid)
    if tr["Endothelial"] is None or tr["Epithelial"] is None: continue
    for s in subs:
        n=z[z.subtype==s]
        if not len(n): continue
        de=tr["Endothelial"].query(n[["x","y"]])[0]; dp=tr["Epithelial"].query(n[["x","y"]])[0]; score=de/(de+dp+1e-6)
        for v in score: axis.append([pid,n.Diagnosis2.iloc[0],s,v])
axis=pd.DataFrame(axis,columns=["PatientID","Diagnosis2","subtype","vascular_to_epithelial_axis"]); axis.to_csv(TAB/"07_neutrophil_vascular_epithelial_axis_by_cell.csv",index=False)
fig,axs=plt.subplots(1,3,figsize=(18,5.5),sharex=True,sharey=True)
for ax,g in zip(axs,groups):
    for s in subs:
        q=axis[(axis.Diagnosis2==g)&(axis.subtype==s)]
        if len(q)>5: sns.kdeplot(q,x="vascular_to_epithelial_axis",color=scols[s],label=slabel[s],bw_adjust=1.2,ax=ax)
    ax.set_title(glabel[g]); ax.set_xlim(0,1); ax.set_xlabel("Endothelial → epithelial spatial axis"); ax.set_ylabel("Density" if ax is axs[0] else "")
axs[-1].legend(title="Neutrophil subtype",frameon=False,bbox_to_anchor=(1.02,1),loc="upper left")
fig.suptitle("Neutrophil states along the vascular-to-epithelial axis",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_N07_vascular_to_epithelial_gradient")

# 8. Marker states in fibroblast-proximal versus fibroblast-distant neutrophils
markers=["PADI4","CXCR4","MX.1","OSM","CD66b","CD15","CD16","HLA.DR","CD11b","Ki67","HLA.ABC"]
distmap=pd.DataFrame(cell_dist,columns=["cell_ID","PatientID","subtype","distance"])
n=d[d.subtype.notna()].merge(distmap[["cell_ID","distance"]],on="cell_ID",how="inner"); mr=[]
for (pid,g,s),q in n.groupby(["PatientID","Diagnosis2","subtype"]):
    near=q[q.distance<=50]; far=q[q.distance>100]
    if len(near)<5 or len(far)<5: continue
    for m in markers:
        mr.append([pid,g,s,m,len(near),len(far),np.log2((near[m].mean()+.1)/(far[m].mean()+.1))])
mr=pd.DataFrame(mr,columns=["PatientID","Diagnosis2","subtype","marker","n_proximal","n_distant","log2_proximal_distant"]); mr.to_csv(TAB/"08_marker_expression_fibroblast_proximal_vs_distant.csv",index=False)
fig,axs=plt.subplots(1,3,figsize=(20,5.5),sharey=True)
for ax,g in zip(axs,groups):
    m=mr[mr.Diagnosis2==g].groupby(["subtype","marker"]).log2_proximal_distant.mean().unstack().reindex(index=subs,columns=markers)
    sns.heatmap(m,cmap="vlag",center=0,vmin=-1.5,vmax=1.5,ax=ax,cbar=ax is axs[-1],cbar_kws={"label":"Mean log2 proximal/distant"})
    ax.set_title(glabel[g]); ax.set_xlabel("Protein marker"); ax.set_ylabel("Neutrophil subtype" if ax is axs[0] else "")
    ax.set_yticks(np.arange(len(subs))+.5); ax.set_yticklabels([slabel[s] for s in subs],rotation=0); ax.tick_params(axis="x",rotation=45,labelsize=9)
fig.suptitle("Neutrophil protein states near fibroblasts (≤50 versus >100 µm)",fontweight="bold"); fig.tight_layout(); save(fig,"Figure_N08_fibroblast_proximity_marker_states")

print(f"Saved 8 figure sets and analysis tables to {OUT}")
