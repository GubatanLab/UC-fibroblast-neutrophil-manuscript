from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu, kruskal
from sklearn.cluster import MiniBatchKMeans, KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns

SEED=20260710
rng=np.random.default_rng(SEED)
GROUPS=["Control","UC_Noninflamed","UC_Inflamed"]
COLORS={"Control":"#4C78A8","UC_Noninflamed":"#E6BE4A","UC_Inflamed":"#D95F5F"}
OUT=Path("giotto_codex_results/additional_analyses"); OUT.mkdir(parents=True,exist_ok=True)
TAB=OUT/"tables"; FIG=Path("giotto_codex_results/figures"); TAB.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
cells=pd.read_csv(OUT/"codex_extended_cells.csv")
nh=pd.read_csv("giotto_codex_results/nolan_neighborhoods/cells_with_nolan_neighborhoods.csv",usecols=["cell_ID","neighborhood","neighborhood_id"])
cells=cells.merge(nh,on="cell_ID",how="left")
cells["cell_type"]=cells.cell_type.fillna("Unclassified").astype(str)
meta=cells[["PatientID","Diagnosis2","Medication"]].drop_duplicates("PatientID")
cell_types=sorted(cells.cell_type.unique()); ct_i={x:i for i,x in enumerate(cell_types)}

def bh(x):
    x=np.asarray(x,float); out=np.full(len(x),np.nan); ok=np.isfinite(x)
    if ok.any(): out[ok]=multipletests(x[ok],method="fdr_bh")[1]
    return out
def pstats(df,value,keys):
    rows=[]
    for key,sub in df.groupby(keys,dropna=False):
        if not isinstance(key,tuple): key=(key,)
        base=dict(zip(keys,key))
        for a,b in [("UC_Noninflamed","Control"),("UC_Inflamed","Control"),("UC_Inflamed","UC_Noninflamed")]:
            va=sub.loc[sub.Diagnosis2==a,value].dropna().values; vb=sub.loc[sub.Diagnosis2==b,value].dropna().values
            if len(va) and len(vb):
                rows.append({**base,"contrast":f"{a} vs {b}","group_1":a,"group_2":b,
                    "n_group_1":len(va),"n_group_2":len(vb),"median_group_1":np.median(va),
                    "median_group_2":np.median(vb),"median_difference":np.median(va)-np.median(vb),
                    "p_value":mannwhitneyu(va,vb,alternative="two-sided").pvalue})
    z=pd.DataFrame(rows)
    if len(z): z["p_adj_BH"]=z.groupby("contrast")["p_value"].transform(bh)
    return z

# Shared spatial indices and local neighbor matrices.
spatial={}; knn10={}; knn30={}; row_by_patient={}
for pid,idx in cells.groupby("PatientID",sort=False).groups.items():
    idx=np.asarray(list(idx),int); row_by_patient[pid]=idx
    xy=cells.loc[idx,["x","y"]].values; tree=cKDTree(xy); spatial[pid]=(idx,xy,tree)
    knn10[pid]=tree.query(xy,k=min(11,len(idx)))[1][:,1:11]
    knn30[pid]=tree.query(xy,k=min(31,len(idx)))[1][:,1:31]

# 1. Independent composition-adjusted KNN colocalization.
pairs=[("Fibroblast","Neutrophil"),("Fibroblast","Neutrophil_PADI4"),("Fibroblast","Neutrophil_CXCR4"),
       ("Fibroblast","Neutrophil_MX1"),("Fibroblast","Neutrophil_OSM"),("Epithelial Cell","Neutrophil"),
       ("Endothelial Cell","Neutrophil"),("Macrophage","Neutrophil"),("Epithelial Cell","CD8 T"),("Epithelial Cell","CD4 T")]
interaction=[]; multiscale=[]; radii=[20,50,100,200]
for pid,(idx,xy,tree) in spatial.items():
    typ=cells.loc[idx,"cell_type"].values; group=cells.loc[idx,"Diagnosis2"].iloc[0]; nbr=knn10[pid]
    for a,b in pairs:
        ia=np.where(typ==a)[0]; ib=np.where(typ==b)[0]
        if len(ia)<5 or len(ib)<5: continue
        obs=(typ[nbr[ia]]==b).mean(); exp=(typ==b).mean()
        interaction.append([pid,group,a,b,len(ia),len(ib),obs,exp,np.log2((obs+1e-5)/(exp+1e-5))])
        tb=cKDTree(xy[ib]); dist=tb.query(xy[ia],k=1)[0]
        for r in radii: multiscale.append([pid,group,a,b,r,np.median(dist),np.mean(dist<=r)])
interaction=pd.DataFrame(interaction,columns=["PatientID","Diagnosis2","celltype_1","celltype_2","n_1","n_2","observed_neighbor_fraction","expected_fraction","log2_enrichment"])
multiscale=pd.DataFrame(multiscale,columns=["PatientID","Diagnosis2","celltype_1","celltype_2","radius_um","median_nearest_distance","fraction_within_radius"])
interaction.to_csv(TAB/"01_independent_knn_colocalization_by_patient.csv",index=False)
pstats(interaction,"log2_enrichment",["celltype_1","celltype_2"]).to_csv(TAB/"01_independent_knn_colocalization_statistics.csv",index=False)
multiscale.to_csv(TAB/"02_multiscale_distance_by_patient.csv",index=False)
pstats(multiscale,"fraction_within_radius",["celltype_1","celltype_2","radius_um"]).to_csv(TAB/"02_multiscale_distance_statistics.csv",index=False)

# Build local cell-composition windows for stability and density-adjusted regions.
type_idx=cells.cell_type.map(ct_i).values
def windows_for(k):
    w=np.zeros((len(cells),len(cell_types)),np.uint8)
    for pid,idx in row_by_patient.items():
        local=knn30[pid][:,:min(k,knn30[pid].shape[1])]
        vals=type_idx[idx[local]]
        for j,g in enumerate(idx): w[g]=np.bincount(vals[j],minlength=len(cell_types))
    return w
w10=windows_for(10)
fit_idx=rng.choice(len(cells),min(60000,len(cells)),replace=False)
basekm=MiniBatchKMeans(10,random_state=SEED,batch_size=4096,n_init=10).fit(w10[fit_idx]); base=basekm.predict(w10)

# 3. Parameter and leave-one-patient-out stability.
stab=[]
for k in [5,10,20,30]:
    w=w10 if k==10 else windows_for(k)
    for nc in [8,10,12,15]:
        km=MiniBatchKMeans(nc,random_state=SEED+k+nc,batch_size=4096,n_init=5).fit(w[fit_idx])
        lab=km.predict(w)
        stab.append(["parameter",k,nc,"all",adjusted_rand_score(base,lab)])
for pid,idx in row_by_patient.items():
    train=np.setdiff1d(fit_idx,idx,assume_unique=False)
    km=MiniBatchKMeans(10,random_state=SEED,batch_size=4096,n_init=5).fit(w10[train])
    lab=km.predict(w10[idx]); stab.append(["leave_one_patient_out",10,10,pid,adjusted_rand_score(base[idx],lab)])
stab=pd.DataFrame(stab,columns=["analysis","k_neighbors","n_clusters","omitted_patient","adjusted_rand_index"])
stab.to_csv(TAB/"03_neighborhood_stability.csv",index=False)

# 4. Patient-aware binomial GLMs for neighborhood frequency.
freq=pd.read_csv("giotto_codex_results/nolan_neighborhoods/neighborhood_frequency_by_patient.csv")
glmrows=[]
for n,sub in freq.groupby("neighborhood"):
    sub=sub.copy(); sub["UC_Noninflamed"]=(sub.Diagnosis2=="UC_Noninflamed").astype(int); sub["UC_Inflamed"]=(sub.Diagnosis2=="UC_Inflamed").astype(int)
    X=sm.add_constant(sub[["UC_Noninflamed","UC_Inflamed"]]); y=sub.frequency
    try:
        # Each patient contributes one independent proportion; do not weight by cell count.
        fit=sm.GLM(y,X,family=sm.families.Binomial()).fit(cov_type="HC3")
        for term in ["UC_Noninflamed","UC_Inflamed"]: glmrows.append([n,term,fit.params[term],fit.bse[term],fit.pvalues[term],np.exp(fit.params[term])])
    except Exception: pass
glm=pd.DataFrame(glmrows,columns=["neighborhood","term","log_odds","robust_se","p_value","odds_ratio"]); glm["p_adj_BH"]=glm.groupby("term").p_value.transform(bh)
glm.to_csv(TAB/"04_neighborhood_binomial_GLM.csv",index=False)

# 5. Distances to epithelium.
targets=["Neutrophil","Neutrophil_PADI4","Neutrophil_CXCR4","Neutrophil_MX1","Neutrophil_OSM","Fibroblast","CD8 T","CD4 T"]
epi=[]
for pid,(idx,xy,tree) in spatial.items():
    typ=cells.loc[idx,"cell_type"].values; ie=np.where(typ=="Epithelial Cell")[0]
    if len(ie)<5: continue
    et=cKDTree(xy[ie]); group=cells.loc[idx,"Diagnosis2"].iloc[0]
    for t in targets:
        it=np.where(typ==t)[0]
        if len(it)<5: continue
        d=et.query(xy[it])[0]; epi.append([pid,group,t,len(it),np.median(d),np.mean(d<=20),np.mean(d<=50),np.mean(d<=100)])
epi=pd.DataFrame(epi,columns=["PatientID","Diagnosis2","cell_type","n_cells","median_distance_to_epithelium","fraction_20um","fraction_50um","fraction_100um"])
epi.to_csv(TAB/"05_epithelial_distance_by_patient.csv",index=False); pstats(epi,"median_distance_to_epithelium",["cell_type"]).to_csv(TAB/"05_epithelial_distance_statistics.csv",index=False)

# 6. Neutrophil subtype local composition (20 nearest neighbors).
subtypes=["Neutrophil_PADI4","Neutrophil_CXCR4","Neutrophil_MX1","Neutrophil_OSM","Neutrophil"]
niche=[]
for pid,idx in row_by_patient.items():
    typ=cells.loc[idx,"cell_type"].values; group=cells.loc[idx,"Diagnosis2"].iloc[0]; nbr=knn30[pid][:,:20]
    for s in subtypes:
        ii=np.where(typ==s)[0]
        if len(ii)<5: continue
        vals=typ[nbr[ii]].ravel()
        for ct in cell_types: niche.append([pid,group,s,ct,len(ii),np.mean(vals==ct)])
niche=pd.DataFrame(niche,columns=["PatientID","Diagnosis2","neutrophil_subtype","neighbor_cell_type","n_index_cells","neighbor_fraction"])
niche.to_csv(TAB/"06_neutrophil_subtype_niche_by_patient.csv",index=False); pstats(niche,"neighbor_fraction",["neutrophil_subtype","neighbor_cell_type"]).to_csv(TAB/"06_neutrophil_subtype_niche_statistics.csv",index=False)

# 7. Cell-type-specific marker states at patient level.
marker_cols=[c for c in cells.columns if c in ["PADI4","CXCR4","MX.1","OSM","CD66b","CD15","CD16","FAPa","CD140a","aSMA","Vimentin","Collagen.4","Podoplanin","CD34","CD68","CD163","CD14","HLA.DR","CD11b","CD4","CD8","FoxP3","CD25","CD279","TIGIT","LAG.3","EpCAM","CK7","P53","Ki67","HLA.ABC"]]
state=cells.groupby(["PatientID","Diagnosis2","cell_type"])[marker_cols].median().reset_index()
state.to_csv(TAB/"07_marker_state_patient_medians.csv",index=False)
long=state.melt(id_vars=["PatientID","Diagnosis2","cell_type"],var_name="marker",value_name="median_intensity").dropna()
pstats(long,"median_intensity",["cell_type","marker"]).to_csv(TAB/"07_marker_state_statistics.csv",index=False)

# 8. Fibroblast state clustering and abundance.
fib=cells[cells.cell_type=="Fibroblast"].copy(); fib_mark=[c for c in ["FAPa","CD140a","aSMA","Vimentin","Collagen.4","Podoplanin","CD34"] if c in fib]
X=np.log1p(fib[fib_mark].clip(lower=0).fillna(0).values); X=StandardScaler().fit_transform(X)
fkm=KMeans(4,random_state=SEED,n_init=30).fit(X); fib["fibroblast_cluster"]=fkm.labels_+1
cent=pd.DataFrame(fkm.cluster_centers_,columns=fib_mark,index=range(1,5)); fn={i:f"F{i}: "+"/".join(cent.loc[i].nlargest(2).index) for i in cent.index}
fib["fibroblast_state"]=fib.fibroblast_cluster.map(fn); fib[["cell_ID","PatientID","Diagnosis2","fibroblast_state"]].to_csv(TAB/"08_fibroblast_state_assignments.csv",index=False)
fc=fib.groupby(["PatientID","Diagnosis2","fibroblast_state"]).size().rename("count").reset_index(); ft=fib.groupby(["PatientID","Diagnosis2"]).size().rename("total").reset_index(); fc=fc.merge(ft); fc["frequency"]=fc["count"]/fc["total"]
fc.to_csv(TAB/"08_fibroblast_state_frequency_by_patient.csv",index=False); pstats(fc,"frequency",["fibroblast_state"]).to_csv(TAB/"08_fibroblast_state_statistics.csv",index=False); cent.to_csv(TAB/"08_fibroblast_state_marker_zscores.csv")

# 9. Neighborhood boundary mixing on KNN edges.
mix=[]
for pid,idx in row_by_patient.items():
    lab=cells.loc[idx,"neighborhood_id"].astype(int).values; group=cells.loc[idx,"Diagnosis2"].iloc[0]; nbr=knn10[pid]
    for a in range(1,11):
        ia=np.where(lab==a)[0]
        if not len(ia): continue
        for b in range(a,11): mix.append([pid,group,a,b,np.mean(lab[nbr[ia]]==b)])
mix=pd.DataFrame(mix,columns=["PatientID","Diagnosis2","neighborhood_1","neighborhood_2","boundary_fraction"])
mix.to_csv(TAB/"09_neighborhood_mixing_by_patient.csv",index=False); pstats(mix,"boundary_fraction",["neighborhood_1","neighborhood_2"]).to_csv(TAB/"09_neighborhood_mixing_statistics.csv",index=False)

# 10. Density-residualized local tissue regions.
w30=windows_for(30).astype(float)/30
density=np.zeros(len(cells))
for pid,(idx,xy,tree) in spatial.items(): density[idx]=tree.query(xy,k=min(31,len(idx)))[0][:,-1]
ld=np.log1p(density); residual=np.empty_like(w30)
for j in range(w30.shape[1]): residual[:,j]=w30[:,j]-np.polyval(np.polyfit(ld,w30[:,j],2),ld)
rfit=rng.choice(len(cells),min(60000,len(cells)),replace=False); rkm=MiniBatchKMeans(8,random_state=SEED,n_init=15,batch_size=4096).fit(residual[rfit]); rl=rkm.predict(residual)
top=np.argsort(rkm.cluster_centers_,axis=1)[:,-2:][:,::-1]; rnames={i:f"DR{i+1}: {cell_types[top[i,0]]}/{cell_types[top[i,1]]}" for i in range(8)}
cells["density_adjusted_region"]=[rnames[i] for i in rl]
cells[["cell_ID","PatientID","Diagnosis2","density_adjusted_region"]].to_csv(TAB/"10_density_adjusted_region_assignments.csv",index=False)
rc=cells.groupby(["PatientID","Diagnosis2","density_adjusted_region"]).size().rename("count").reset_index(); rt=cells.groupby(["PatientID","Diagnosis2"]).size().rename("total").reset_index(); rc=rc.merge(rt); rc["frequency"]=rc["count"]/rc.total
rc.to_csv(TAB/"10_density_adjusted_region_frequency.csv",index=False); pstats(rc,"frequency",["density_adjusted_region"]).to_csv(TAB/"10_density_adjusted_region_statistics.csv",index=False)

# Compact publication figures.
sns.set_theme(style="whitegrid",context="talk")
plt.rcParams.update({"lines.linewidth":2.6,"lines.markersize":8,"axes.linewidth":1.4,
                     "font.weight":"medium","axes.titleweight":"bold"})
display={"Control":"Control","UC_Noninflamed":"UC noninflamed","UC_Inflamed":"UC inflamed"}
q=multiscale[multiscale.celltype_1=="Fibroblast"].groupby(["Diagnosis2","celltype_2","radius_um"]).fraction_within_radius.median().reset_index(); q["Group"]=q.Diagnosis2.map(display)
fig,axs=plt.subplots(2,3,figsize=(15,9),sharex=True,sharey=True)
for ax,(ct,sub) in zip(axs.flat,q.groupby("celltype_2",sort=False)):
    sns.lineplot(sub,x="radius_um",y="fraction_within_radius",hue="Group",palette={display[k]:v for k,v in COLORS.items()},marker="o",linewidth=3,markersize=9,ax=ax)
    ax.set_title(ct.replace("Neutrophil_","")); ax.set_xlabel("Radius (µm)"); ax.set_ylabel("Median patient fraction")
    if ax is not axs.flat[0]: ax.get_legend().remove()
for ax in axs.flat[len(q.celltype_2.unique()):]: ax.axis("off")
fig.suptitle("Fibroblast proximity to neutrophil states across spatial scales",fontweight="bold"); fig.tight_layout(); fig.savefig(FIG/"Figure_13_multiscale_spatial_proximity.png",dpi=300,bbox_inches="tight"); plt.close(fig)
ep=epi.copy(); ep["Group"]=ep.Diagnosis2.map(display)
fig,ax=plt.subplots(figsize=(12,7)); sns.boxplot(ep,x="cell_type",y="median_distance_to_epithelium",hue="Group",palette={display[k]:v for k,v in COLORS.items()},linewidth=1.8,ax=ax); sns.stripplot(ep,x="cell_type",y="median_distance_to_epithelium",hue="Group",dodge=True,color="black",size=6,linewidth=.4,edgecolor="white",legend=False,ax=ax); ax.tick_params(axis="x",rotation=40); ax.set(title="Distance to epithelial compartment",xlabel="",ylabel="Median distance (µm)"); fig.tight_layout(); fig.savefig(FIG/"Figure_14_epithelial_distance.png",dpi=300,bbox_inches="tight"); plt.close(fig)
heat=niche.groupby(["neutrophil_subtype","neighbor_cell_type"]).neighbor_fraction.mean().unstack(fill_value=0); fig,ax=plt.subplots(figsize=(14,5)); sns.heatmap(heat,cmap="mako",ax=ax); ax.set(title="Neutrophil subtype local niche composition",xlabel="Neighbor cell type",ylabel="Neutrophil subtype"); fig.tight_layout(); fig.savefig(FIG/"Figure_15_neutrophil_subtype_niches.png",dpi=300,bbox_inches="tight"); plt.close(fig)
fc["Group"]=fc.Diagnosis2.map(display); fig,axs=plt.subplots(1,2,figsize=(16,6)); sns.heatmap(cent,cmap="vlag",center=0,ax=axs[0]); axs[0].set_title("Fibroblast-state marker z scores"); sns.boxplot(fc,x="fibroblast_state",y="frequency",hue="Group",palette={display[k]:v for k,v in COLORS.items()},ax=axs[1]); axs[1].tick_params(axis="x",rotation=35); axs[1].set_title("Fibroblast-state abundance"); fig.tight_layout(); fig.savefig(FIG/"Figure_16_fibroblast_states.png",dpi=300,bbox_inches="tight"); plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5)); sns.boxplot(stab,x="analysis",y="adjusted_rand_index",linewidth=1.8,ax=ax); sns.stripplot(stab,x="analysis",y="adjusted_rand_index",color="black",size=6,linewidth=.4,edgecolor="white",ax=ax); ax.set(title="Cellular-neighborhood robustness",xlabel="",ylabel="Adjusted Rand index"); fig.tight_layout(); fig.savefig(FIG/"Figure_17_neighborhood_stability.png",dpi=300,bbox_inches="tight"); plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(16,5))
for ax,g in zip(axs,GROUPS):
    z=mix[mix.Diagnosis2==g].groupby(["neighborhood_1","neighborhood_2"]).boundary_fraction.median().unstack(fill_value=0)
    sns.heatmap(z,cmap="magma",vmin=0,vmax=.25,ax=ax,cbar=ax is axs[-1]); ax.set_title(display[g]); ax.set_xlabel("Neighbor CN"); ax.set_ylabel("Index CN")
fig.suptitle("Neighborhood boundary mixing",fontweight="bold"); fig.tight_layout(); fig.savefig(FIG/"Figure_18_neighborhood_mixing.png",dpi=300,bbox_inches="tight"); plt.close(fig)
rc["Group"]=rc.Diagnosis2.map(display); fig,ax=plt.subplots(figsize=(14,6)); sns.boxplot(rc,x="density_adjusted_region",y="frequency",hue="Group",palette={display[k]:v for k,v in COLORS.items()},ax=ax); ax.tick_params(axis="x",rotation=35); ax.set(title="Density-adjusted tissue regions",xlabel="",ylabel="Patient-level frequency"); fig.tight_layout(); fig.savefig(FIG/"Figure_19_density_adjusted_regions.png",dpi=300,bbox_inches="tight"); plt.close(fig)

with open(OUT/"analysis_manifest.txt","w") as f:
    f.write("Analyses 1-10 completed\nInput cells: %d\nPatients: %d\nSeed: %d\n"%(len(cells),cells.PatientID.nunique(),SEED))
print("Completed analyses 1-10 for",len(cells),"cells")
