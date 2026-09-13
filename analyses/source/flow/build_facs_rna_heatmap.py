from pathlib import Path
import sys,json,itertools,hashlib
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
MAN=Path('.')
SRC=MAN/'output/Nature_Extended_Data_Consolidated_2026-09-03/source_data/updated_analyses'
OUT=ROOT/'outputs/11-FACS-RNA-module-heatmap';OUT.mkdir(exist_ok=True)
mods=['OSM_inflammation','Type_I_interferon','Degranulation','Oxidative_burst','Chemotaxis','Hypoxia_glycolysis']
names=['Inflammatory cytokines','Interferon response','Granule/protease program','NADPH oxidase machinery','Chemotaxis','Hypoxia/glycolysis']
cons=['control_fibroblasts_vs_alone','UC_vs_alone','UC_vs_control_fibroblasts','alpha_alone','alpha_in_UC']
D=pd.read_csv(MAN/'Figure 3/source_data/Figure_3_program_source_data.csv')
print('Contrasts',D.contrast.unique())
assert set(cons)<=set(D.contrast)
R=D[D.program.isin(mods)].copy();assert len(R)==30
R['q_bh_30']=stats.false_discovery_control(R.p.to_numpy(),method='bh')
don=pd.read_csv(SRC/'population_program_donor_differences.csv')
checks=[]
for idx,r in R.iterrows():
    sub=don[(don.program==r.program)&(don.contrast==r.contrast)&(don.population=='Resolved')&(don.version=='full_gene')&(don.min_cells==20)]
    x=sub.delta.to_numpy();assert len(x)==r.n
    assert np.isclose(x.mean(),r.effect,atol=1e-10)
    assert np.isclose(stats.ttest_1samp(x,0).pvalue,r.p,atol=1e-10)
    signs=np.array(list(itertools.product([-1,1],repeat=len(x))))
    R.loc[idx,'p_exact_signflip']=np.mean(np.abs((signs*x).mean(axis=1))>=abs(x.mean())-1e-12)
R['q_signflip_bh_30']=stats.false_discovery_control(R.p_exact_signflip.to_numpy())
coverage=pd.read_csv(SRC/'program_gene_coverage.csv');coverage=coverage[coverage.program.isin(mods)]
coverage.to_csv(OUT/'RNA_module_gene_membership.csv',index=False)
R.to_csv(OUT/'RNA_module_effects_BH30.csv',index=False)
don[(don.program.isin(mods))&(don.population=='Resolved')&(don.min_cells==20)&don.contrast.isin(cons)].to_csv(OUT/'RNA_donor_differences.csv',index=False)
F=pd.read_csv(ROOT/'outputs/10-single-BH-heatmap/BH_effect_source_data.csv')
F.to_csv(OUT/'FACS_effects_BH45.csv',index=False)
markers=['OSM','CXCR4','MPO','MX1','PADI4'];pairs=list(dict.fromkeys(zip(F.comparison,F.reference)))
labels=['Ctrl FB\nvs alone','UC FB\nvs alone','UC FB\nvs Ctrl FB','NAMPTi\nvs untreated','α5β1i\nvs untreated','NAMPTi\nvs untreated','α5β1i\nvs untreated','Dual blockade\nvs untreated','Dual blockade\nvs NAMPTi']
pos=[0,1,2,3.25,4.25,5.5,6.5,7.5,8.5];rna_cols={0:cons[0],1:cons[1],2:cons[2],4:cons[3],6:cons[4]}
plt.rcParams.update({'font.family':'Arial','font.size':7,'svg.fonttype':'none'})
fig=plt.figure(figsize=(210/25.4,133/25.4));cmap=plt.get_cmap('RdBu_r')
def stars(q):return '***' if q<.001 else '**' if q<.01 else '*' if q<.05 else ''
def panel(bounds,rows,data,rna=False):
    ax=fig.add_axes(bounds);limit=.6 if rna else 3;norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit)
    for j,p in enumerate(pos):
        for i,m in enumerate(rows):
            if rna and j not in rna_cols:
                ax.add_patch(plt.Rectangle((p,i),1,1,fc='#eeeeee',ec='white',lw=.6));ax.text(p+.5,i+.5,'—',ha='center',va='center',color='#888',fontsize=7);continue
            r=data[(data.program==m)&(data.contrast==rna_cols[j])].iloc[0] if rna else data[(data.marker==m)&(data.comparison==pairs[j][0])&(data.reference==pairs[j][1])].iloc[0]
            v=r.effect if rna else r.difference_mean_z;q=r.q_bh_30 if rna else r.q_bh_45
            ax.add_patch(plt.Rectangle((p,i),1,1,fc=cmap(norm(v)),ec='white',lw=.6))
            ax.text(p+.5,i+.5,f'{v:+.2f}{stars(q)}',ha='center',va='center',fontsize=6,color='white' if abs(v)/limit>.58 else '#222')
        if not rna:ax.text(p+.5,-.24,labels[j],ha='center',va='bottom',fontsize=6,clip_on=False)
    ax.set(xlim=(0,9.5),ylim=(len(rows),0),xticks=[],yticks=np.arange(len(rows))+.5,yticklabels=names if rna else rows)
    ax.tick_params(length=0,pad=5)
    for s in ax.spines.values():s.set_visible(False)
    cax=fig.add_axes([.914,bounds[1]+.015,.012,bounds[3]-.03]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,extend='both' if rna else 'neither');cb.set_ticks([-limit,0,limit]);cb.ax.tick_params(labelsize=6,length=2)
    cb.set_label('Δ RNA program score' if rna else 'Δ protein z score',fontsize=6)
    return ax
ax=panel([.25,.52,.64,.255],markers,F)
for x,t in [(1.5,'Fibroblast exposure'),(4.25,'Neutrophils alone'),(7.5,'Blockade in UC FB coculture')]:ax.text(x,-1.38,t,ha='center',va='bottom',fontsize=6.5,fontweight='bold',clip_on=False)
panel([.25,.16,.64,.306],mods,R,True)
fig.text(.02,.968,'Neutrophil protein expression and RNA response programs',fontsize=11,fontweight='bold')
fig.text(.02,.78,'FACS proteins',fontweight='bold',fontsize=8)
fig.text(.02,.477,'RNA modules',fontweight='bold',fontsize=8)
fig.text(.25,.477,'Paired donor effects; n = 5–6',fontsize=6.5)
fig.text(.02,.109,'Stars: BH q < 0.05 (*), < 0.01 (**), < 0.001 (***). Separate families: 45 FACS tests; 30 RNA tests.',fontsize=6.7)
fig.text(.02,.075,'Gray: RNA comparison unavailable. RNA and protein use separate scales. RNA scores infer programs, not measured function.',fontsize=6.7)
fig.text(.02,.04,'Ctrl FB, control fibroblasts; UC FB, UC fibroblasts; NAMPTi, FK-866; α5β1i, ATN-161; dual blockade, both inhibitors.',fontsize=6.7)
for ext in ['png','svg','tiff']:
    fig.savefig(OUT/f'FACS_RNA_BH_heatmap.{ext}',**({'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {'dpi':300}))
plt.close(fig)
report='''# Combined FACS and RNA heatmap

Six existing Figure 3 full-gene programs are added beneath the unchanged FACS BH results. These are existing custom source definitions, not newly downloaded Hallmark or Reactome sets. The exact gene membership is provided. Chemotaxis is shown without an adhesion claim; the source granule program is used intact. Oxidative_burst includes MPO as well as NADPH oxidase components.

RNA: resolved neutrophils, >=20 cells per recorded donor-condition arm. Log-normalized genes were standardized using the source equally weighted untreated donor-condition reference, averaged within programs and then within samples. Displayed values are source mean paired differences, independently reproduced from donor contrasts. Source two-sided paired t-test P values were independently reproduced. BH was recalculated across all 30 displayed RNA tests (six modules x five contrasts); this differs from the prior Figure 3 family of 14 programs separately per contrast. RNA and FACS are separate testing families and have separate color scales. RNA colors saturate at +/-0.6 with exact tile values retained. FACS values and BH45 stars are unchanged. Missing RNA comparisons are gray, not zero.

An exact two-sided paired sign-flip sensitivity analysis was also computed from all 2^n donor sign assignments and BH adjusted across 30 tests; see RNA_module_effects_BH30.csv. Primary RNA stars retain the source paired-t framework rather than choosing a test based on significance. Recorded RNA donor matching remains unverified. FACS pairing remains unconfirmed and uses the existing unpaired permutation tests; the two modalities are not assumed matched. Source identity-markers-removed donor contrasts are included where available, but a new exclusion of all five FACS markers was not performed. No new maturation/contamination adjustment was performed.

The top panel retains the prior FACS limitations: MPO is antibody abundance, not enzyme activity; MX1/PADI4 source assignments conflict with protocol annotations; parent gates and compensation retained; most co-culture samples lack viability measurement. RNA programs do not establish migration, ROS generation, granule release, NET formation, NAD depletion, or fibroblast activation. Adding these rows provides mechanistic context, not independent functional validation. Canonical assembled Figure 3 was not overwritten.
'''
(OUT/'Figure_legend_and_methods.md').write_text(report,encoding='utf8')
audit={'rna_tests':len(R),'rna_bh_significant':int((R.q_bh_30<.05).sum()),'rna_signflip_bh_significant':int((R.q_signflip_bh_30<.05).sum()),'facs_tests':len(F),'facs_bh_significant':int((F.q_bh_45<.05).sum()),'source_effects_and_p_reproduced':True,'sources':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [MAN/'Figure 3/source_data/Figure_3_program_source_data.csv',SRC/'population_program_donor_differences.csv',SRC/'program_gene_coverage.csv']}}
(OUT/'verification.json').write_text(json.dumps(audit,indent=2));print(json.dumps(audit,indent=2))
