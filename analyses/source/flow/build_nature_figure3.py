from pathlib import Path
import sys,shutil,json
ROOT=Path(__file__).resolve().parents[1]
old=(ROOT/'analysis/assemble_figure3_revision.py').read_text(encoding='utf8')
old=old.split("base=OUT/'assembly_base.pdf'")[0]
# Explicit physical layout, with text re-typeset at final size rather than scaling a page.
changes={
"OUT=ROOT/'outputs/15-Figure3-revised'":"OUT=ROOT/'outputs/16-Figure3-Nature'",
"W,H=210,282":"W,H=183,170",
"fontsize=10,fontweight='bold');text(x+5,y,title,fontsize=8":"fontsize=7,fontweight='bold');text(x+4,y,title,fontsize=7",
"heading('b','Neutrophil protein distributions',104,3)":"heading('b','Neutrophil protein distributions',77,2)",
"x=110+(k%2)*47;y=10+(k//2)*4":"x=83+(k%2)*49;y=7+(k//2)*3.3",
"axbox(109+i*19.3,25,16.4,25)":"axbox(83+i*19.5,19,16.2,18)",
"fontsize=4.2":"fontsize=5",
"labelsize=4.5":"labelsize=5",
"text(112,54,'Compensated fluorescence (arcsinh); representative samples',fontsize=5.7)":"text(92,41,'Compensated fluorescence (arcsinh)',fontsize=5)",
"heading('c','Fibroblast exposure and blockade effects on neutrophil proteins',3,62)":"heading('c','Fibroblast exposure and blockade effects on neutrophil proteins',2,47)",
"axbox(25,80,164,29)":"axbox(16,63,151,21)",
"axbox(193,81,2,27)":"axbox(171,63,1.7,21)",
"text(25,112,'BH across 45 FACS tests; n = 6/arm, except Ctrl FB n = 5',fontsize=6)":"",
"heading('d','NETosis Assay (Elastase Release)',3,122)":"heading('d','NETosis Assay (Elastase Release)',2,89)",
"axbox(18,132,57,50)":"axbox(14,98,45,25)",
"text(18,191,'Mean ± SD; exact tests; Holm across 3 tests',fontsize=5.5)":"",
"heading('e','Neutrophil RNA response programs',82,122)":"heading('e','Neutrophil RNA response programs',66,89)",
"axbox(122,145,73,38)":"axbox(105,101,67,24)",
"axbox(198,146,2,36)":"axbox(175,101,1.5,24)",
"text(122,185,'Δ RNA program score; BH across 30 tests',fontsize=5.5);text(122,189,'Paired t-tests; n = 5–6 recorded donors',fontsize=5.5)":"",
"heading('g','Predicted ligand–target networks',111,199)":"heading('g','Predicted ligand–target networks',93,133)",
"text(116,205,'Top 3 prior-ranked targets per ligand',fontsize=6,color='#555')":"",
"T=pd.read_csv(OUT/'Figure3g_prior_targets.csv')":"T=pd.read_csv(ROOT/'outputs/15-Figure3-revised/Figure3g_prior_targets.csv')",
"axbox(113,215,92,53)":"axbox(96,143,84,22)",
"ax.text(.01,i+.77,'3f receptor: '+receptor,fontsize=4.6,color='#666',va='center')":"",
"text(116,270,'Dashed edges: ligand–target priors, not receptor-specific',fontsize=5.4)":"text(97,167,'Dashed edges: prior predictions; not measured activation',fontsize=5)",
"text(116,274,'signaling or measured treatment responses.',fontsize=5.4)":"",
}
for a,b in changes.items():
    assert a in old,a
    old=old.replace(a,b)
old=old.replace("ax.text(p,-1.17,t", "ax.text(p,-1.9,t")
old=old.replace("ax.tick_params(length=0,labelsize=5.4,pad=3)", "ax.tick_params(length=0,labelsize=5.4,pad=3);ax.set_ylabel('Neutrophil RNA modules',fontsize=5.5,labelpad=6)")
old=old.replace("y=i+.15+k*.29","y=i+.40")
old=old.replace("ax.annotate('',xy=(.67,y),xytext=(.30,i+.40),arrowprops=dict(arrowstyle='->',linestyle='--',color='#8899A8',lw=.6))\n        ax.text(.69,y,r.target,fontsize=5.8,va='center')", "ax.annotate('',xy=(.39+k*.205,y),xytext=(.25,i+.40),arrowprops=dict(arrowstyle='->',linestyle='--',color='#8899A8',lw=.5,connectionstyle='arc3,rad='+str((k-1)*.16)))\n        ax.text(.40+k*.205,y,r.target,fontsize=5,va='center',bbox=dict(fc='white',ec='none',pad=.2))")
start=old.index('    for k,(_,r) in enumerate(sub.iterrows()):')
end=old.index("text(97,167,",start)
old=old[:start]+"    ax.annotate('',xy=(.35,i+.40),xytext=(.21,i+.40),arrowprops=dict(arrowstyle='->',linestyle='--',color='#8899A8',lw=.6))\n    ax.text(.37,i+.40,'   '.join(sub.target),fontsize=5,va='center',bbox=dict(boxstyle='round,pad=.30',fc='#F7F9FB',ec='#CCD3DA',lw=.4))\n"+old[end:]
old=old.replace('i+.40','i+.50')
old=old[:old.index("heading('g',")]
ns={'__file__':__file__};exec(old,ns)
globals().update({k:ns[k] for k in ['fig','axbox','text','heading','plt','pd','np','OUT','cmap','TwoSlopeNorm']})
# Compact culture design, same experimental arms and independent baseline.
heading('a','Culture and blockade design',2,2)
ax=axbox(4,8,69,33);ax.set(xlim=(0,1),ylim=(1,0));ax.axis('off')
ax.text(.01,.05,'Peripheral-blood neutrophils',fontsize=6,fontweight='bold')
for i,(title,desc,col) in enumerate([('Alone','Untreated / NAMPTi / α5β1i','#F0F1F5'),('Ctrl FB','Untreated','#F0F5F7'),('UC FB','Untreated / NAMPTi / α5β1i / dual','#FAEEEE')]):
    y=.19+i*.20
    ax.add_patch(plt.Rectangle((0,y),1,.18,fc=col,ec='#B9C2CA',lw=.5))
    ax.text(.025,y+.06,title,fontsize=5.5,fontweight='bold',va='center');ax.text(.025,y+.135,desc,fontsize=5,va='center')
ax.text(.01,.87,'Independent preculture reference; RNA profiling',fontsize=5)
# Redraw f to retain its route order, labels and values at readable final type size.
heading('f','Predicted fibroblast–neutrophil routes',2,133)
S=Path('./output/Canonical_Nature_Revision_2026-09-03/source_data/figure3/Figure_3J_fibroblast_neutrophil_signaling_reversals.csv')
D=pd.read_csv(S);ligs=['COL1A2','FN1','IL1B','CEACAM1','CXCL12','C5'];route=['COL1A2 → CD44','FN1 → CD44','IL1B → IL1R2','CEACAM1 → CEACAM8','CXCL12 → CXCR4','C5 → C5AR1'];state=['ECM-FB → OSM','ECM-FB → CXCR4','CCN-FB → MX1/ISG','ECM-FB → PADI4','CCN-FB → CXCR4','ECM-FB → CXCR4']
ax=axbox(36,143,48,22);norm=TwoSlopeNorm(vmin=-.9,vcenter=0,vmax=.9)
for i,l in enumerate(ligs):
    for j,con in enumerate(['UF_vs_CF','UA5_vs_UF']):
        v=D[(D.ligand==l)&(D.contrast_id==con)].consensus_change.iloc[0]
        ax.add_patch(plt.Rectangle((j,i),1,1,fc=cmap(norm(v)),ec='white',lw=.4));ax.text(j+.5,i+.5,f'{v:+.2f}',fontsize=5.5,ha='center',va='center',color='white' if abs(v)>.5 else '#222')
    ax.text(-.05,i+.29,str(i+1)+'  '+route[i],ha='right',va='center',fontsize=5)
    ax.text(-.05,i+.75,state[i],ha='right',va='center',fontsize=5,color='#555')
ax.text(.5,-.17,'UC vs\nCtrl FB',ha='center',va='bottom',fontsize=5);ax.text(1.5,-.17,'α5β1i\nvs UC FB',ha='center',va='bottom',fontsize=5)
ax.set(xlim=(0,2),ylim=(6,0));ax.axis('off')
text(36,167,'LR consensus change (−0.9 to +0.9)',fontsize=5)
heading('g','Module support for ligands in f',93,133)
G=pd.read_csv(ROOT/'outputs/17-Figure3g-module-preview/Ligand_module_scores.csv')
modules=['OSM_inflammation','Type_I_interferon','Degranulation','Oxidative_burst','Chemotaxis','Hypoxia_glycolysis']
gl=['Inflam.\ncytokines','IFN\nresponse','Granule/\nprotease','NADPH\noxidase','Chemo-\ntaxis','Hypoxia/\nglycolysis']
ax=axbox(111,143,59,22);gcmap=plt.get_cmap('YlGnBu')
for i,l in enumerate(ligs):
    for j,m in enumerate(modules):
        r=G[(G.ligand==l)&(G.module==m)].iloc[0];v=(r.score+1)/2;n=int(r.top_decile_targets)
        ax.scatter(j,i,s=25+4*n,c=[gcmap(v)],ec='#49616D',lw=.3)
        ax.text(j,i,str(n),ha='center',va='center',fontsize=5,color='white' if v>.65 else '#222')
ax.set(xlim=(-.5,5.5),ylim=(5.5,-.5),xticks=range(6),yticks=range(6),yticklabels=[str(i+1)+'  '+l for i,l in enumerate(ligs)])
ax.set_xticklabels(gl,fontsize=5);ax.xaxis.tick_top();ax.tick_params(length=0,pad=3,labelsize=5)
ax.set_xticks(np.arange(-.5,6,1),minor=True);ax.set_yticks(np.arange(-.5,6,1),minor=True);ax.grid(which='minor',color='#E9EDF0',lw=.35);ax.tick_params(which='minor',length=0);ax.set_axisbelow(True)
for s in ax.spines.values():s.set_visible(False)
from matplotlib.colors import Normalize
cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=gcmap),cax=axbox(174,143,1.5,22));cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=5,length=1.5,pad=2)
from matplotlib.lines import Line2D
for i in range(6):
    y=143+(i+.5)*22/6
    fig.add_artist(Line2D([85/183,95/183],[1-y/170]*2,transform=fig.transFigure,color='#B4BDC5',lw=.45,ls=':'))
text(95,167,'Matched ligands; color: prior rank; dot size/count: top decile',fontsize=5)
fig.savefig(OUT/'Figure_3_Nature.pdf');fig.savefig(OUT/'Figure_3_Nature.svg');fig.savefig(OUT/'Figure_3_Nature.png',dpi=300)
plt.close(fig)
for p in (ROOT/'outputs/15-Figure3-revised').glob('*.csv'):shutil.copy2(p,OUT/p.name)
for p in (ROOT/'outputs/17-Figure3g-module-preview').glob('*.csv'):shutil.copy2(p,OUT/p.name)
legend=(ROOT/'outputs/15-Figure3-revised/Figure_3_legend.md').read_text(encoding='utf8')
legend=legend.replace('Original recorded culture design, preserved from the source Figure 3.','Recorded culture design, compactly redrawn from the source Figure 3.').replace('Original six selected predicted fibroblast-to-neutrophil routes, preserved without changing values, labels, colors or source artwork.','Six selected predicted fibroblast-to-neutrophil routes, redrawn at final size with the same route order, receiver states, numeric values and color limits.').replace('Receptor names are cross-references to panel f;','Receptor identities remain in panel f;')
legend+='\nLayout revision: 183 × 170 mm, following https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/. Detailed statistical and interpretive notes are in this legend. All plotted observations, heatmap values and significance annotations, and the 18 predicted targets are retained. Panel g groups targets horizontally to reduce height.\n'
start=legend.index('(g)');end=legend.index('For c and e:',start)
legend=legend[:start]+'''(g) Predicted ligand–module support. Matching row numbers 1–6 and dotted horizontal guides link each route in f to its ligand in g; guides indicate identity, not receptor-specific causality. Six ligands from panel f are compared with the exact six custom full-gene RNA modules in panel e using the local NicheNet ligand_target_matrix_nsga2r_final.rds prior. The background is the union of genes tested in Resolved-neutrophil pseudobulk analysis and present in the prior matrix. Each gene receives its within-ligand midrank percentile across that background, including zero weights and average ranks for ties. Color is the mean percentile across available module genes (0–1); it is not observed RNA expression or an activation direction. Dot area is 25 + 4 times the count of module genes at percentile >=0.90. Printed numbers give exact top-decile counts; a dot labeled 0 is not missing data. Counts depend on module size. Gene membership, coverage and per-gene weights are supplied. No significance or pathway-enrichment test was applied. Broad patterns can reflect shared prior structure and gene connectivity. This panel does not establish receptor-specific signaling, module activation, or inhibitor effects. It replaces the former display of 18 individual targets.\n\n'''+legend[end:]
legend=legend.replace('and the 18 predicted targets are retained. Panel g groups targets horizontally to reduce height.','in panels a–f are retained. Panel g now summarizes 36 ligand–module combinations as approved by the user.')
(OUT/'Figure_3_legend.md').write_text(legend,encoding='utf8')
from pypdf import PdfReader
page=PdfReader(OUT/'Figure_3_Nature.pdf').pages[0]
assert abs(float(page.mediabox.width)*25.4/72-183)<.01
assert abs(float(page.mediabox.height)*25.4/72-170)<.01
(OUT/'layout_verification.json').write_text(json.dumps({'width_mm':183,'height_mm':170,'facs_tests':45,'rna_tests':30,'assay_tests':3,'ligand_module_cells':36,'source_data_changed':False,'canonical_overwritten':False,'guidance':'https://research-figure-guide.nature.com/figures/building-and-exporting-figure-panels/'},indent=2))
print('Nature layout exported: 183 x 170 mm.')



