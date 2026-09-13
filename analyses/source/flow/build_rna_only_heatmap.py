from pathlib import Path
import sys,shutil
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
SRC=ROOT/'outputs/11-FACS-RNA-module-heatmap'
OUT=ROOT/'outputs/12-RNA-module-heatmap';OUT.mkdir(exist_ok=True)
D=pd.read_csv(SRC/'RNA_module_effects_BH30.csv')
mods=['OSM_inflammation','Type_I_interferon','Degranulation','Oxidative_burst','Chemotaxis','Hypoxia_glycolysis']
names=['Inflammatory cytokines','Interferon response','Granule/protease program','NADPH oxidase machinery','Chemotaxis','Hypoxia/glycolysis']
cons=['control_fibroblasts_vs_alone','UC_vs_alone','UC_vs_control_fibroblasts','alpha_alone','alpha_in_UC']
labels=['Ctrl FB\nvs alone','UC FB\nvs alone','UC FB\nvs Ctrl FB','α5β1i\nvs untreated','α5β1i\nvs untreated']
positions=[0,1,2,3.25,4.5]
plt.rcParams.update({'font.family':'Arial','font.size':7,'svg.fonttype':'none','pdf.fonttype':42})
fig=plt.figure(figsize=(183/25.4,100/25.4));ax=fig.add_axes([.30,.235,.57,.49])
cmap=plt.get_cmap('RdBu_r');norm=TwoSlopeNorm(vmin=-.6,vcenter=0,vmax=.6)
for j,(con,pos) in enumerate(zip(cons,positions)):
    for i,m in enumerate(mods):
        r=D[(D.program==m)&(D.contrast==con)].iloc[0]
        v=r.effect;q=r.q_bh_30;s='***' if q<.001 else '**' if q<.01 else '*' if q<.05 else ''
        ax.add_patch(plt.Rectangle((pos,i),1,1,fc=cmap(norm(v)),ec='white',lw=.7))
        ax.text(pos+.5,i+.5,f'{v:+.2f}{s}',ha='center',va='center',fontsize=7,color='white' if abs(v)>.348 else '#222')
    ax.text(pos+.5,-.18,labels[j],ha='center',va='bottom',fontsize=6.7,clip_on=False)
for center,title in [(1.5,'Fibroblast exposure'),(3.75,'Neutrophils\nalone'),(5,'UC FB\ncoculture')]:
    ax.text(center,-.95,title,ha='center',va='bottom',fontsize=7,fontweight='bold',clip_on=False)
ax.set(xlim=(0,5.5),ylim=(6,0),xticks=[],yticks=np.arange(6)+.5,yticklabels=names)
ax.tick_params(length=0,pad=6)
for s in ax.spines.values():s.set_visible(False)
cax=fig.add_axes([.902,.25,.015,.46]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_ticks([-.6,0,.6]);cb.ax.tick_params(labelsize=6,length=2);cb.set_label('Δ RNA program score',fontsize=7)
fig.text(.025,.955,'Neutrophil RNA response programs',fontsize=11,fontweight='bold')
fig.text(.025,.904,'Fibroblast exposure and α5β1 blockade',fontsize=8,color='#444')
fig.text(.025,.166,'Mean paired donor differences (n = 5–6). BH across 30 tests: * q < 0.05; ** q < 0.01; *** q < 0.001.',fontsize=6.6)
fig.text(.025,.117,'Stars use paired t-tests. Exact paired sign-flip sensitivity: no comparisons pass BH q < 0.05.',fontsize=6.6)
fig.text(.025,.068,'Ctrl FB, control fibroblasts; UC FB, UC fibroblasts; α5β1i, ATN-161. Scores describe RNA programs, not measured function.',fontsize=6.4)
for ext in ['png','svg','tiff']:
    fig.savefig(OUT/f'RNA_modules_BH_heatmap.{ext}',**({'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {'dpi':300}))
plt.close(fig)
for name in ['RNA_module_effects_BH30.csv','RNA_module_gene_membership.csv','RNA_donor_differences.csv']:
    shutil.copy2(SRC/name,OUT/name)
(OUT/'Figure_legend.md').write_text('''Neutrophil RNA response programs. Six existing Figure 3 full-gene programs are shown for the five available contrasts. Values are unchanged mean paired donor differences in source RNA program scores. Resolved neutrophils with at least 20 cells per donor-condition arm were included; recorded donor n is 5, 6, 5, 5 and 5 across columns. Scores use the existing log-normalization and equally weighted untreated donor-condition reference standardization. Gene membership is supplied; these are custom source programs. Oxidative-burst machinery includes MPO alongside NADPH oxidase components.\n\nStars retain BH-adjusted two-sided paired t-tests across all 30 displayed tests (six modules by five contrasts): * q<0.05, ** q<0.01, *** q<0.001. Sixteen comparisons pass this threshold. Exact paired sign-flip tests with BH across 30 tests yield no discoveries. Recorded donor matching remains unverified. No values, tests, correction family or gene definitions were changed for this standalone display; the four unavailable RNA columns were removed.\n\nCtrl FB, control fibroblasts; UC FB, UC fibroblasts; alpha5beta1i, ATN-161. Alone denotes untreated neutrophils alone. Untreated in the fourth column denotes untreated neutrophils alone; in the fifth it denotes untreated UC fibroblast coculture. Red indicates increased and blue decreased RNA program scores. Scores describe transcriptional programs and do not establish functional capacity. The assembled manuscript Figure 3 and prior FACS figures were not overwritten.\n''',encoding='utf8')
assert len(D)==30 and (D.q_bh_30<.05).sum()==16
print('Standalone RNA heatmap exported: six modules, five complete columns; 30 effects unchanged.')
