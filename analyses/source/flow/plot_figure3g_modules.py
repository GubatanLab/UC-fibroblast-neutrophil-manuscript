from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
OUT=ROOT/'outputs/17-Figure3g-module-preview'
D=pd.read_csv(OUT/'Ligand_module_scores.csv');genes=pd.read_csv(OUT/'Ligand_module_gene_support.csv')
mods=['OSM_inflammation','Type_I_interferon','Degranulation','Oxidative_burst','Chemotaxis','Hypoxia_glycolysis']
labels=['Inflammatory\ncytokines','Interferon\nresponse','Granule/\nprotease','NADPH\noxidase','Chemotaxis','Hypoxia/\nglycolysis'];ligs=['COL1A2','FN1','IL1B','CEACAM1','CXCL12','C5']
plt.rcParams.update({'font.family':'Arial','font.size':7,'svg.fonttype':'none'})
fig=plt.figure(figsize=(132/25.4,96/25.4));ax=fig.add_axes([.15,.31,.69,.47]);cmap=plt.get_cmap('YlGnBu');norm=Normalize(0,1)
for i,lig in enumerate(ligs):
    for j,m in enumerate(mods):
        r=D[(D.ligand==lig)&(D.module==m)].iloc[0];g=genes[(genes.ligand==lig)&(genes.module==m)]
        assert np.isclose(r.score,2*g.percentile.mean()-1) and r.top_decile_targets==int((g.percentile>=.9).sum())
        ax.scatter(j,i,s=18+11*r.top_decile_targets,c=[cmap((r.score+1)/2)],edgecolors='#49616D',linewidths=.35)
        ax.text(j,i,str(int(r.top_decile_targets)),ha='center',va='center',fontsize=5,color='white' if (r.score+1)/2>.65 else '#222')
ax.set(xticks=range(6),yticks=range(6),yticklabels=ligs,xlim=(-.5,5.5),ylim=(5.5,-.5))
ax.set_xticklabels(labels,fontsize=6);ax.xaxis.tick_top();ax.tick_params(length=0,pad=5);ax.set_ylabel('Ligands from panel 3f',fontsize=7,labelpad=8)
ax.set_xticks(np.arange(-.5,6,1),minor=True);ax.set_yticks(np.arange(-.5,6,1),minor=True);ax.grid(which='minor',color='#E9EDF0',lw=.5);ax.tick_params(which='minor',length=0);ax.set_axisbelow(True)
for s in ax.spines.values():s.set_visible(False)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.89,.33,.018,.42]));cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6,length=2);cb.set_label('Mean target-prior percentile',fontsize=6)
fig.text(.015,.966,'g',fontsize=10,fontweight='bold');fig.text(.065,.966,'Predicted ligand–module support',fontsize=10,fontweight='bold')
fig.text(.15,.88,'Neutrophil RNA modules from panel 3e',fontsize=7)
handles=[ax.scatter([],[],s=18+11*n,c='#8BAFC1',edgecolors='#49616D',lw=.35,label=str(n)) for n in [0,2,5,10]]
fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.16,.155),ncol=4,frameon=False,title='Dot size and number: targets in the top prior decile',title_fontsize=6,fontsize=6,handletextpad=.5,columnspacing=1.6)
fig.text(.05,.10,'Descriptive prior support; no activation direction or significance test.',fontsize=6.5)
fig.text(.05,.05,'Module size and shared network structure can influence this pattern.',fontsize=6.5)
for ext in ['png','svg','tiff']:
    fig.savefig(OUT/f'Figure_3g_module_preview.{ext}',**({'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {'dpi':300}))
plt.close(fig)
lines=['# Recommended Figure 3g: predicted ligand–module support','', 'Six ligands from panel 3f are compared against the exact six custom full-gene RNA modules in panel 3e. This is a descriptive prior-summary panel, not a pathway-enrichment significance analysis. It uses all ligand-target weights in the local NicheNet matrix, restricted to the union of genes tested in resolved-neutrophil pseudobulk analysis. Each gene receives its midrank percentile within this full common background for each ligand, including zero weights and average ranks for ties. Color shows the mean percentile across available module genes. The CSV score is 2 x that mean percentile minus 1; the plot displays the original 0–1 mean percentile.','', 'Dot area is 18 + 11 x the number of available module genes with percentile >=0.90. Printed numbers give exact counts; a small dot labeled 0 means no genes meet the top-decile threshold, not missing data. Module sizes differ; counts are not normalized and should be interpreted with membership coverage below. No genes were selected using observed treatment-effect direction or P values. No P values, BH correction, enrichment claims or activation signs are assigned. Predicted support is neither measured gene-module activity nor evidence that the indicated receptor in panel f mediates these targets. Broad shared patterns can reflect gene connectivity and prior architecture. This preview does not replace canonical Figure 3g.','', '## Available genes and highest-ranked examples','']
for m,l in zip(mods,labels):
    r=D[D.module==m].iloc[0];lines.append(f'- {l.replace(chr(10)," ")}: {int(r.n_module_genes)}/{int(r.n_requested)} genes available.')
lines+=['','Two highest-prior module genes per ligand/module are recorded in Ligand_module_scores.csv; these are descriptive examples, not independently validated downstream targets.']
(OUT/'Figure_3g_legend.md').write_text('\n'.join(lines),encoding='utf8')
print('36 ligand-module cells verified and exported.')
