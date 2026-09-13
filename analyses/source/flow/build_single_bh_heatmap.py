from pathlib import Path
import sys,json,shutil
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
OUT=ROOT/'outputs/10-single-BH-heatmap';OUT.mkdir(exist_ok=True)
D=pd.read_csv(ROOT/'outputs/08-compact-facs-panel/Culture_and_blockade_effect_summary.csv')
markers=['OSM','CXCR4','MPO','MX1','PADI4'];pairs=list(dict.fromkeys(zip(D.comparison,D.reference)))
labels=['Ctrl FB\nvs alone','UC FB\nvs alone','UC FB\nvs Ctrl FB','NAMPTi\nvs untreated','α5β1i\nvs untreated','NAMPTi\nvs untreated','α5β1i\nvs untreated','Dual blockade\nvs untreated','Dual blockade\nvs NAMPTi']
positions=[0,1,2,3.25,4.25,5.5,6.5,7.5,8.5]
plt.rcParams.update({'font.family':'Arial','font.size':7,'svg.fonttype':'none','pdf.fonttype':42})
fig=plt.figure(figsize=(183/25.4,76/25.4));ax=fig.add_axes([.085,.25,.80,.43]);cmap=plt.get_cmap('RdBu_r');norm=TwoSlopeNorm(vmin=-3,vcenter=0,vmax=3)
for j,((a,b),pos) in enumerate(zip(pairs,positions)):
    for i,m in enumerate(markers):
        r=D[(D.marker==m)&(D.comparison==a)&(D.reference==b)].iloc[0];q=r.q_bh_45;stars='****' if q<.0001 else '***' if q<.001 else '**' if q<.01 else '*' if q<.05 else '';value=r.difference_mean_z
        ax.add_patch(plt.Rectangle((pos,i),1,1,facecolor=cmap(norm(value)),edgecolor='white',lw=.5));ax.text(pos+.5,i+.5,f'{value:+.2f}{stars}',ha='center',va='center',fontsize=6.5,color='white' if abs(value)>1.74 else '#222222')
    ax.text(pos+.5,-.22,labels[j],ha='center',va='bottom',fontsize=6.2,clip_on=False)
for center,title in [(1.5,'Fibroblast exposure'),(4.25,'Neutrophils alone'),(7.5,'Blockade in UC FB coculture')]:ax.text(center,-1.1,title,ha='center',va='bottom',fontsize=7,fontweight='bold',clip_on=False)
ax.set_xlim(0,9.5);ax.set_ylim(5,0);ax.set_yticks(np.arange(5)+.5,markers);ax.set_xticks([]);ax.tick_params(length=0,pad=5)
for s in ax.spines.values():s.set_visible(False)
cax=fig.add_axes([.91,.275,.018,.375]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax);cb.set_ticks([-3,0,3]);cb.set_label('Difference in mean z score',fontsize=6.5);cb.ax.tick_params(labelsize=6,length=2)
fig.text(.02,.955,'Fibroblast exposure and blockade effects on neutrophil proteins',fontsize=10,fontweight='bold',va='top')
fig.text(.085,.17,'Values: comparison minus reference. Stars: BH FDR across 45 tests; * q < 0.05, ** q < 0.01.',fontsize=6.4)
fig.text(.085,.10,'Ctrl FB, control fibroblasts; UC FB, UC fibroblasts; NAMPTi, FK-866; α5β1i, ATN-161.',fontsize=6.4)
fig.text(.085,.04,'Dual blockade, NAMPTi + α5β1i. All-neutrophil sample medians; n = 6/group except Ctrl FB (n = 5).',fontsize=6.4)
for ext in ['pdf','svg','png','tiff']:
    kwargs={'dpi':300} if ext=='png' else {'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(OUT/f'BH_FACS_effect_heatmap.{ext}',**kwargs)
plt.close(fig)
D.to_csv(OUT/'BH_effect_source_data.csv',index=False)
pd.DataFrame([{'column':i+1,'heading':labels[i].replace('\n',' '),'comparison_source':a,'reference_source':b,'context':'Fibroblast exposure' if i<3 else 'Neutrophils alone' if i<5 else 'UC fibroblast coculture'} for i,(a,b) in enumerate(pairs)]).to_csv(OUT/'Column_label_key.csv',index=False)
(OUT/'Heatmap_legend.md').write_text('''**Fibroblast exposure and blockade effects on neutrophil proteins.** Tiles show the mean of biological-sample median compensated fluorescence in the comparison condition minus that in the reference condition, divided by each marker's sample SD across all 47 biological samples. This equals a difference in group mean within-marker z scores. Red denotes increased and blue decreased fluorescence. Numbers show standardized effects; stars indicate BH FDR-adjusted exact two-sided unpaired permutation P values across all 45 tests (* q < 0.05; ** q < 0.01). Tests used unstandardized biological-sample medians. All neutrophils in the retained parent gate were included, without CXCR4 restriction. n = 6 per condition except control FB coculture (n = 5). All 31 BH-significant comparisons are retained.

Terminology follows the UC neutrophil-fibroblast Figure 3: Ctrl FB, control fibroblasts (source non-inflamed fibroblasts); UC FB, UC fibroblasts (source inflamed fibroblasts); NAMPTi, FK-866; α5β1i, ATN-161; dual blockade, FK-866 plus ATN-161. The block headings define context. In the first block, 'alone' means untreated neutrophils alone. In the middle block, 'untreated' means untreated neutrophils alone. In the right block, 'untreated' means untreated UC fibroblast coculture. The last column compares dual blockade with NAMPTi in UC fibroblast coculture. Column_label_key.csv preserves the exact source conditions for every column.

This is a display revision of the canonical BH results, with no change to measurements, contrasts, adjustment family or P values. It does not establish receptor-specific causality, enzymatic activity or NET formation. MX1/PADI4 use the FCS assignments; the conflicting protocol annotation remains unresolved. Original parent gates and compensation were retained; 44/47 files lack viability measurements. Matching across conditions is unconfirmed and treatment/plate confounding remains. BH assumptions and BY sensitivity results are documented in the FDR comparison report. No canonical assembled Figure 3 was overwritten.
''',encoding='utf-8')
assert len(D)==45 and (D.q_bh_45<.05).sum()==31
print('Single heatmap exported; 45 effects and 31 BH discoveries preserved.')
