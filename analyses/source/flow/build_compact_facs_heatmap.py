from pathlib import Path
import sys,re,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
OUT=ROOT/'outputs/08-compact-facs-panel';OUT.mkdir(exist_ok=True)
R=pd.read_csv(ROOT/'outputs/07-figure3-preview/FACS_sample_expression.csv')
markers=['OSM','CXCR4','MPO','MX1','PADI4'];groups=['N','N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK','N+FK','N+ATN']
def natural(s):return tuple(int(x) if x.isdigit() else x for x in re.split(r'(\d+)',s))
samples=[];blocks=[];positions=[];x=0
for g in groups:
    ss=sorted(R[R.group==g]['sample'].unique(),key=natural);start=x
    for s in ss:samples.append(s);positions.append(x);x+=1
    blocks.append((g,start,x,len(ss)));x+=.8
raw=R.pivot(index='marker',columns='sample',values='median_fi').loc[markers,samples]
assert raw.shape==(5,47) and raw.notna().all().all()
z=raw.sub(raw.mean(axis=1),axis=0).div(raw.std(axis=1,ddof=1),axis=0)
raw.to_csv(OUT/'Raw_sample_median_fluorescence.csv');z.to_csv(OUT/'Within_marker_z_scores.csv')
pd.DataFrame({'sample':samples,'group':[R[R['sample']==s].group.iloc[0] for s in samples],'position':range(1,48)}).to_csv(OUT/'Sample_column_key.csv',index=False)
limit=float(np.ceil(abs(z.to_numpy()).max()*2)/2)
plt.rcParams.update({'font.family':'Arial','font.size':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.5})
fig=plt.figure(figsize=(180/25.4,115/25.4));ax=fig.add_axes([.09,.63,.79,.225])
cmap=plt.get_cmap('RdBu_r');norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit)
for j,pos in enumerate(positions):
    for i,m in enumerate(markers):ax.add_patch(plt.Rectangle((pos,i),1,1,facecolor=cmap(norm(z.iloc[i,j])),edgecolor='white',linewidth=.25))
ax.set_xlim(0,x-.8);ax.set_ylim(5,0);ax.set_yticks(np.arange(5)+.5,markers);ax.set_xticks([]);ax.tick_params(axis='y',length=0,pad=5)
for sp in ax.spines.values():sp.set_visible(False)
labels=['N','N + NF','N + IF','N + IF\n+ FK','N + IF\n+ ATN','N + IF\n+ ATN + FK','N + FK','N + ATN']
for (g,start,end,n),label in zip(blocks,labels):
    center=(start+end)/2
    ax.text(center,-.32,label,ha='center',va='bottom',fontsize=6.5,linespacing=1.15,clip_on=False)
    ax.text(center,5.38,f'n = {n}',ha='center',va='top',fontsize=6.2,clip_on=False)
cax=fig.add_axes([.91,.645,.018,.195]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax);cb.set_label('Within-marker z score',fontsize=6.5,labelpad=4);cb.ax.tick_params(labelsize=6,length=2);cb.set_ticks([-limit,0,limit])
fig.text(.02,.97,'Neutrophil protein expression: samples and culture effects',fontsize=10,fontweight='bold',va='top')
fig.text(.02,.91,'a',fontsize=10,fontweight='bold')
contrasts=[('N+NF','N'),('N+IF','N'),('N+IF','N+NF'),('N+FK','N'),('N+ATN','N'),('N+IF+FK','N+IF'),('N+IF+ATN','N+IF'),('N+IF+ATN+FK','N+IF'),('N+IF+ATN+FK','N+IF+FK')]
effect_rows=[];effect=np.zeros((5,len(contrasts)))
for i,m in enumerate(markers):
    for j,(a,b) in enumerate(contrasts):
        va=R[(R.marker==m)&(R.group==a)].median_fi.to_numpy();vb=R[(R.marker==m)&(R.group==b)].median_fi.to_numpy()
        difference=float(va.mean()-vb.mean());scale=float(raw.loc[m].std(ddof=1));effect[i,j]=difference/scale
        za=z.loc[m,R[(R.marker==m)&(R.group==a)]['sample']].mean();zb=z.loc[m,R[(R.marker==m)&(R.group==b)]['sample']].mean();assert np.isclose(effect[i,j],za-zb)
        p=float(stats.permutation_test((va,vb),lambda a,b:np.mean(a)-np.mean(b),vectorized=False,n_resamples=np.inf,alternative='two-sided').pvalue)
        effect_rows.append(dict(marker=m,comparison=a,reference=b,n_comparison=len(va),n_reference=len(vb),mean_sample_median_comparison=va.mean(),mean_sample_median_reference=vb.mean(),difference_fluorescence=difference,all_sample_sd=scale,difference_mean_z=effect[i,j],p_exact=p))
effect_table=pd.DataFrame(effect_rows);ordered=effect_table.sort_values('p_exact',kind='stable').index;running=0
for rank,idx in enumerate(ordered):
    running=max(running,min(1,effect_table.at[idx,'p_exact']*(len(ordered)-rank)));effect_table.at[idx,'p_holm_45']=running
effect_table['q_bh_45']=stats.false_discovery_control(effect_table.p_exact.to_numpy(),method='bh')
effect_table['stars']=effect_table.q_bh_45.map(lambda p:'****' if p<.0001 else '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else '')
effect_table['displayed_significance']='BH FDR-adjusted exact permutation P across all 45 tests'
effect_table.to_csv(OUT/'Culture_and_blockade_effect_summary.csv',index=False)
ax2=fig.add_axes([.09,.205,.79,.245]);lim2=float(np.ceil(abs(effect).max()*2)/2);norm2=TwoSlopeNorm(vmin=-lim2,vcenter=0,vmax=lim2)
xx=[0,1,2,3.25,4.25,5.5,6.5,7.5,8.5]
for j,pos in enumerate(xx):
    for i,m in enumerate(markers):
        value=effect[i,j];stars=effect_table.iloc[i*len(contrasts)+j].stars;ax2.add_patch(plt.Rectangle((pos,i),1,1,facecolor=cmap(norm2(value)),edgecolor='white',lw=.5));ax2.text(pos+.5,i+.5,f'{value:+.2f}'+stars,ha='center',va='center',fontsize=6.3,color='white' if abs(value)>.58*lim2 else '#20252A')
ax2.set_xlim(0,9.5);ax2.set_ylim(5,0);ax2.set_yticks(np.arange(5)+.5,markers);ax2.tick_params(axis='y',length=0,pad=5);ax2.set_xticks([])
for sp in ax2.spines.values():sp.set_visible(False)
effect_labels=['N+NF\n− N','N+IF\n− N','N+IF\n− N+NF','N+FK\n− N','N+ATN\n− N','N+IF+FK\n− N+IF','N+IF+ATN\n− N+IF','N+IF+ATN+FK\n− N+IF','N+IF+ATN+FK\n− N+IF+FK']
for pos,label in zip(xx,effect_labels):ax2.text(pos+.5,-.25,label,ha='center',va='bottom',fontsize=5.7,clip_on=False)
for center,label in [(1.5,'Fibroblast exposure'),(4.25,'Blockade alone'),(7.5,'Blockade in inflamed co-culture')]:ax2.text(center,-1.1,label,ha='center',va='bottom',fontsize=6.2,fontweight='bold',clip_on=False)
fig.text(.02,.56,'b',fontsize=10,fontweight='bold');fig.text(.09,.56,'Effect summary: comparison minus reference',fontsize=8,fontweight='bold')
cax2=fig.add_axes([.91,.225,.018,.205]);cb2=fig.colorbar(plt.cm.ScalarMappable(norm=norm2,cmap=cmap),cax=cax2);cb2.set_label('Difference in mean z score',fontsize=6.5,labelpad=4);cb2.ax.tick_params(labelsize=6,length=2);cb2.set_ticks([-lim2,0,lim2])
any_significant=bool((effect_table.q_bh_45<.05).any())
fig.text(.09,.15,'BH FDR across 45 tests: * q < 0.05, ** q < 0.01, *** q < 0.001, **** q < 0.0001.',fontsize=6.2)
fig.text(.09,.105,'Top: 47 biological samples. Bottom: group mean differences, divided by each marker’s all-sample SD.',fontsize=6.2)
fig.text(.09,.06,'N, neutrophils; NF/IF, non-inflamed/inflamed fibroblasts; FK, FK-866; ATN, ATN-161.',fontsize=6.2)
for ext in ['pdf','svg','png','tiff']:
    kwargs={'dpi':300} if ext=='png' else {'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
    fig.savefig(OUT/f'Compact_FACS_expression_panel.{ext}',**kwargs)
plt.close(fig)
(OUT/'Panel_legend.md').write_text('''**Neutrophil protein expression across culture conditions.** Each column represents one biological sample and each row one marker. Tiles show the sample median compensated fluorescence across all gated neutrophils, standardized within each marker by subtracting the mean of all 47 sample medians and dividing by their sample standard deviation (ddof = 1). Samples receive equal weight; no values are clipped. The color scale is symmetric around zero and spans all observed z scores. Columns are ordered by recorded culture condition and natural sample ID, with the same order for every marker; no clustering or donor pairing is implied. n = 6 per condition except N+NF (n = 5; P5-5 absent). Raw fluorescence and a sample-column key are supplied separately. Colors compare samples within a marker and do not compare absolute protein abundance between markers.

N, neutrophils alone; NF/IF, non-inflamed/inflamed fibroblasts; FK, FK-866; ATN, ATN-161. Markers use the FCS assignments (MX1 B710, PADI4 B515); the protocol's reverse assignment remains documented in the analysis. Original parent gates and compensation are retained; 44/47 files lack viability measurements. This panel summarizes antibody fluorescence, not enzymatic activity, positive fractions or RNA-state abundance. No new statistical tests were performed.
''',encoding='utf-8')
assert np.allclose(z.mean(axis=1),0,atol=1e-12) and np.allclose(z.std(axis=1,ddof=1),1)
with (OUT/'Panel_legend.md').open('a',encoding='utf-8') as f:f.write('\n**Lower panel: fibroblast-culture and blockade effects.** Each tile is the group mean sample-median fluorescence in the first listed condition minus that in the second, divided by the marker-specific sample SD across all 47 biological samples. Equivalently, it is the difference between group mean z scores from the upper panel. Values are annotated to two decimals and use a separate symmetric, unclipped color scale. Red indicates an increase and blue a decrease. The nine explicit contrasts cover three fibroblast-exposure comparisons, two blockade-alone comparisons and four blockade comparisons in inflamed co-culture. The final contrast isolates addition of ATN-161 in the FK-866 background. No non-inflamed-fibroblast blockade arm is available. Effects are descriptive; they are not fold changes, interaction tests, calibrated abundance differences or measures of statistical significance. No tests were added. All significance symbols have been removed.\n')
legend_path=OUT/'Panel_legend.md';legend=legend_path.read_text(encoding='utf-8').replace('No new statistical tests were performed.','').replace('No tests were added. All significance symbols have been removed.','Significance was tested on unstandardized biological-sample medians using exact two-sided unpaired permutation tests of group mean differences. All allocations were enumerated (924 for 6 versus 6; 462 for 5 versus 6). Holm correction covers all 45 tests displayed in the lower panel. Stars denote adjusted P < 0.05 (*), < 0.01 (**), < 0.001 (***) or < 0.0001 (****). Unmarked tiles do not meet adjusted P < 0.05, which does not establish equivalence. This family differs from the narrower five-marker family used previously. Pairing is unconfirmed; the tests assume exchangeability, and plate confounding remains unresolved.');legend_path.write_text(legend,encoding='utf-8')
legend=legend_path.read_text(encoding='utf-8').replace('Holm correction covers all 45 tests displayed in the lower panel. Stars denote adjusted P < 0.05 (*), < 0.01 (**), < 0.001 (***) or < 0.0001 (****). Unmarked tiles do not meet adjusted P < 0.05, which does not establish equivalence. This family differs from the narrower five-marker family used previously.','At the user\'s request, displayed stars use unadjusted exact P values, with no multiple-testing adjustment: P < 0.05 (*), < 0.01 (**), < 0.001 (***) or < 0.0001 (****). Unmarked tiles do not meet unadjusted P < 0.05, which does not establish equivalence. These are exploratory nominal significance annotations across 45 comparisons. The previously calculated Holm values remain in Source Data for traceability but are not used for the stars.');legend_path.write_text(legend,encoding='utf-8')
legend=legend_path.read_text(encoding='utf-8');start=legend.index("At the user's request, displayed stars")
end=legend.index(' Pairing is unconfirmed;',start)
legend=legend[:start]+'Displayed stars use Benjamini-Hochberg FDR-adjusted exact P values across all 45 tests: q < 0.05 (*), < 0.01 (**), < 0.001 (***) or < 0.0001 (****). Unmarked tiles do not meet q < 0.05, which does not establish equivalence. BH was selected by the user after comparison of correction methods; this is an exploratory post hoc analysis, not a prospectively specified testing plan. BH requires independence or appropriate positive dependence among valid P values; arbitrary-dependence BY sensitivity results are retained in the separate FDR comparison report. The prior Holm results remain in Source Data but are not used for annotations.'+legend[end:]
legend_path.write_text(legend,encoding='utf-8')
effect_table[['marker','comparison','reference','p_exact','q_bh_45','stars','displayed_significance']].to_csv(OUT/'Significance_annotations.csv',index=False)
print('BH significant tiles:',int((effect_table.q_bh_45<.05).sum()),'minimum BH q:',effect_table.q_bh_45.min())
print(json.dumps({'samples':len(samples),'markers':len(markers),'sample_tiles':raw.size,'effect_tiles':effect.size,'color_limits':[-limit,limit],'effect_color_limits':[-lim2,lim2],'size_mm':[180,115],'checks_passed':True}))
