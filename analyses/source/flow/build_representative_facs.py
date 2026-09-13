from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT=ROOT/'outputs/13-representative-FACS';OUT.mkdir(exist_ok=True)
CACHE=ROOT/'analysis/regating/events'
channels={'OSM':'R670-A','CXCR4':'V610-A','MX1':'B710-A','PADI4':'B515-A','MPO':'YG610-A'}
M=pd.read_csv(ROOT/'outputs/08-compact-facs-panel/Raw_sample_median_fluorescence.csv',index_col=0).loc[list(channels)]
groups={'Neutrophils alone':[f'P1-{i}' for i in range(1,7)],'Ctrl FB':[f'P5-{i}' for i in [1,2,3,4,6]],'UC FB':[f'P3-{i}' for i in range(1,7)],'UC FB + α5β1i':[f'P4-{i}' for i in range(1,7)]}
selection=[];selected={}
# One sample per arm selected jointly across all five markers, not separately per marker.
for g,samples in groups.items():
    x=M[samples].T;center=x.median();scale=x.std(ddof=1)
    distances=(((x-center)/scale)**2).mean(axis=1)
    chosen=distances.idxmin();selected[g]=chosen
    for s in samples:selection.append({'group':g,'sample':s,'standardized_distance_squared':distances[s],'selected':s==chosen})
pd.DataFrame(selection).to_csv(OUT/'Representative_sample_selection.csv',index=False)
arrays={};audit=[]
for g,s in selected.items():
    meta=json.loads((CACHE/f'co_{s}.json').read_text());z=np.load(CACHE/f'co_{s}.npz')
    raw=ROOT/meta['file'];digest=hashlib.sha256(raw.read_bytes()).hexdigest()
    inventory=pd.read_csv(ROOT/'analysis/results/raw_fcs_inventory.csv').set_index('file')
    assert digest==inventory.loc[meta['file'],'sha256']
    for marker,ch in channels.items():
        v=z['events'][z['old_3'],meta['channels'].index(ch)]
        assert np.isfinite(v).all() and np.isclose(np.median(v),M.loc[marker,s],atol=1e-8)
        arrays[g,marker]=v
        audit.append({'group':g,'sample':s,'marker':marker,'channel':ch,'events':len(v),'median_fluorescence':float(np.median(v)),'raw_sha256':digest,'source_file':str(raw)})
pd.DataFrame(audit).to_csv(OUT/'Representative_event_audit.csv',index=False)
# Retain the original complete nine-contrast BH heatmap without re-estimating its statistics.
src=(ROOT/'analysis/build_single_bh_heatmap.py').read_text(encoding='utf-8').split("for ext in ['pdf'")[0]
src=src.replace("markers=['OSM','CXCR4','MPO','MX1','PADI4']","markers=['OSM','CXCR4','MX1','PADI4','MPO']")
ns={'__file__':str(ROOT/'analysis/build_single_bh_heatmap.py')};exec(src,ns)
fig=ns['fig'];fig.set_size_inches(183/25.4,148/25.4)
factor=76/148
for ax in fig.axes:
    p=ax.get_position();ax.set_position([p.x0,p.y0*factor,p.width,p.height*factor])
for t in fig.texts:
    x,y=t.get_position();t.set_position((x,y*factor))
colors=['#555555','#B58B24','#B44B45','#247BA0'];histrows=[]
for i,marker in enumerate(channels):
    ax=fig.add_axes([.07+i*.183,.635,.151,.205])
    transformed=[np.arcsinh(arrays[g,marker]/150) for g in groups]
    allv=np.concatenate(transformed);lo,hi=np.quantile(allv,[.001,.999]);margin=(hi-lo)*.025;edges=np.linspace(lo-margin,hi+margin,65)
    for (g,col),v in zip(zip(groups,colors),transformed):
        counts,_=np.histogram(v,bins=edges);percent=100*counts/len(v)
        ax.stairs(percent,edges,color=col,lw=1.1,label=g)
        for k,count in enumerate(counts):histrows.append({'marker':marker,'group':g,'sample':selected[g],'bin_left_asinh':edges[k],'bin_right_asinh':edges[k+1],'count':int(count),'percent_all_parent_events':percent[k]})
        ax.axvline(np.arcsinh(np.median(arrays[g,marker])/150),color=col,ls=':',lw=.65,alpha=.7)
    ticks=[-1000,0,1000,10000,100000];ticks=[t for t in ticks if edges[0]<=np.arcsinh(t/150)<=edges[-1]]
    ax.set_xticks(np.arcsinh(np.array(ticks)/150),[str(t) if abs(t)<10000 else f'{t//1000}k' for t in ticks])
    ax.tick_params(labelsize=5.5,length=2,pad=2);ax.set_xlim(edges[0],edges[-1]);ax.set_ylim(bottom=0)
    ax.set_title(marker,fontsize=8,fontweight='bold');ax.set_xlabel(channels[marker],fontsize=6)
    if i==0:ax.set_ylabel('Neutrophils / bin (%)',fontsize=6)
    ax.spines[['top','right']].set_visible(False)
fig.text(.025,.963,'Representative FACS: fibroblast exposure and α5β1 blockade',fontsize=10,fontweight='bold')
fig.legend(handles=[Line2D([0],[0],color=c,lw=1.5,label=f'{g} ({selected[g]})') for g,c in zip(groups,colors)],loc='upper center',bbox_to_anchor=(.5,.936),ncol=2,frameon=False,fontsize=7)
fig.text(.07,.565,'All gated neutrophils. Dotted lines: sample medians. Compensated fluorescence; arcsinh scale (cofactor 150).',fontsize=6.3)
fig.text(.07,.538,'One representative sample per arm, selected jointly across five markers; cohort inference below uses all samples.',fontsize=6.3)
for ext in ['png','svg','tiff']:
    fig.savefig(OUT/f'Representative_FACS_and_BH_heatmap.{ext}',**({'dpi':600,'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {'dpi':300}))
plt.close(fig)
pd.DataFrame(histrows).to_csv(OUT/'Histogram_source_data.csv',index=False)
legend=f'''Representative FACS distributions with the complete BH heatmap. Neutrophils alone ({selected['Neutrophils alone']}), control fibroblast coculture ({selected['Ctrl FB']}), UC fibroblast coculture untreated ({selected['UC FB']}) and UC fibroblast coculture with ATN-161 ({selected['UC FB + α5β1i']}). Each arm's representative sample minimizes the mean squared standardized distance of its five sample medians from the within-arm componentwise median; standardization uses within-arm sample SD. Selection is joint across OSM, CXCR4, MPO, MX1 and PADI4. No matched-donor pairing is assumed. All candidate samples and selection distances are supplied (n=6 per arm except Ctrl FB n=5).\n\nHistograms include all events in the retained scatter/singlet/CD16/CD11b neutrophil gate (old_3), without a CXCR4 or marker-positive restriction. Cached compensated event medians were verified against the source table, and original FCS SHA256 values were verified. Counts are normalized to percent of all parent events, using 64 identical bins per marker across the four displayed samples in arcsinh(fluorescence/150) space. X-axis labels are untransformed compensated fluorescence. The displayed range uses the pooled 0.1–99.9th transformed percentiles with 2.5% range padding; extreme tails can fall outside the display. No smoothing, positive thresholds or event subsampling were used. Dotted lines mark sample median fluorescence. Different marker panels have their own fluorescence and frequency ranges; the four conditions share both axes within each marker.\n\nThe full nine-contrast heatmap and its 45-test BH correction are unchanged, with 31 significant comparisons. Histograms are illustrative, not biological replication; inferential statistics use biological-sample medians. MPO antibody intensity is not enzyme activity. FCS MX1/PADI4 labels conflict with protocol assignments; this remains unresolved. Parent gates/compensation were retained and most co-culture files lack viability measurements. Treatment/plate confounding remains. No new claim of direct receptor specificity or NET formation is implied.\n'''
(OUT/'Figure_legend.md').write_text(legend,encoding='utf8')
print(json.dumps({'selected':selected,'events_verified':True,'facs_bh_tests':len(ns['D'])},indent=2))


