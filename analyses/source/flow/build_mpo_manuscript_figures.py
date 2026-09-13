from pathlib import Path
import sys,json,shutil,zipfile,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
SRC=ROOT/'outputs/03-mpo-blockade';OUT=ROOT/'outputs/04-mpo-manuscript';OUT.mkdir(exist_ok=True)
R=pd.read_csv(SRC/'MPO_sample_measurements.csv');T=pd.read_csv(SRC/'MPO_blockade_comparisons.csv');T=T[T.metric=='median_mpo']
CACHE=ROOT/'analysis/regating/events';ORDER=['N','N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK','N+FK','N+ATN']
PAIRS=[('N+IF','N+IF+ATN'),('N+IF+FK','N+IF+ATN+FK'),('N','N+ATN')]
TITLES=['Inflamed co-culture','Inflamed co-culture + FK-866','Neutrophils alone']
COLORS=['#0072B2','#D55E00'];EVENTS={};REP={};selection=[]
for row in R.itertuples():
    m=json.loads((CACHE/f'co_{row.sample}.json').read_text());z=np.load(CACHE/f'co_{row.sample}.npz')
    v=z['events'][z['old_3'],m['channels'].index('YG610-A')]
    assert len(v)==row.neutrophil_events and np.isclose(np.median(v),row.median_mpo)
    EVENTS[row.sample]=v
for g in ORDER:
    x=R[R.group==g].copy();target=x.median_mpo.median();x['distance']=abs(x.median_mpo-target)
    row=x.sort_values(['distance','sample']).iloc[0];REP[g]=row['sample']
    selection.append(dict(group=g,sample=row['sample'],sample_median=row.median_mpo,group_median_of_sample_medians=target,neutrophil_events=int(row.neutrophil_events),rule='Minimum absolute distance to group median of sample medians; sample ID breaks ties'))
pd.DataFrame(selection).to_csv(OUT/'Representative_samples.csv',index=False)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.7,'lines.linewidth':1.2,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
# Identical bins and display limits across every fluorescence panel, containing all events.
allv=np.concatenate([np.arcsinh(v/150) for v in EVENTS.values()]);edges=np.linspace(allv.min()-.02,allv.max()+.02,150)
centers=(edges[1:]+edges[:-1])/2;ticks=np.array([-1000,0,1000,10000,100000]);tickpos=np.arcsinh(ticks/150)
def histogram(ax,v,color,label=None,alpha=1,lw=1.2):
    h=np.histogram(np.arcsinh(v/150),bins=edges)[0].astype(float);h=100*h/h.max()
    ax.step(centers,h,where='mid',color=color,label=label,alpha=alpha,lw=lw)
def histaxes(ax):
    ax.set_xlim(edges[0],edges[-1]);ax.set_ylim(0,108);ax.set_yticks([0,50,100]);ax.set_xticks(tickpos,['−10³','0','10³','10⁴','10⁵']);ax.set_xlabel('MPO fluorescence (a.u.)');ax.set_ylabel('Normalized count (% of peak)')
def letter(ax,s):ax.text(-.22,1.12,s,transform=ax.transAxes,fontsize=12,fontweight='bold')
def export(fig,name):
    fig.savefig(OUT/f'{name}.svg')
    fig.savefig(OUT/f'{name}.png',dpi=300)
    fig.savefig(OUT/f'{name}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
with PdfPages(OUT/'MPO_manuscript_figures.pdf') as pdf:
    fig,axes=plt.subplots(2,3,figsize=(7.2,5.8))
    for i,((ref,trt),title) in enumerate(zip(PAIRS,TITLES)):
        ax=axes[0,i]
        for g,c,label in [(ref,COLORS[0],'No ATN-161'),(trt,COLORS[1],'+ ATN-161')]:histogram(ax,EVENTS[REP[g]],c,label)
        histaxes(ax);ax.set_title(title,pad=27);ax.legend(loc='upper center',bbox_to_anchor=(.5,1.14),frameon=False,ncol=1,handlelength=1.6,labelspacing=.25)
        ax.text(.5,-.33,f'{REP[ref]} / {REP[trt]}',transform=ax.transAxes,ha='center',fontsize=7,color='#555555');letter(ax,chr(65+i))
        ax=axes[1,i]
        for j,g in enumerate([ref,trt]):
            values=R[R.group==g].sort_values('sample').median_mpo.to_numpy()
            ax.scatter(j+np.linspace(-.09,.09,len(values)),values,s=20,color=COLORS[j],edgecolors='white',linewidths=.4,zorder=3)
            ax.errorbar(j,values.mean(),yerr=values.std(ddof=1),fmt='_',color='black',markersize=12,capsize=3,lw=.8)
        row=T[(T.reference==ref)&(T.treated==trt)].iloc[0]
        ax.plot([0,0,1,1],[5000,5150,5150,5000],color='black',lw=.7);ax.text(.5,5250,f'Adjusted P = {row.p_holm:.4f}',ha='center',fontsize=8)
        ax.set_ylim(-650,5700);ax.set_yticks([0,2000,4000]);ax.set_xlim(-.45,1.45);ax.set_xticks([0,1],['No ATN-161','+ ATN-161']);ax.set_ylabel('Median MPO fluorescence (a.u.)');letter(ax,chr(68+i))
    fig.subplots_adjust(left=.10,right=.98,top=.86,bottom=.10,wspace=.65,hspace=.70)
    pdf.savefig(fig);export(fig,'Figure_MPO_blockade');plt.close(fig)
    fig=plt.figure(figsize=(7.2,7.0));gs=fig.add_gridspec(3,3,height_ratios=[1.5,1,1],left=.11,right=.98,top=.95,bottom=.08,hspace=.65,wspace=.65)
    ax=fig.add_subplot(gs[0,:]);labels=['N','N + NF','N + IF','N + IF\n+ FK','N + IF\n+ ATN','N + IF\n+ ATN + FK','N + FK','N + ATN']
    for i,g in enumerate(ORDER):
        x=R[R.group==g].sort_values('sample').median_mpo.to_numpy();color=COLORS[1] if 'ATN' in g else COLORS[0]
        ax.scatter(i+np.linspace(-.09,.09,len(x)),x,s=19,color=color,edgecolors='white',linewidths=.4,zorder=3)
        ax.errorbar(i,x.mean(),yerr=x.std(ddof=1),fmt='_',color='black',capsize=3,markersize=10,lw=.8)
        ax.text(i,5100,f'n = {len(x)}',ha='center',fontsize=7)
    ax.set_ylim(-650,5500);ax.set_xticks(range(8),labels);ax.set_ylabel('Median MPO fluorescence (a.u.)');ax.text(-.065,1.025,'A',transform=ax.transAxes,fontsize=12,fontweight='bold')
    for i,(ref,trt) in enumerate(PAIRS):
        for j,g in enumerate([ref,trt]):
            ax=fig.add_subplot(gs[1+j,i]);values=R[R.group==g]
            for sample in values['sample']:histogram(ax,EVENTS[sample],COLORS[j],alpha=.3,lw=.7)
            histogram(ax,EVENTS[REP[g]],COLORS[j],lw=1.4)
            histaxes(ax);ax.set_title(g,pad=5);letter(ax,chr(66+j*3+i))
    pdf.savefig(fig);export(fig,'Supplement_MPO_all_samples');plt.close(fig)
legend='''# Figure legend: ATN-161 and MPO staining in neutrophils

**Figure. ATN-161 treatment is associated with lower MPO antibody fluorescence in co-cultured neutrophils.** (A-C) Representative flow-cytometry histograms of MPO antibody staining in all gated neutrophils from inflamed fibroblast co-culture (A), inflamed co-culture with FK-866 (B), and neutrophils cultured alone (C), without (blue) or with (orange) ATN-161. Sample IDs are shown beneath each panel. Histograms display compensated YG610-A fluorescence using an inverse-hyperbolic-sine transform with cofactor 150; ticks indicate fluorescence in arbitrary units. All panels use identical bin edges and fluorescence limits. Each histogram is normalized independently to its maximum bin count (100%); histogram height does not represent cell abundance. Representative samples were selected by the smallest absolute distance between their sample median MPO fluorescence and the group median of sample medians; sample ID breaks ties. Representative samples across conditions are not asserted to be paired. (D-F) Median MPO fluorescence for all biological samples in the corresponding conditions. Each point represents one independent biological sample; bars show mean ± SD, n = 6 per condition. All neutrophils in the retained scatter/singlet/CD16/CD11b parent gate were included, without restriction to CXCR4-positive or MPO-positive cells. Exact two-sided unpaired permutation tests compared group means of sample medians (924 allocations per comparison); P values were adjusted by Holm's method across the three ATN-161 contrasts. Negative compensated fluorescence values were retained. MPO was measured by antibody staining and is an abundance-related readout, not a measurement of MPO catalytic activity.

**Supplementary figure. MPO staining across all conditions and individual biological samples.** (A) Median MPO fluorescence across all eight available conditions, with each point representing one sample and bars showing mean ± SD. n = 6 per condition except N + NF (n = 5; P5-5 absent). N, neutrophils alone; NF, non-inflamed fibroblasts; IF, inflamed fibroblasts; FK, FK-866; ATN, ATN-161. (B-G) All individual sample histograms for the six groups in the primary blockade comparisons. Thin transparent lines show every sample; the thick line identifies the representative sample used in the main figure. Histogram transforms, bins and normalization match the main figure. No hypothesis tests are added for the remaining descriptive conditions.

## Internal submission notes

These are manuscript-layout figures derived from the reviewed analysis, not a certification that staining or parent gates have been validated. A viability channel is absent in 44 of 47 FCS files, so the population cannot consistently be described as live neutrophils. The original parent gates and compensation matrix were retained. MPO clone/conjugate, staining conditions, inhibitor doses, culture duration and matching/experimental blocks should be completed in the manuscript methods. No matching MPO FMO was found; no MPO-positive gate is drawn. The neutrophil-only effect and treatment/plate blocks limit attribution exclusively to fibroblasts or a target-specific mechanism. The three-comparison correction is a focused exploratory follow-up and is not a correction across all earlier manuscript endpoints.

## File specifications and provenance

Main figure: 7.2 × 5.8 inches; supplementary figure: 7.2 × 7.0 inches. Each is provided as editable SVG, 300 dpi PNG and 600 dpi LZW-compressed TIFF. The two-page PDF contains the main and supplementary figures with embedded fonts. Plots were generated from FCS-derived compensated events, not screenshots from FlowJo. No smoothing was applied. All events contribute to each histogram and sample statistic. Event counts, representative sample identities, numerical source data and the existing exact-test results are supplied alongside the figures. No statistical results were changed during figure preparation.
'''
(OUT/'Figure_legends_and_methods.md').write_text(legend,encoding='utf-8')
for name in ['MPO_sample_measurements.csv','MPO_group_summary.csv','MPO_blockade_comparisons.csv','verification.json']:shutil.copy2(SRC/name,OUT/name)
shutil.copy2(__file__,OUT/Path(__file__).name)
manifest=[{'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='manifest.json']
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
with zipfile.ZipFile(ROOT/'outputs/MPO_manuscript_figure_package.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.iterdir():
        if p.is_file():z.write(p,'MPO_manuscript/'+p.name)
print(json.dumps(selection,indent=2));print('Completed manuscript figure package')
