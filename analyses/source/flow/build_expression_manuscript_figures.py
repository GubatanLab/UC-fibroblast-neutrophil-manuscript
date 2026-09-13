"""Create four-marker manuscript histograms and biological-replicate plots."""
from pathlib import Path
import sys,json,shutil,zipfile,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
OUT=ROOT/'outputs/05-expression-blockade';CACHE=ROOT/'analysis/regating/events'
clean_labels='--manuscript-labels' in sys.argv
if clean_labels:
    previous=OUT;OUT=ROOT/'outputs/06-expression-manuscript';OUT.mkdir(exist_ok=True)
    for name in ['Expression_sample_measurements.csv','Expression_group_summary.csv','Expression_blockade_comparisons.csv','verification.json','Expression_blockade_report.md','marker_assignments.json']:
        shutil.copy2(previous/name,OUT/name)
    metadata=json.loads((OUT/'marker_assignments.json').read_text())
    metadata['figure_label_policy']='User requested removal of provisional figure flags; display FCS marker labels. This is not independent confirmation of the antibody panel.'
    (OUT/'marker_assignments.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
R=pd.read_csv(OUT/'Expression_sample_measurements.csv');T=pd.read_csv(OUT/'Expression_blockade_comparisons.csv');T=T[T.metric=='median_fi']
config=json.loads((OUT/'marker_assignments.json').read_text());MARKERS=config['markers']
ORDER=['N','N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK','N+FK','N+ATN'];PAIRS=[('N+IF','N+IF+ATN'),('N+IF+FK','N+IF+ATN+FK'),('N','N+ATN')]
TITLES=['Inflamed co-culture','Inflamed co-culture + FK-866','Neutrophils alone'];COLORS=['#0072B2','#D55E00'];selection=[]
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.7,'lines.linewidth':1.2,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'ps.fonttype':42})
def export(fig,name):
    fig.savefig(OUT/f'{name}.svg');fig.savefig(OUT/f'{name}.png',dpi=300);fig.savefig(OUT/f'{name}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
def letter(ax,s):ax.text(-.22,1.12,s,transform=ax.transAxes,fontsize=12,fontweight='bold')
with PdfPages(OUT/'Expression_manuscript_figures.pdf') as pdf:
    for marker,ch in MARKERS.items():
        rr=R[R.marker==marker];events={};rep={};short=marker if clean_labels or config['confirmed'] or marker in ['OSM','CXCR4'] else ch.split('-')[0]
        title=marker+' antibody fluorescence' if short==marker else f'{short} fluorescence ({marker} FCS label; provisional identity)'
        for row in rr.itertuples():
            m=json.loads((CACHE/f'co_{row.sample}.json').read_text());z=np.load(CACHE/f'co_{row.sample}.npz');v=z['events'][z['old_3'],m['channels'].index(ch)]
            assert len(v)==row.neutrophil_events and np.isclose(np.median(v),row.median_fi,rtol=0,atol=1e-8);events[row.sample]=v
        for g in ORDER:
            x=rr[rr.group==g].copy();target=x.median_fi.median();x['distance']=abs(x.median_fi-target)
            # Treat numerically equivalent distances as ties; select by stable ID.
            candidates=x[np.isclose(x.distance,x.distance.min(),rtol=1e-12,atol=1e-9)];row=candidates.sort_values('sample').iloc[0];rep[g]=row['sample']
            selection.append(dict(marker=marker,channel=ch,group=g,sample=row['sample'],sample_median=row.median_fi,group_median=target,neutrophil_events=int(row.neutrophil_events)))
        allv=np.concatenate([np.arcsinh(v/150) for v in events.values()]);edges=np.linspace(allv.min()-.02,allv.max()+.02,150);centers=(edges[1:]+edges[:-1])/2
        ticks=np.array([-100000,-10000,-1000,0,1000,10000,100000,1000000]);positions=np.arcsinh(ticks/150);keep=(positions>=edges[0])&(positions<=edges[-1]);ticks=ticks[keep];positions=positions[keep]
        ticklabels=[('−' if v<0 else '')+'10'+str(abs(int(v))).replace('1','',1) if False else ('0' if v==0 else ('−' if v<0 else '')+{3:'10³',4:'10⁴',5:'10⁵',6:'10⁶'}[int(np.log10(abs(v)))]) for v in ticks]
        def hist(ax,v,c,label=None,alpha=1,lw=1.2):
            counts=np.histogram(np.arcsinh(v/150),bins=edges)[0];ax.step(centers,100*counts/counts.max(),where='mid',color=c,label=label,alpha=alpha,lw=lw)
        def histaxes(ax):
            ax.set_xlim(edges[0],edges[-1]);ax.set_ylim(0,108);ax.set_yticks([0,50,100]);ax.set_xticks(positions,ticklabels);ax.set_xlabel(short+' fluorescence (a.u.)');ax.set_ylabel('Normalized count (% of peak)')
        means=rr.groupby('group').median_fi.agg(['mean','std']);low=min(rr.median_fi.min(),(means['mean']-means['std']).min(),0);high=max(rr.median_fi.max(),(means['mean']+means['std']).max());span=high-low
        bottom=low-.08*span;top=high+.28*span
        fig,axes=plt.subplots(2,3,figsize=(7.2,6.1));fig.suptitle(title,x=.1,ha='left',fontsize=11,fontweight='bold',y=.985)
        for i,((ref,trt),subtitle) in enumerate(zip(PAIRS,TITLES)):
            ax=axes[0,i]
            for g,c,label in [(ref,COLORS[0],'No ATN-161'),(trt,COLORS[1],'+ ATN-161')]:hist(ax,events[rep[g]],c,label)
            histaxes(ax);ax.set_title(subtitle,pad=27);ax.legend(loc='upper center',bbox_to_anchor=(.5,1.14),frameon=False,handlelength=1.6,labelspacing=.25);ax.text(.5,-.33,f'{rep[ref]} / {rep[trt]}',transform=ax.transAxes,ha='center',fontsize=7,color='#555555');letter(ax,chr(65+i))
            ax=axes[1,i]
            for j,g in enumerate([ref,trt]):
                values=rr[rr.group==g].sort_values('sample').median_fi.to_numpy();ax.scatter(j+np.linspace(-.09,.09,len(values)),values,s=20,color=COLORS[j],edgecolors='white',linewidths=.4,zorder=3);ax.errorbar(j,values.mean(),yerr=values.std(ddof=1),fmt='_',color='black',markersize=12,capsize=3,lw=.8)
            p=T[(T.marker==marker)&(T.reference==ref)].iloc[0].p_holm_12;y=high+.12*span
            ax.plot([0,0,1,1],[y-.025*span,y,y,y-.025*span],color='black',lw=.7);ax.text(.5,y+.025*span,f'Adjusted P = {p:.4f}',ha='center',fontsize=8)
            ax.set_ylim(bottom,top);ax.set_xlim(-.45,1.45);ax.set_xticks([0,1],['No ATN-161','+ ATN-161']);ax.set_ylabel('Median '+short+' fluorescence (a.u.)');ax.ticklabel_format(axis='y',style='plain',useOffset=False);ax.yaxis.set_major_locator(plt.MaxNLocator(4));letter(ax,chr(68+i))
        fig.subplots_adjust(left=.11,right=.98,top=.83,bottom=.095,wspace=.70,hspace=.70);pdf.savefig(fig);export(fig,'Figure_'+marker);plt.close(fig)
        fig=plt.figure(figsize=(7.2,7.3));fig.suptitle(title+' - all samples',x=.11,ha='left',fontsize=10,fontweight='bold',y=.985)
        gs=fig.add_gridspec(3,3,height_ratios=[1.5,1,1],left=.11,right=.98,top=.93,bottom=.08,hspace=.65,wspace=.65);ax=fig.add_subplot(gs[0,:]);labels=['N','N + NF','N + IF','N + IF\n+ FK','N + IF\n+ ATN','N + IF\n+ ATN + FK','N + FK','N + ATN']
        for i,g in enumerate(ORDER):
            x=rr[rr.group==g].sort_values('sample').median_fi.to_numpy();color=COLORS[1] if 'ATN' in g else COLORS[0];ax.scatter(i+np.linspace(-.09,.09,len(x)),x,s=19,color=color,edgecolors='white',linewidths=.4,zorder=3);ax.errorbar(i,x.mean(),yerr=x.std(ddof=1),fmt='_',color='black',capsize=3,markersize=10,lw=.8);ax.text(i,high+.12*span,f'n = {len(x)}',ha='center',fontsize=7)
        ax.set_ylim(bottom,top);ax.set_xticks(range(8),labels);ax.set_ylabel('Median '+short+' fluorescence (a.u.)');ax.ticklabel_format(axis='y',style='plain',useOffset=False);ax.text(-.065,1.025,'A',transform=ax.transAxes,fontsize=12,fontweight='bold')
        for i,(ref,trt) in enumerate(PAIRS):
            for j,g in enumerate([ref,trt]):
                ax=fig.add_subplot(gs[1+j,i])
                for sample in rr[rr.group==g]['sample']:hist(ax,events[sample],COLORS[j],alpha=.3,lw=.7)
                hist(ax,events[rep[g]],COLORS[j],lw=1.4);histaxes(ax);ax.set_title(g,pad=5);letter(ax,chr(66+j*3+i))
        pdf.savefig(fig);export(fig,'Supplement_'+marker);plt.close(fig)
pd.DataFrame(selection).to_csv(OUT/'Representative_samples.csv',index=False)
legend='''# Four-marker manuscript figure legends

**Main figures: OSM, CXCR4 and the provisionally assigned MX1/PADI4 signals following ATN-161 treatment.** A separate main figure is provided for each detector/marker. (A-C) Representative flow-cytometry histograms from all gated neutrophils in inflamed fibroblast co-culture (A), inflamed co-culture with FK-866 (B), and neutrophils alone (C), without (blue) or with (orange) ATN-161. Sample IDs are printed beneath each panel. Histograms display compensated fluorescence using an inverse-hyperbolic-sine transform with cofactor 150, with ticks in original fluorescence units. Bin edges and fluorescence limits are identical across all panels for a given marker. Each histogram is normalized independently to its maximum bin count (100%); heights do not measure cell abundance. No smoothing was applied. For each marker and group, the representative sample was selected by minimum absolute distance from the group median of sample medians, with numerically equivalent distances resolved by sample ID. Samples across conditions are not asserted to be paired. (D-F) Median fluorescence in each independent biological sample for the corresponding conditions; points show all replicates and bars show mean ± SD, n = 6 per condition. Exact two-sided unpaired label-permutation tests (924 allocations) compare means of sample medians. Figure P values use Holm adjustment across all 12 primary tests: four markers by three ATN-161 contrasts. All neutrophils in the retained scatter/singlet/CD16/CD11b parent are included, without CXCR4 or marker-positivity restrictions. Negative compensated values were retained.

**Supplementary figures: all conditions and all-sample distributions.** (A) Median fluorescence across all eight available conditions, with one point per biological sample and mean ± SD. n = 6 except N + NF (n = 5; P5-5 absent). (B-G) Individual histograms for the six groups in the primary comparisons. Thin transparent lines show every sample and thick lines identify the representatives used in the main figure. Transform, bins and normalization match the corresponding main figure. No new tests are applied to the additional descriptive conditions. N = neutrophils; NF/IF = non-inflamed/inflamed fibroblasts; FK = FK-866; ATN = ATN-161.

## Panel identity and submission notes

Consult marker_assignments.json for detector identities and confirmation status. OSM is R670-A and CXCR4 is V610-A. Where unresolved, B710/B515 are identified by detector on axes and their FCS marker assignments are explicitly provisional in figure titles. Do not adopt definitive MX1/PADI4 claims until the antibody panel is reconciled. Figure headings and this legend must follow the confirmed identities before submission.

These plots describe antibody fluorescence, not biological activity or a validated positive fraction. The parent gates and original compensation matrix were retained. A viability channel is absent from 44/47 files, so the population cannot consistently be called live neutrophils. Matching marker FMOs are unavailable. Treatment/plate blocks and unconfirmed pairing constrain inference; complete antibody clones/conjugates, inhibitor doses, staining conditions, culture duration and experimental blocks in the manuscript methods. Primary median and secondary mean analyses can differ; the results report documents both. The primary family is a focused exploratory follow-up, not a prospectively registered or manuscript-wide correction.

## Files

The eight-page PDF alternates main and supplementary figures for OSM, CXCR4, MX1-labeled signal and PADI4-labeled signal. Main figures are 7.2 × 6.1 inches; supplements are 7.2 × 7.3 inches. Each is also supplied as editable SVG, 300 dpi PNG and 600 dpi LZW-compressed TIFF. PDF fonts are embedded. Source-data CSVs and representative sample identities accompany the figures. All panels were generated from FCS-derived compensated events, not FlowJo screenshots. All acquired events within the neutrophil parent contribute to counts, histograms and sample statistics.
'''
if config['confirmed']:legend=legend.replace('the provisionally assigned MX1/PADI4 signals','MX1 and PADI4 expression')
if clean_labels:
    prefix,notes=legend.split('## Panel identity and submission notes',1)
    notes,files=notes.split('## Files',1)
    prefix=prefix.replace('the provisionally assigned MX1/PADI4 signals','MX1 and PADI4 expression')
    legend=prefix+'Detector assignments follow the FCS annotations: OSM, R670-A; CXCR4, V610-A; MX1, B710-A; PADI4, B515-A.\n\n## Files'+files.replace('MX1-labeled signal and PADI4-labeled signal','MX1 and PADI4')
    note_text='# Analysis notes\n\nFigure labels follow the FCS annotations at the user\'s request. The protocol reverses the B515/B710 marker assignments; removing figure flags does not resolve that source discrepancy. Marker-assignment metadata and the original analysis report retain this provenance.\n\n'+notes.split('\n\nThese plots',1)[1]
    (OUT/'Analysis_notes.md').write_text(note_text,encoding='utf-8')
(OUT/'Figure_legends_and_methods.md').write_text(legend,encoding='utf-8');shutil.copy2(__file__,OUT/Path(__file__).name)
manifest=[{'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='manifest.json'];(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
archive='Expression_manuscript_figures.zip' if clean_labels else 'Expression_blockade_manuscript_package.zip'
with zipfile.ZipFile(ROOT/'outputs'/archive,'w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.iterdir():
        if p.is_file():z.write(p,'Expression_blockade/'+p.name)
print('Created eight-page PDF, eight SVG/PNG/TIFF figures, legends and source tables.')
