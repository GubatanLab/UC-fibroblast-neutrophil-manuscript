from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/mouse_packages'))
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
OUT=ROOT/'outputs/07-mouse-extensive';WORK=ROOT/'analysis/mouse_extensive'
R=pd.read_csv(OUT/'Marker_measurements.csv');F=pd.read_csv(OUT/'FMO_frequencies.csv');P=pd.read_csv(OUT/'Population_frequencies.csv');T=pd.read_csv(OUT/'Statistical_comparisons.csv');Q=pd.read_csv(OUT/'Sample_QC.csv');J=pd.read_csv(OUT/'Joint_frequencies.csv');A=pd.read_csv(OUT/'Mouse_level_associations.csv');H=pd.read_csv(OUT/'FMO_thresholds.csv')
ORDER=['WD','FD','FC','FG'];MARKERS=['CXCR4','OSM','PADI4','MX1','NAMPT','CD177'];COLORS={'WD':'#34699a','FD':'#d47934','FC':'#427f66','FG':'#9566a4'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
manifest=[];representatives=[]
footer='WD: WT + DSS | FD: FAP-TK + GCV + DSS | FC: FAP-TK + PBS + water | FG: FAP-TK + GCV + water\nPoints: individual mice (n=6/group; FD n=5); bars: mean +/- SD. P: exact unpaired, Holm within stated family.\nExploratory: no main comparison passes correction across all 52 primary tests. Immune and stromal viable denominators differ.'
def setup(title,nrow=2,ncol=3):
    fig,axes=plt.subplots(nrow,ncol,figsize=(11,8.2),squeeze=False)
    fig.suptitle(title,x=.08,ha='left',y=.975,fontsize=18)
    fig.subplots_adjust(left=.085,right=.96,top=.87,bottom=.19,hspace=.50,wspace=.43)
    return fig,axes.ravel()
def dot(ax,f,col,title,ylabel,family=None,endpoint=None):
    for i,g in enumerate(ORDER):
        rr=f[f.group==g].sort_values('sample');v=rr[col].to_numpy();x=i+np.linspace(-.09,.09,len(v))
        ax.scatter(x,v,color=COLORS[g],s=23,zorder=3)
        if len(v):ax.errorbar(i,v.mean(),yerr=v.std(ddof=1),fmt='_',color='black',capsize=3,markersize=11)
    ax.set_xticks(range(4),ORDER);ax.set_ylabel(ylabel,fontsize=8)
    if family is not None:
        tt=T[(T.family==family)&(T.endpoint==endpoint)].set_index('reference')
        title+='\n'+ ' | '.join(f'{trt}-{ref}: P={tt.loc[ref,"p_holm_family"]:.3g}' for ref,trt in [('WD','FD'),('FC','FG')])
    ax.set_title(title,fontsize=10,pad=10);ax.tick_params(labelsize=8)
def save(fig,name,caption,pdf):
    fig.text(.08,.04,footer,fontsize=8,linespacing=1.45)
    pdf.savefig(fig);fig.savefig(OUT/(name+'.png'),dpi=200);fig.savefig(OUT/(name+'.svg'));plt.close(fig)
    manifest.append(dict(name=name,caption=caption))
with PdfPages(OUT/'Extensive_mouse_figures.pdf') as pdf:
    fig,axs=setup('01 | Immune population composition')
    eps=list(P.endpoint.unique())[:6]
    for ax,ep in zip(axs,eps):dot(ax,P[P.endpoint==ep],'value',ep,'Percent of stated parent',None if ep=='Viable / singlets' else 'Immune composition (10)',ep)
    save(fig,'Figure_01_immune_composition','Corrected saved immune gates. F4/80 gate is myeloid-enriched, not an exclusive macrophage definition. Acquired events are not absolute tissue counts. Viability is descriptive.',pdf)
    for comp,family,num,title,markers in [('neutrophil','Neutrophil phenotype medians (12)',2,'Neutrophil phenotype',MARKERS),('F4/80 myeloid','F4/80 myeloid phenotype medians (12)',3,'F4/80-gated myeloid phenotype',MARKERS)]:
        fig,axs=setup(f'{num:02d} | {title}')
        for ax,k in zip(axs,markers):dot(ax,R[(R.compartment==comp)&(R.marker==k)],'median',k,'Sample median (compensated a.u.)',family,k)
        save(fig,f'Figure_{num:02d}_'+('neutrophil' if comp=='neutrophil' else 'myeloid')+'_phenotype','All retained cells in the indicated parent contribute to each mouse median. Negative compensated values are retained. Neutrophil family is primary; F4/80 phenotype is a separate exploratory secondary family.',pdf)
    fig,axs=setup('04 | Neutrophil FMO-based frequencies')
    for ax,k in zip(axs,MARKERS[:-1]):
        f=F[(F.compartment=='neutrophil')&(F.marker==k)&(F['quantile']==.99)]
        dot(ax,f,'percent',k,'Above FMO threshold (%)','neutrophil FMO q=0.99 (10)',k+' / neutrophils')
    ax=axs[-1];ax.axis('off');hh=H[(H.compartment=='neutrophil')&(H['quantile']==.99)]
    text='FMO neutrophil counts\n\n'+'\n'.join(f'{r.marker}: {r.control_n} events' for r in hh.itertuples())+'\n\nCD177 FMO is absent.\nR710 omission identity is provisional.\nSparse control tails limit precision.\n98th/99.5th-percentile sensitivity\nis supplied in the data tables.'
    ax.text(0,.95,text,va='top',fontsize=10,linespacing=1.6)
    save(fig,'Figure_04_neutrophil_FMO','Single-marker thresholds use the conservative 99th percentile of available FMOs within retained neutrophils. These are review thresholds, not validated biological positivity.',pdf)
    fig,axs=setup('05 | Stromal candidates and the FAP depletion question')
    for ax,k in zip(axs[:4],['FAPa','PDPN','a5b1','NAMPT']):dot(ax,R[(R.compartment=='stromal')&(R.marker==k)],'median',k,'Sample median (compensated a.u.)','Stromal candidate medians (8)',k)
    for ax,den in zip(axs[4:],['CD45-negative candidates','stromal viable']):dot(ax,F[(F.compartment=='stromal')&(F.marker=='FAPa')&(F['quantile']==.99)&(F.denominator==den)],'percent','FAPa / '+den,'Above FMO threshold (%)','stromal FMO q=0.99 (12)','FAPa / '+den)
    save(fig,'Figure_05_stromal_depletion','Expression in the existing live CD45-negative candidate gate; CD31/EpCAM exclusion is unavailable. FAPa positivity uses its FMO and a shared threshold. Candidate and viable denominators are explicitly separated.',pdf)
    fig,axs=setup('06 | Stromal combinations and threshold sensitivity')
    for ax,k in zip(axs[:2],['PDPN','a5b1']):dot(ax,F[(F.compartment=='stromal')&(F.marker==k)&(F['quantile']==.99)&(F.denominator=='CD45-negative candidates')],'percent',k,'Above FMO threshold (%)','stromal FMO q=0.99 (12)',k+' / CD45-negative candidates')
    for ax,ep in zip(axs[2:5],['FAPa & PDPN','FAPa & a5b1','PDPN & a5b1']):dot(ax,J[(J.compartment=='stromal')&(J.endpoint==ep)],'value',ep,'Both above FMO thresholds (%)','stromal joint FMO (6)',ep)
    ax=axs[-1]
    f=F[(F.compartment=='stromal')&(F.marker=='FAPa')&(F.denominator=='CD45-negative candidates')]
    for g in ORDER:
        s=f[f.group==g].groupby('quantile').percent.mean();ax.plot([98,99,99.5],s.values,'o-',label=g,color=COLORS[g])
    ax.set_xlabel('Control percentile used as cutoff');ax.set_ylabel('Mean FAPa above threshold (%)');ax.set_title('FAPa is threshold-sensitive');ax.legend(frameon=False,ncol=2,fontsize=8)
    save(fig,'Figure_06_stromal_combinations','Joint frequencies are exact intersections on one common candidate parent. The FAPa water-group direction changes between 98th and 99th percentile thresholds; none establishes depletion.',pdf)
    fig,axs=setup('07 | Mouse-level patterns, not new cell subtypes',2,2)
    ax=axs[0];N=R[(R.compartment=='neutrophil')&R.marker.isin(MARKERS)].pivot(index='sample',columns='marker',values='median')[MARKERS]
    corr=N.corr(method='spearman');im=ax.imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r');ax.set_xticks(range(6),MARKERS,rotation=45,ha='right');ax.set_yticks(range(6),MARKERS);ax.set_title('Pooled Spearman correlations')
    matrix=np.eye(6)
    for rr in A[A.kind=='neutrophil markers'].itertuples():i,j=MARKERS.index(rr.a),MARKERS.index(rr.b);matrix[i,j]=matrix[j,i]=rr.group_centered_rank_r
    ax=axs[1];ax.imshow(matrix,vmin=-1,vmax=1,cmap='RdBu_r');ax.set_xticks(range(6),MARKERS,rotation=45,ha='right');ax.set_yticks(range(6),MARKERS);ax.set_title('Ranks after group centering')
    fig.colorbar(im,ax=list(axs[:2]),shrink=.75,fraction=.025,pad=.04,label='Correlation',ticks=[-1,0,1])
    scores=pd.read_csv(OUT/'PCA_scores.csv');details=json.loads((OUT/'PCA_details.json').read_text());ax=axs[2]
    for g in ORDER:
        f=scores[scores.group==g];ax.scatter(f.PC1,f.PC2,color=COLORS[g],label=g,s=30)
    ax.set_xlabel(f'PC1 ({details["variance_fraction"][0]*100:.1f}%)');ax.set_ylabel(f'PC2 ({details["variance_fraction"][1]*100:.1f}%)');ax.legend(frameon=False,ncol=2,fontsize=7);ax.set_title('One point per mouse')
    ax=axs[3];L=pd.read_csv(OUT/'PCA_loadings.csv');ax.barh(np.arange(6)-.17,L.PC1,.33,label='PC1',color='#34699a');ax.barh(np.arange(6)+.17,L.PC2,.33,label='PC2',color='#d47934');ax.set_yticks(range(6),L.marker);ax.legend(frameon=False,fontsize=8);ax.set_title('PCA loadings');ax.set_xlabel('Loading (standardized asinh signals)')
    save(fig,'Figure_07_patterns','PCA uses standardized asinh-transformed per-mouse medians (cofactor 150). Correlations can reflect compensation as well as biology. Group-centered rank associations use 9,999 within-group permutations and BH across 21 tests; cross-compartment tests are in the data tables.',pdf)
    fig,axs=setup('08 | Parent populations and sample quality',2,2)
    ax=axs[0];qq=Q.copy();qq['sort']=qq.group.map({g:i for i,g in enumerate(ORDER)});qq=qq.sort_values(['sort','sample']);ax.bar(range(23),qq.neutrophil_events,color=[COLORS[g] for g in qq.group]);ax.axhline(100,color='black',ls='--',lw=.8);ax.set_xticks(range(23),qq['sample'],rotation=90,fontsize=7);ax.set_ylabel('Retained neutrophil events');ax.set_title('FC3 and FC4 have <100 events')
    dot(axs[1],Q,'viable_percent','Immune viability yield','Percent of scatter-selected singlets')
    dot(axs[2],Q,'overlap_percent_immune','Overlap of immune/stromal live gates','Overlap / immune viable (%)')
    c=pd.read_csv(OUT/'Legacy_CD177_gate_audit.csv');c=c.groupby(['sample','group']).percent.sum().reset_index();c['gap']=100-c.percent;dot(axs[3],c,'gap','Legacy CD177 positive/negative gate gap','Neutrophils in neither gate (%)')
    save(fig,'Figure_08_sample_quality','The two workspaces select different event populations from each specimen. Legacy CD177 boxes do not form an exact complementary partition; gap shown here does not identify a biological subset.',pdf)
    fig,axs=setup('09 | Compensation constrains interpretation',2,2)
    audit=pd.read_csv(OUT/'Single_stain_compensation_audit.csv');a=audit.pivot(index='primary',columns='detector',values='residual_in_negative_IQR')
    vals=a.to_numpy().copy()
    for i,k in enumerate(a.index):
        if k in a.columns:vals[i,list(a.columns).index(k)]=np.nan
    ax=axs[0];ax.imshow(np.clip(vals,-10,10),vmin=-10,vmax=10,cmap='RdBu_r');ax.set_xticks(range(len(a.columns)),a.columns,rotation=90,fontsize=6);ax.set_yticks(range(len(a)),a.index,fontsize=6);ax.set_title('Immune residuals (clipped +/-10 IQR)')
    audit2=pd.read_csv(OUT/'Stromal_single_stain_audit.csv');a=audit2.pivot(index='primary',columns='detector',values='residual_in_negative_IQR');vals=a.to_numpy().copy()
    for i,k in enumerate(a.index):
        if k in a.columns:vals[i,list(a.columns).index(k)]=np.nan
    ax=axs[1];ax.imshow(np.clip(vals,-10,10),vmin=-10,vmax=10,cmap='RdBu_r',aspect='auto');ax.set_xticks(range(len(a.columns)),a.columns,rotation=90,fontsize=6);ax.set_yticks(range(len(a)),a.index,fontsize=6);ax.set_title('Stromal residuals (clipped +/-10 IQR)')
    fig.colorbar(ax.images[0],ax=list(axs[:2]),shrink=.75,fraction=.025,pad=.04,label='Residual / negative IQR',ticks=[-10,0,10])
    ax=axs[2];d=pd.read_csv(OUT/'Compensation_direction_diagnostic.csv');a=d.pivot(index='marker',columns='reference',values='direction_agrees').reindex(MARKERS);ax.imshow(a.astype(int),cmap='RdYlGn',vmin=0,vmax=1,aspect='auto');ax.set_xticks([0,1],['FG-FC','FD-WD']);ax.set_yticks(range(6),MARKERS);ax.set_title('Raw vs compensated direction');
    for i in range(6):
        for j in range(2):ax.text(j,i,'Same' if a.iloc[i,j] else 'Reversed',ha='center',va='center',fontsize=9)
    ax=axs[3];ax.axis('off');ax.text(0,1,'Diagnostic, not a replacement correction\n\nImmune matrix: 11 of 14 fluorescence channels.\nEffective stromal matrix: 6 of 14.\nSome single stains have poor separation.\nRaw signals include spillover and are not biology.\n\nRed/blue heatmap cells show off-target residuals.\nMode fitting is exploratory; weak controls can\nmake residual estimates unreliable.\n\nA full, validated matrix and matched controls\nare needed for mechanistic conclusions.',va='top',linespacing=1.5,fontsize=9)
    save(fig,'Figure_09_compensation','Off-diagonal residual is positive-minus-negative control median in units of negative-control IQR using exploratory two-mode fits. Diagonals are masked. Sign reversals are QC diagnostics on identical retained neutrophils, not evidence that raw fluorescence is preferable.',pdf)
    hh=json.loads((WORK/'histograms.json').read_text())
    for num,ref,trt,label in [(10,'FC','FG','Water'),(11,'WD','FD','DSS')]:
        fig,axs=setup(f'{num:02d} | Neutrophil distributions: {label} context')
        for ax,k in zip(axs,MARKERS):
            for g in [ref,trt]:
                rr=R[(R.compartment=='neutrophil')&(R.marker==k)&(R.group==g)].copy();rr['distance']=abs(rr['median']-rr['median'].median());rep=rr.sort_values(['distance','sample']).iloc[0]['sample']
                representatives.append(dict(marker=k,group=g,sample=rep,rule='Closest sample median to group median; ID breaks ties'))
                hist=[h for h in hh if h['compartment']=='neutrophil' and h['marker']==k and h['group']==g]
                for h in hist:
                    c=np.array(h['counts']);ax.stairs(c/c.max()*100,h['edges'],color=COLORS[g],lw=1.8 if h['sample']==rep else .7,alpha=1 if h['sample']==rep else .2,label=f'{g}: {rep}' if h['sample']==rep else None)
            edges=hist[0]['edges'];ticks=[x for x in [-100000,-1000,0,1000,100000] if edges[0]<=np.arcsinh(x/150)<=edges[-1]]
            ax.set_xticks(np.arcsinh(np.array(ticks)/150),[f'{x:,}' for x in ticks],fontsize=7);ax.set_ylabel('Mode-normalized (%)');ax.set_xlabel('Compensated a.u. (asinh scale)',fontsize=8);ax.set_title(k);ax.legend(frameon=False,fontsize=7)
        save(fig,f'Figure_{num:02d}_'+label.lower()+'_histograms','Unsmoothed, all-event histograms with common bins per marker, asinh cofactor 150, independently peak-normalized. Thick lines are deterministic representatives, thin lines all mice. Representatives are not asserted to be paired. Plot height is not abundance.',pdf)
    fig,axs=setup('12 | Acquisition-time sensitivity examples',2,2)
    worst=Q.sort_values('time_quartile_asinh_median_range_OSM',ascending=False).head(4)
    for ax,r in zip(axs,worst.itertuples()):
        m=json.loads((ROOT/'analysis/regating/events'/f'mouse_{r.sample}.json').read_text())
        with np.load(ROOT/'analysis/regating/events'/(m['key']+'.npz')) as z:v=z['events'][z['old_5']]
        t=v[:,m['channels'].index('Time')];x=np.arcsinh(v[:,m['channels'].index('B515-A')]/150);ax.scatter(t,x,s=8,alpha=.6,color=COLORS[r.group]);ix=np.argsort(t);blocks=np.array_split(ix,4);ax.plot([np.median(t[b]) for b in blocks],[np.median(x[b]) for b in blocks],'ko-',lw=1)
        ax.set_title(f'{r.sample}: {len(x)} neutrophils');ax.set_xlabel('Recorded Time (instrument units)');ax.set_ylabel('asinh(OSM fluorescence / 150)')
    save(fig,'Figure_12_acquisition_time','Four samples with the largest OSM acquisition-quartile median ranges. All retained neutrophils are shown; black lines join quartile medians. Small quartile counts and multimodality can create apparent drift; no acquisition-time exclusion was applied.',pdf)
    fig,axs=setup('13 | Lineage-marker fluorescence diagnostics',2,2)
    for ax,k in zip(axs,['CD45','CD11b','Ly6G','F4/80']):dot(ax,R[(R.compartment=='neutrophil')&(R.marker==k)],'median',k,'Sample median (compensated a.u.)','Lineage marker diagnostic (8)',k)
    save(fig,'Figure_13_lineage_diagnostics','Lineage-marker medians are conditional on the retained neutrophil gates. CD45, CD11b and Ly6G gate truncation limits interpretation as activation or maturity. F4/80 fluorescence alone does not reassign lineage.',pdf)
    fig,axs=setup('14 | The saved myeloid gates overlap',2,2)
    dot(axs[0],Q,'f480_overlap_percent_neutrophils','Neutrophils also inside F4/80 gate','Percent of retained neutrophils')
    for ax,ep in zip(axs[1:3],['F4/80-only / viable','F4/80-only / CD45-CD11b gate']):dot(ax,P[P.endpoint==ep],'value',ep,'Percent of stated parent','Exclusive F4/80 sensitivity (4)',ep)
    ax=axs[3];bottom=np.zeros(4)
    for ep,color,label in [('F4/80-only / CD45-CD11b gate','#e39854','F4/80 gate only'),('Neutrophil-only / CD45-CD11b gate','#5a8db6','Neutrophil gate only'),('Overlap / CD45-CD11b gate','#ad81ac','Both gates'),('Neither / CD45-CD11b gate','#dddddd','Neither gate')]:
        values=P[P.endpoint==ep].groupby('group').value.mean().reindex(ORDER).to_numpy();ax.bar(range(4),values,bottom=bottom,color=color,label=label);bottom+=values
    assert np.allclose(bottom,100);ax.set_xticks(range(4),ORDER);ax.set_ylabel('Mean percent of CD45-CD11b gate');ax.legend(frameon=False,fontsize=7,loc='upper left');ax.set_title('An exclusive partition, including overlap')
    save(fig,'Figure_14_gate_overlap','F4/80 and neutrophil gates are not mutually exclusive. Only means outside the other complete gate, not a validated marker-negative assignment. The exclusive F4/80 comparison is an additional four-test exploratory sensitivity family, not added to the 52 primary tests.',pdf)
(OUT/'Figure_manifest.json').write_text(json.dumps(manifest,indent=2));pd.DataFrame(representatives).to_csv(OUT/'Representative_histograms.csv',index=False)
print('Created',len(manifest),'figure pages',flush=True)

