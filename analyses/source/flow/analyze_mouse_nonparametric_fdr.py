"""Per-mouse neutrophil fluorescence: exact rank tests and explicit FDR families."""
from pathlib import Path
import sys,json,itertools,hashlib,shutil,zipfile,html
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/mouse_packages'))
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'outputs/08-mouse-MFI-FDR';OUT.mkdir(exist_ok=True)
SRC=ROOT/'outputs/07-mouse-extensive/Marker_measurements.csv'
R=pd.read_csv(SRC);R=R[R.compartment=='neutrophil'].copy()
ORDER=['WD','FD','FC','FG'];MARKERS=['OSM','CXCR4','MX1','PADI4','NAMPT','CD177','CD11b','CD45','Ly6G','F4/80'];PHENOTYPE=MARKERS[:6]
GROUPS={'WD':'FAP-TK + PBS + DSS','FD':'FAP-TK + GCV + DSS','FC':'FAP-TK + PBS + water','FG':'FAP-TK + GCV + water'}
pd.DataFrame([{'group':g,'annotation':label,'n':R[R.group==g]['sample'].nunique()} for g,label in GROUPS.items()]).to_csv(OUT/'Mouse_group_annotations.csv',index=False)
PAIRS=list(itertools.combinations(ORDER,2));COLORS=['#34699a','#d47934','#427f66','#9566a4']
assert len(R)==230 and R.groupby('marker').size().eq(23).all()
assert R.groupby('group')['sample'].nunique().to_dict()=={'FC':6,'FD':5,'FG':6,'WD':6}
R.to_csv(OUT/'Neutrophil_sample_MFI.csv',index=False)
rng=np.random.default_rng(20260908);perms=np.argsort(rng.random((99999,23)),axis=1)
sizes=[6,5,6,6];bounds=np.cumsum([0]+sizes);alloc={};pairrows=[];omni=[];summary=[];checks=0
def adjust(p):
    p=np.asarray(p);ix=np.argsort(p);ordered=p[ix]*len(p)/np.arange(1,len(p)+1)
    corrected=np.minimum(1,np.minimum.accumulate(ordered[::-1])[::-1]);out=np.zeros(len(p));out[ix]=corrected
    assert np.allclose(out,stats.false_discovery_control(p,method='bh'),atol=1e-12,rtol=0)
    return out
for metric in ['median','mean']:
    for marker in MARKERS:
        fr=R[R.marker==marker];groupvalues={g:fr[fr.group==g].sort_values('sample')[metric].to_numpy() for g in ORDER}
        for g,x in groupvalues.items():summary.append(dict(metric=metric,marker=marker,group=g,n=len(x),group_median=np.median(x),group_q25=np.quantile(x,.25),group_q75=np.quantile(x,.75),group_mean=x.mean(),group_sd=x.std(ddof=1)))
        pooled=np.concatenate(list(groupvalues.values()));ranks=stats.rankdata(pooled);_,counts=np.unique(pooled,return_counts=True);tie=1-np.sum(counts**3-counts)/(23**3-23)
        stat=stats.kruskal(*groupvalues.values());calc=(12/(23*24)*sum(ranks[bounds[j]:bounds[j+1]].sum()**2/sizes[j] for j in range(4))-3*24)/tie
        assert np.isclose(calc,stat.statistic,atol=1e-12)
        shuffled=ranks[perms];hsum=np.zeros(len(perms))
        for j,n in enumerate(sizes):hsum+=shuffled[:,bounds[j]:bounds[j+1]].sum(axis=1)**2/n
        null=(12/(23*24)*hsum-3*24)/tie;p=(1+np.sum(null>=calc-1e-12))/100000
        omni.append(dict(metric=metric,marker=marker,H=calc,df=3,p_permutation=p,p_chisquare=stat.pvalue,permutations=99999,mc_standard_error=np.sqrt(p*(1-p)/100000),n=23))
        for ref,trt in PAIRS:
            x=groupvalues[trt];y=groupvalues[ref];nx,ny=len(x),len(y);z=np.r_[x,y];rank=stats.rankdata(z);u=rank[:nx].sum()-nx*(nx+1)/2
            if (nx,ny) not in alloc:alloc[nx,ny]=np.array(list(itertools.combinations(range(nx+ny),nx)))
            uu=rank[alloc[nx,ny]].sum(axis=1)-nx*(nx+1)/2;center=nx*ny/2;tol=1e-12
            p=float(np.mean(abs(uu-center)>=abs(u-center)-tol));has_ties=len(np.unique(z))<len(z)
            if not has_ties:
                independent=stats.mannwhitneyu(x,y,alternative='two-sided',method='exact');assert np.isclose(p,independent.pvalue,atol=1e-12)
            else:
                independent=stats.permutation_test((x,y),lambda a,b:abs(stats.mannwhitneyu(a,b,method='asymptotic').statistic-center),n_resamples=np.inf,alternative='greater',vectorized=False);assert np.isclose(p,independent.pvalue,atol=1e-12)
            checks+=1
            pairrows.append(dict(metric=metric,marker=marker,role='phenotype' if marker in PHENOTYPE else 'gate-conditioned lineage diagnostic',reference=ref,comparison=trt,n_reference=ny,n_comparison=nx,median_reference=np.median(y),median_comparison=np.median(x),difference_group_medians=np.median(x)-np.median(y),hodges_lehmann_shift=np.median((x[:,None]-y).ravel()),U_comparison=u,rank_biserial=2*u/(nx*ny)-1,p_exact=p,allocations=len(uu),ties=has_ties,contrast_context='within recorded context' if (ref,trt) in [('WD','FD'),('FC','FG')] else 'cross-context; preparation day and design confounding'))
T=pd.DataFrame(pairrows);O=pd.DataFrame(omni);S=pd.DataFrame(summary)
for metric in ['median','mean']:
    ix=T[T.metric==metric].index;assert len(ix)==60;T.loc[ix,'q_BH_60']=adjust(T.loc[ix,'p_exact'])
    ix=O[O.metric==metric].index;O.loc[ix,'q_BH_10']=adjust(O.loc[ix,'p_permutation'])
T.to_csv(OUT/'Pairwise_Mann_Whitney_FDR.csv',index=False);O.to_csv(OUT/'Omnibus_Kruskal_Wallis_FDR.csv',index=False);S.to_csv(OUT/'MFI_group_summary.csv',index=False)
primary=T[T.metric=='median'];focus=primary[((primary.reference=='WD')&(primary.comparison=='FD'))|((primary.reference=='FC')&(primary.comparison=='FG'))]
focus.to_csv(OUT/'Within_context_MFI_FDR.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
caption='WD: FAP-TK + PBS + DSS | FD: FAP-TK + GCV + DSS\nFC: FAP-TK + PBS + water | FG: FAP-TK + GCV + water\nEach point is one mouse (n=6/group; FD n=5). Bars: group median and IQR of per-mouse median fluorescence.\nExact two-sided Mann-Whitney U tests; BH FDR across all 60 comparisons. Negative compensated values retained.'
for name,markers,shape in [('Neutrophil_phenotype_MFI_FDR',PHENOTYPE,(2,3)),('Neutrophil_lineage_MFI_FDR',MARKERS[6:],(2,2))]:
    fig,axs=plt.subplots(*shape,figsize=(11,8));fig.suptitle('Mouse neutrophil median fluorescence | Non-parametric FDR',x=.085,ha='left',fontsize=16,y=.96)
    for ax,k in zip(axs.ravel(),markers):
        for j,g in enumerate(ORDER):
            v=R[(R.marker==k)&(R.group==g)].sort_values('sample')['median'].to_numpy();q25,mid,q75=np.quantile(v,[.25,.5,.75]);ax.scatter(j+np.linspace(-.09,.09,len(v)),v,s=26,color=COLORS[j],zorder=3);ax.errorbar(j,mid,yerr=[[mid-q25],[q75-mid]],fmt='_',color='black',markersize=13,capsize=4)
        ff=focus[focus.marker==k].set_index('reference');ax.set_title(f'{k}\nDSS q={ff.loc["WD","q_BH_60"]:.3g} | Water q={ff.loc["FC","q_BH_60"]:.3g}',fontsize=11);ax.set_xticks(range(4),ORDER);ax.set_ylabel('Per-mouse median fluorescence (a.u.)',fontsize=8)
    fig.subplots_adjust(left=.09,right=.96,top=.84,bottom=.20,hspace=.45,wspace=.40);fig.text(.085,.04,caption,fontsize=8,linespacing=1.5);fig.savefig(OUT/(name+'.png'),dpi=200);fig.savefig(OUT/(name+'.svg'));plt.close(fig)
pairlabels=[trt+' - '+ref for ref,trt in PAIRS];qmat=np.array([[primary[(primary.marker==k)&(primary.reference==ref)&(primary.comparison==trt)].iloc[0].q_BH_60 for ref,trt in PAIRS] for k in MARKERS])
fig,ax=plt.subplots(figsize=(10,7));im=ax.imshow(-np.log10(qmat),cmap='Blues',vmin=0,vmax=max(2,float(-np.log10(qmat).min())),aspect='auto')
for i,k in enumerate(MARKERS):
    for j in range(6):ax.text(j,i,f'{qmat[i,j]:.3f}'+('*' if qmat[i,j]<.05 else ''),ha='center',va='center',fontsize=9)
ax.set_xticks(range(6),pairlabels,rotation=30,ha='right');ax.set_yticks(range(10),MARKERS);ax.set_title('All model comparisons | BH-adjusted Mann-Whitney P values',pad=20);fig.colorbar(im,ax=ax,label='-log10(FDR-adjusted P)');fig.text(.1,.025,'* q < 0.05. Cross-context comparisons mix DSS status with preparation day; they are descriptive associations.',fontsize=8);fig.subplots_adjust(left=.13,right=.96,bottom=.20,top=.90);fig.savefig(OUT/'All_model_MFI_FDR.png',dpi=200);fig.savefig(OUT/'All_model_MFI_FDR.svg');plt.close(fig)
def mdtable(headers,rows):return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)
table=[]
for k in MARKERS:
    ff=focus[focus.marker==k].set_index('reference');table.append([k,f'{ff.loc["WD","q_BH_60"]:.4f}',f'{ff.loc["FC","q_BH_60"]:.4f}'])
text='''# Mouse neutrophil fluorescence: non-parametric tests with FDR

All 23 mice were retained. Primary MFI here means per-mouse median compensated fluorescence in the original neutrophil parent; arithmetic mean fluorescence is supplied as a separate secondary analysis. Group summaries are medians and IQRs of mouse-level measurements, not pooled events. No marker-positive gate or positivity-threshold sensitivity is used. Ten markers are measured: OSM, CXCR4, MX1, PADI4, NAMPT, CD177, CD11b, CD45, Ly6G and F4/80. MPO was not recorded. Alpha5beta1 blockade is not a recorded mouse treatment.

## Tests and correction scope

Pairwise comparisons use exact two-sided Mann-Whitney U rank tests, enumerating all allocations; tied values, if present, use the exact label-permutation distribution of centered U. All six group pairs are tested for all ten markers (60 tests), with Benjamini-Hochberg FDR across the complete 60-test family. The two within-context contrasts shown below retain their adjustment from that full family; they are not readjusted after selection. A separate 60-test family is used for arithmetic-mean MFI. No pairwise comparison is selected based on a significant omnibus result.

The four-group omnibus test uses the tie-corrected Kruskal-Wallis H statistic, with 99,999 Monte Carlo label permutations (seed 20260908; plus-one correction) and BH correction across ten markers. Conventional chi-square P values and Monte Carlo standard errors are also supplied. The mean-MFI omnibus tests form a separate ten-test family. Omnibus significance indicates a difference somewhere among groups, not specifically an effect of GCV during DSS.

Mann-Whitney tests compare distributions; interpreting them purely as tests of population medians requires comparable distribution shapes. The biological mouse is the independent unit. No pairing is assumed. BH has its usual independence/positive-dependence conditions; shared biological signals and experimental design limit confirmatory interpretation. This requested rank/FDR analysis differs from the earlier mean-difference/Holm analysis and is not chosen to maximize significance.

## Within-context comparisons: primary median MFI

'''+mdtable(['Marker','FD vs WD: BH q (60 tests)','FG vs FC: BH q (60 tests)'],table)+'''

WD = FAP-TK+PBS+DSS (n=6); FD = FAP-TK+GCV+DSS (n=5); FC = FAP-TK+PBS+water (n=6); FG = FAP-TK+GCV+water (n=6). WD annotation was corrected by the investigator on 2026-09-07. FD3 is absent; WD3 retains the corrected identity. The DSS comparison now represents GCV versus PBS within FAP-TK mice receiving DSS. No marker reaches q < 0.05 in that comparison. Cross-context comparisons mix DSS status with preparation day, so their significance does not establish a DSS or blockade mechanism.

CD11b, CD45 and Ly6G distributions are conditional on selection by lineage gates. F4/80 and neutrophil gates overlap, so F4/80 signal does not prove macrophage identity. Negative compensated values and incomplete compensation coverage remain material limitations; rank-based tests do not correct measurement error. Lower antibody fluorescence does not by itself establish lower secretion, enzyme activity or cell recruitment. All source measurements and previous analyses are unchanged.

## Verification and files

All 120 pairwise P values (median and mean analyses) were checked against independent SciPy calculations. Custom BH calculations match scipy.stats.false_discovery_control. Kruskal-Wallis observed statistics match SciPy. Exact source data, group median/IQR summaries, omnibus tests, all pairwise tests, rank-biserial effect sizes and Hodges-Lehmann pairwise shifts are supplied. The source measurement-table SHA-256 and software versions are recorded in verification.json.

Methods: [Mann-Whitney U documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html) and [BH FDR documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html).
'''
(OUT/'MFI_FDR_report.md').write_text(text,encoding='utf-8')
body='<h1>Mouse neutrophil MFI | Non-parametric FDR</h1><p>Primary: median fluorescence per mouse. Exact Mann-Whitney tests; BH correction across 60 pairwise comparisons.</p>'+pd.DataFrame(table,columns=['Marker','DSS: FD vs WD q','Water: FG vs FC q']).to_html(index=False)+''.join('<img src="'+f+'.png" style="width:100%">' for f in ['Neutrophil_phenotype_MFI_FDR','Neutrophil_lineage_MFI_FDR','All_model_MFI_FDR'])+'<h2>Full methods</h2><pre style="white-space:pre-wrap;font:15px/1.6 system-ui">'+html.escape(text)+'</pre><h2>Omnibus results</h2>'+O[O.metric=='median'].to_html(index=False,float_format=lambda x:f'{x:.5g}')
(OUT/'MFI_FDR_report.html').write_text('<!doctype html><meta charset="utf-8"><title>Mouse MFI FDR</title><body style="max-width:1100px;margin:40px auto;font:16px/1.6 system-ui;color:#183647">'+body,encoding='utf-8')
verification=dict(source=str(SRC.relative_to(ROOT)),source_sha256=hashlib.sha256(SRC.read_bytes()).hexdigest(),biological_mice=23,primary_marker_measurements=230,pairwise_tests_checked=checks,primary_pairwise_family=60,primary_omnibus_family=10,tied_pairwise_tests=int(T.ties.sum()),all_passed=True,versions=dict(numpy=np.__version__,pandas=pd.__version__,scipy=__import__('scipy').__version__))
(OUT/'verification.json').write_text(json.dumps(verification,indent=2));shutil.copy2(__file__,OUT/Path(__file__).name)
manifest=[dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='manifest.json'];(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2))
with zipfile.ZipFile(ROOT/'outputs/Mouse_neutrophil_MFI_FDR.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.iterdir():
        if p.is_file():z.write(p,'Mouse_MFI_FDR/'+p.name)
print('WITHIN CONTEXT',focus[['marker','reference','comparison','median_reference','median_comparison','p_exact','q_BH_60']].round(5).to_string(index=False),flush=True)
print('OMNIBUS',O[O.metric=='median'][['marker','H','p_permutation','q_BH_10']].round(6).to_string(index=False),flush=True)
print('SIGNIFICANT ALL PAIRS',primary[primary.q_BH_60<.05][['marker','reference','comparison','p_exact','q_BH_60']].round(6).to_string(index=False),flush=True)
print(verification,flush=True)

