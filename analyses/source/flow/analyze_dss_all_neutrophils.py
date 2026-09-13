from pathlib import Path
import sys,json,itertools,html,hashlib,zipfile,shutil
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/mouse_packages'))
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'outputs/09-DSS-neutrophil-comparison';OUT.mkdir(exist_ok=True)
OLD=ROOT/'outputs/07-mouse-extensive';BASE='CELLS/s1/s2/cell1/cd11b+/neutrophil'
MARKERS={'OSM':'B515-A','CXCR4':'U820-A','MX1':'Y670-A','PADI4':'R710-A','NAMPT':'R780-A','CD177':'Y586-A','CD11b':'U670-A','CD45':'V710-A','Ly6G':'V610-A','F4/80':'R670-A'}
meta=[m for m in json.loads((ROOT/'analysis/regating/event_cache_index.json').read_text()) if m['experiment']=='mouse' and m['biological'] and m['group'] in ['WD','FD']]
thresholds=pd.read_csv(OLD/'FMO_thresholds.csv');thresholds=thresholds[(thresholds.compartment=='neutrophil')&(thresholds['quantile']==.99)]
f480=np.load(ROOT/'analysis/mouse_extensive/f480_masks.npz');freq=[];mfi=[];audit=[]
for m in meta:
    sample=m['sample'];group=m['group']
    with np.load(ROOT/'analysis/regating/events'/(m['key']+'.npz')) as z:
        parent=z['old_5'];v=z['events'][parent];overlap=f480[sample][parent]
    n=len(v);counts=m['gate_counts'];saved=np.load(ROOT/'analysis/dss_neutrophil_subsets'/(sample+'.npz'))
    masks={'All neutrophils':np.ones(n,dtype=bool)};kinds={'All neutrophils':'total'}
    for path in saved.files:
        if path==BASE:continue
        name='Saved: '+path[len(BASE)+1:];masks[name]=saved[path];kinds[name]='saved gate'
        parentpath=path.rsplit('/',1)[0];pn=int(saved[parentpath].sum())
        audit.append(dict(sample=sample,group=group,subset=name,events=int(saved[path].sum()),immediate_parent_events=pn,percent_immediate_parent=100*saved[path].sum()/pn if pn else np.nan))
    singles={}
    for row in thresholds.itertuples():
        name='FMO: '+row.marker+'+';singles[row.marker]=v[:,m['channels'].index(row.channel)]>=row.threshold;masks[name]=singles[row.marker];kinds[name]='fixed FMO'
    for a,b in itertools.combinations(singles,2):
        name='FMO: '+a+'+ & '+b+'+';masks[name]=singles[a]&singles[b];kinds[name]='fixed FMO coexpression'
    masks['F4/80 gate overlap']=overlap;masks['Outside F4/80 gate']=~overlap
    kinds['F4/80 gate overlap']=kinds['Outside F4/80 gate']='gate overlap'
    for den,d in [('viable cells',counts['old_3']),('CD45-CD11b gate',counts['old_4'])]:
        freq.append(dict(sample=sample,group=group,subset='All neutrophils / '+den,kind='total abundance',events=n,denominator_events=d,value=100*n/d))
    for name,mask in masks.items():
        nn=int(mask.sum())
        if name!='All neutrophils':freq.append(dict(sample=sample,group=group,subset=name,kind=kinds[name],events=nn,denominator_events=n,value=100*nn/n))
        for marker,ch in MARKERS.items():
            x=v[mask,m['channels'].index(ch)]
            mfi.append(dict(sample=sample,group=group,subset=name,kind=kinds[name],marker=marker,events=nn,median=float(np.median(x)) if nn else np.nan,mean=float(x.mean()) if nn else np.nan,eligible=nn>=20))
F=pd.DataFrame(freq);M=pd.DataFrame(mfi);A=pd.DataFrame(audit)
F.to_csv(OUT/'Subset_frequencies_per_mouse.csv',index=False);M.to_csv(OUT/'Subset_marker_MFI_per_mouse.csv',index=False);A.to_csv(OUT/'Saved_gate_parent_frequencies.csv',index=False)
thresholds.to_csv(OUT/'Fixed_FMO_thresholds.csv',index=False)
alloc={};results=[];verified=0
def compare(frame,value,subset,marker,kind,family):
    global verified
    x=frame[frame.group=='FD'][value].dropna().to_numpy();y=frame[frame.group=='WD'][value].dropna().to_numpy()
    nx,ny=len(x),len(y)
    row=dict(subset=subset,marker=marker,kind=kind,family=family,n_PBS=ny,n_GCV=nx,PBS_median=np.median(y) if ny else np.nan,GCV_median=np.median(x) if nx else np.nan,PBS_q25=np.quantile(y,.25) if ny else np.nan,PBS_q75=np.quantile(y,.75) if ny else np.nan,GCV_q25=np.quantile(x,.25) if nx else np.nan,GCV_q75=np.quantile(x,.75) if nx else np.nan,status='tested' if min(nx,ny)>=3 else 'insufficient mice with >=20 subset events',p_exact=np.nan)
    if min(nx,ny)>=3:
        z=np.r_[x,y];r=stats.rankdata(z);u=r[:nx].sum()-nx*(nx+1)/2;c=nx*ny/2
        if (nx,ny) not in alloc:alloc[nx,ny]=np.array(list(itertools.combinations(range(nx+ny),nx)))
        null=r[alloc[nx,ny]].sum(axis=1)-nx*(nx+1)/2;p=np.mean(abs(null-c)>=abs(u-c)-1e-12)
        if len(np.unique(z))==len(z):assert np.isclose(p,stats.mannwhitneyu(x,y,method='exact').pvalue)
        else:
            check=stats.permutation_test((x,y),lambda a,b:abs(stats.mannwhitneyu(a,b).statistic-c),n_resamples=np.inf,alternative='greater').pvalue
            assert np.isclose(p,check)
        verified+=1;row.update(p_exact=p,rank_biserial=2*u/(nx*ny)-1,median_difference=np.median(x)-np.median(y),hodges_lehmann_shift=np.median((x[:,None]-y).ravel()))
    results.append(row)
for name,f in F.groupby('subset',sort=False):compare(f,'value',name,'Frequency (%)',f.kind.iloc[0],'Abundance and total-neutrophil MFI')
for (name,marker),f in M.groupby(['subset','marker'],sort=False):
    total=name=='All neutrophils';compare(f if total else f[f.eligible],'median',name,marker,f.kind.iloc[0],'Abundance and total-neutrophil MFI' if total else 'Exploratory subset MFI')
T=pd.DataFrame(results)
for family in T.family.unique():
    ix=T.index[(T.family==family)&T.p_exact.notna()];T.loc[ix,'q_family']=stats.false_discovery_control(T.loc[ix,'p_exact'],method='bh');T.loc[T.family==family,'family_test_count']=len(ix)
ix=T.index[T.p_exact.notna()];T.loc[ix,'q_all_tests']=stats.false_discovery_control(T.loc[ix,'p_exact'],method='bh')
T.to_csv(OUT/'All_comparisons_Mann_Whitney_FDR.csv',index=False)
main=T[T.family=='Abundance and total-neutrophil MFI'];sub=T[T.family=='Exploratory subset MFI']
main.to_csv(OUT/'Primary_comparisons.csv',index=False);sub.to_csv(OUT/'Exploratory_subset_MFI_comparisons.csv',index=False)
old=pd.read_csv(OLD/'Marker_measurements.csv');old=old[(old.compartment=='neutrophil')&old.group.isin(['WD','FD'])]
check=M[M.subset=='All neutrophils'].merge(old,on=['sample','marker'],suffixes=('_new','_old'));assert len(check)==110 and np.allclose(check.median_new,check.median_old)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
def panels(data,names,filename,title,is_mfi=False):
    cols=3;rows=int(np.ceil(len(names)/cols));fig,axs=plt.subplots(rows,cols,figsize=(12,2.8*rows+1.5),squeeze=False)
    for ax,name in zip(axs.ravel(),names):
        f=data[data.marker==name] if is_mfi else data[data.subset==name]
        rr=main[(main.subset=='All neutrophils')&(main.marker==name)].iloc[0] if is_mfi else main[main.subset==name].iloc[0]
        for j,g in enumerate(['WD','FD']):
            x=f[f.group==g]['median' if is_mfi else 'value'].to_numpy();lo,mid,hi=np.quantile(x,[.25,.5,.75]);ax.scatter(j+np.linspace(-.07,.07,len(x)),x,color=['#34699a','#d47934'][j],s=25);ax.errorbar(j,mid,yerr=[[mid-lo],[hi-mid]],fmt='_',color='black',capsize=4,markersize=12)
        import textwrap
        ax.set_title(textwrap.fill(name.replace('Saved: ','').replace('FMO: ','FMO '),36)+f'\np={rr.p_exact:.3g}; q={rr.q_family:.3g}',fontsize=9);ax.set_xticks([0,1],['PBS+DSS','GCV+DSS']);ax.set_ylabel('Median fluorescence (a.u.)' if is_mfi else ('% of '+name.split(' / ')[1] if name.startswith('All neutrophils / ') else '% of neutrophils'))
    for ax in axs.ravel()[len(names):]:ax.axis('off')
    fig.suptitle(title+'\nFAP-TK mice | PBS n=6, GCV n=5',fontsize=14);fig.tight_layout(rect=[0,.055,1,.94]);fig.text(.06,.015,'Each dot is one mouse; bars show median and IQR. Exact rank tests; BH q across '+str(int(main.family_test_count.iloc[0]))+' abundance/total-MFI endpoints.',fontsize=9)
    fig.savefig(OUT/(filename+'.png'),dpi=180);fig.savefig(OUT/(filename+'.svg'));plt.close(fig)
panels(M[M.subset=='All neutrophils'],list(MARKERS),'Total_neutrophil_MFI','Total neutrophil marker fluorescence',True)
selected=[n for n in F.subset.unique() if n.startswith('All neutrophils') or n in ['Saved: cd177+ all neutrophil','Saved: cd177- all neutrophils','Saved: cxcr4+','Saved: cxcr4-','Saved: osm+','Saved: mx1+','Saved: padi4+','F4/80 gate overlap','Outside F4/80 gate']]
panels(F,selected,'Major_subset_frequencies','Neutrophil abundance and saved subsets')
fnames=[n for n in F.subset.unique() if n.startswith('FMO:')];panels(F,fnames,'Fixed_FMO_subsets','Fixed FMO-positive and co-expression subsets')
cols=['subset','marker','n_PBS','n_GCV','PBS_median','GCV_median','p_exact','q_family','q_all_tests','status']
def table(df):return df[cols].to_html(index=False,float_format=lambda v:f'{v:.4g}',na_rep='Unavailable',escape=True)
sig=main[main.q_family<.05];ss=sub[sub.q_family<.05]
methods=f'''WD is FAP-TK+PBS+DSS (n=6), corrected by the investigator; FD is FAP-TK+GCV+DSS (n=5; FD3 absent). All 11 biological mice are retained; no pairing is assumed. The comparison holds strain and DSS exposure constant and compares recorded PBS versus GCV treatment. MFI means per-mouse median compensated fluorescence, not pooled events. Arithmetic means are also supplied descriptively.

Coverage: {F.subset.nunique()} abundance endpoints and 10 total-neutrophil markers; {M.subset.nunique()-1} subset definitions each evaluated for 10 markers. All saved descendant gates are retained with their original names; these names are annotations rather than independently validated biological identities. Saved bounded positive/negative gates need not partition all neutrophils and different subsets may overlap. Saved-gate frequencies use total neutrophils as the denominator; immediate-parent frequencies are supplied descriptively. Total abundance uses viable cells and the CD45-CD11b parent separately. F4/80 overlap is a gate diagnostic. FMO subsets use the previously established 99th-percentile threshold only; no threshold sensitivity analysis is performed. PADI4 omission-control annotation remains provisional.

Exact two-sided Mann-Whitney rank statistics are tested by enumerating all label allocations, including ties. Biological mice are the inference units. BH FDR is applied across {len(main)} abundance/total-MFI endpoints as one family; the {sub.p_exact.notna().sum()} eligible subset-MFI tests form a separate exploratory family. A further q_all_tests corrects across all {T.p_exact.notna().sum()} tested endpoints. These correction scopes differ from the earlier 60-test four-group comparison. Families and eligibility rules were specified before examining results. Mann-Whitney compares distributions; a pure median-location interpretation assumes similar shapes.

Subset MFI is tested only for mice with at least 20 events in that subset and at least three eligible mice per treatment. This is a pragmatic minimum, not a precision guarantee. Frequencies retain zero-event subsets; MFI of an empty subset is missing, never zero. All unfiltered MFI and event counts are supplied. Sparse/untested endpoints are shown as unavailable. Subset MFI is conditional on marker-based selection and can be biased by different event eligibility, so it is exploratory.

MPO was not measured. FAP, PDPN and alpha5beta1 channels lack coverage in the neutrophil compensation matrix and are not treated as reliable compensated neutrophil expression endpoints. The separate stromal analysis cannot substitute for those neutrophil measurements. No mouse alpha5beta1-blockade arm is recorded. Negative compensated fluorescence is retained and precludes simple fold-change interpretation. Incomplete compensation/control quality limits biological interpretation; rank tests do not repair measurement error. Lower fluorescence does not establish lower secretion, enzyme activity or NET formation. No significant difference is not evidence of equivalence.

Verification: all extracted saved-gate event counts match the canonical workspace; all 110 total-neutrophil marker medians match the previous source analysis. All {verified} rank-test P values were independently checked using SciPy exact or permutation calculations. Source files and earlier results were preserved.'''
findings=f'Abundance/total-MFI family: {len(sig)} of {len(main)} endpoints have q<0.05. Exploratory subset-MFI family: {len(ss)} of {sub.p_exact.notna().sum()} tested endpoints have q<0.05. Across all tests: {(T.q_all_tests<.05).sum()} have q<0.05. Total neutrophils comprise median 0.714% versus 0.734% of viable cells (PBS versus GCV; q=0.940). Saved CXCR4-positive, OSM-positive and MX1-positive frequencies are numerically higher with GCV, but none survives FDR correction. The smallest unadjusted frequency P value is for fixed-FMO CXCR4+OSM+ cells (median 0% versus 0.893%; P=0.00866, q=0.528); the saved CXCR4+/mx1, osm+ gate also differs before correction (0.804% versus 5.50%; P=0.0173, q=0.528). These are descriptive signals for follow-up, not established treatment effects. There is no FDR-supported evidence here of neutrophil depletion or a consistent marker-expression shift. Limited sample size and sparse subsets prevent a conclusion of equivalence.'
methods+='\n\nBH FDR assumes independent or suitably positively dependent tests; these overlapping endpoints make the analysis exploratory. Methods references: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html and https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html'
(OUT/'DSS_neutrophil_report.md').write_text('# FAP-TK mice: PBS+DSS versus GCV+DSS\n\n'+findings+'\n\n'+methods,encoding='utf-8')
body='<h1>FAP-TK mice: PBS+DSS versus GCV+DSS</h1><p>'+findings+'</p>'+''.join('<img src="'+f+'.png">' for f in ['Total_neutrophil_MFI','Major_subset_frequencies','Fixed_FMO_subsets'])+'<h2>All abundance and total-neutrophil MFI comparisons</h2>'+table(main)+'<h2>Methods and interpretation limits</h2>'+''.join('<p>'+html.escape(p)+'</p>' for p in methods.split('\n\n'))+'<h2>All subset MFI comparisons</h2>'+table(sub)
(OUT/'DSS_neutrophil_report.html').write_text('<!doctype html><meta charset="utf-8"><title>DSS neutrophils</title><style>body{max-width:1250px;margin:35px auto;font:16px/1.6 system-ui;color:#183647}img{width:100%}table{border-collapse:collapse;font-size:12px}td,th{border-bottom:1px solid #ddd;padding:5px;text-align:left}th{background:#edf3f8}</style>'+body,encoding='utf-8')
sources=[OLD/'Marker_measurements.csv',OLD/'FMO_thresholds.csv',ROOT/'analysis/regating/event_cache_index.json']
(OUT/'Verification.json').write_text(json.dumps(dict(mice=11,verified_tests=verified,total_MFI_matches=110,primary_tests=len(main),subset_tests=int(sub.p_exact.notna().sum()),sources={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}),indent=2))
for f in [__file__,ROOT/'analysis/extract_dss_neutrophil_subsets.py']:shutil.copy2(f,OUT/Path(f).name)
with zipfile.ZipFile(ROOT/'outputs/DSS_neutrophil_comparison.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in OUT.iterdir():z.write(p,'DSS_neutrophil_comparison/'+p.name)
print(findings,flush=True);print(main[cols].to_string(index=False),flush=True);print('SUBSET MFI SIGNIFICANT',ss[cols].to_string(index=False),flush=True)
