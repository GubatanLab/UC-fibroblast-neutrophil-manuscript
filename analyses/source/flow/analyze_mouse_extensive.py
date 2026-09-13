"""Exploratory multi-marker mouse analysis. Biological mouse is the inference unit."""
from pathlib import Path
import sys,json,hashlib,itertools,shutil,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/mouse_packages'))
import numpy as np,pandas as pd
from scipy import stats
OUT=ROOT/'outputs/07-mouse-extensive';OUT.mkdir(exist_ok=True)
WORK=ROOT/'analysis/mouse_extensive';CACHE=ROOT/'analysis/regating/events'
ORDER=['WD','FD','FC','FG'];PAIRS=[('WD','FD'),('FC','FG')]
FUNCTIONAL={'CXCR4':'U820-A','OSM':'B515-A','PADI4':'R710-A','MX1':'Y670-A','NAMPT':'R780-A','CD177':'Y586-A'}
DIAGNOSTIC={'CD45':'V710-A','CD11b':'U670-A','Ly6G':'V610-A','F4/80':'R670-A'}
STROMAL={'FAPa':'U379-A','PDPN':'Y780-A','a5b1':'B710-A','NAMPT':'R780-A'}
meta=json.loads((ROOT/'analysis/regating/event_cache_index.json').read_text());meta=[m for m in meta if m['experiment']=='mouse']
sm=json.loads((WORK/'stromal_index.json').read_text())
inventory=pd.read_csv(ROOT/'analysis/results/raw_fcs_inventory.csv').set_index('file')
bio=sorted([m for m in meta if m['biological']],key=lambda m:(ORDER.index(m['group']),m['sample']))
assert len(bio)==23
myeloid_masks=np.load(WORK/'f480_masks.npz')
controls={};thresholds=[];control_events={}
for compartment,markers,metadata,folder in [('neutrophil',{k:v for k,v in FUNCTIONAL.items() if k!='CD177'},meta,CACHE),('stromal',STROMAL,sm,WORK/'stromal_events')]:
    for marker,ch in markers.items():
        name='samples_FMO '+{'PADI4':'R710','FAPa':'FAPA','a5b1':'A5B1'}.get(marker,marker)
        m=next(m for m in metadata if m['sample']==name)
        with np.load(folder/(m['key']+'.npz')) as z:
            v=z['events'][z['old_5'],m['channels'].index(ch)] if compartment=='neutrophil' else z['events'][:,m['channels'].index(ch)]
        assert len(v)>0 and np.isfinite(v).all()
        control_events[compartment,marker]=v
        for q in [.98,.99,.995]:
            order=float(np.quantile(v,q,method='higher'));cut=order+max(abs(order)*1e-9,1e-6)
            controls[compartment,marker,q]=cut;tail=int((v>=cut).sum())
            thresholds.append(dict(compartment=compartment,marker=marker,channel=ch,quantile=q,threshold=cut,control_n=len(v),tail_n=tail,tail_percent=100*tail/len(v),conditional_tail_upper95_percent=100*stats.beta.ppf(.95,tail+1,len(v)-tail),control_file=m['file'],status='R710 omission provisional' if marker=='PADI4' else 'Review FMO threshold; matrix/parent limitations apply'))
pd.DataFrame(thresholds).to_csv(OUT/'FMO_thresholds.csv',index=False)
rows=[];freq=[];joint=[];pops=[];qcs=[];savedcd=[];distributions={};hashes=[]
for m in bio:
    sample=m['sample'];group=m['group'];digest=hashlib.sha256((ROOT/m['file']).read_bytes()).hexdigest()
    assert digest==inventory.loc[m['file'],'sha256']
    hashes.append(dict(sample=sample,file=m['file'],sha256=digest,cache_sha256=hashlib.sha256((CACHE/(m['key']+'.npz')).read_bytes()).hexdigest()))
    z=np.load(CACHE/(m['key']+'.npz'));ev=z['events'];parent=z['old_5'];v=ev[parent];n=len(v)
    assert np.isfinite(v).all()
    counts=m['gate_counts'];live=counts['old_3'];myeloid=counts['old_4'];singlets=counts['old_2']
    baseline={r['gate']:r for r in m['baseline_comparison']}
    base='CELLS/s1/s2/cell1/cd11b+';f480=baseline[base+'/f480+']['canonical_count']
    mac=myeloid_masks[sample]
    assert int(mac.sum())==f480,(sample,int(mac.sum()),f480)
    for k,ch in FUNCTIONAL.items():
        x=ev[mac,m['channels'].index(ch)];distributions['F4/80 myeloid',sample,k]=x
        rows.append(dict(compartment='F4/80 myeloid',sample=sample,group=group,marker=k,channel=ch,n=len(x),median=np.median(x),mean=x.mean(),q25=np.quantile(x,.25),q75=np.quantile(x,.75),negative_percent=100*np.mean(x<0),raw_reconstructed_median=np.nan,role='Exploratory phenotype in retained F4/80 gate'))
    for endpoint,num,den in [('Neutrophils / viable',n,live),('Neutrophils / CD45-CD11b gate',n,myeloid),('F4/80 gate / viable',f480,live),('F4/80 gate / CD45-CD11b gate',f480,myeloid),('CD45-CD11b gate / viable',myeloid,live),('Viable / singlets',live,singlets)]:
        pops.append(dict(sample=sample,group=group,endpoint=endpoint,numerator=num,denominator=den,value=100*num/den,source='corrected immune workspace'))
    overlap_n=int((mac&parent).sum());exclusive_f=f480-overlap_n;exclusive_n=n-overlap_n;neither=myeloid-(f480+n-overlap_n)
    assert min(exclusive_f,exclusive_n,neither)>=0
    for endpoint,num,den in [('F4/80-only / viable',exclusive_f,live),('F4/80-only / CD45-CD11b gate',exclusive_f,myeloid),('Neutrophil-only / CD45-CD11b gate',exclusive_n,myeloid),('Overlap / CD45-CD11b gate',overlap_n,myeloid),('Neither / CD45-CD11b gate',neither,myeloid)]:
        pops.append(dict(sample=sample,group=group,endpoint=endpoint,numerator=num,denominator=den,value=100*num/den,source='Boolean partition of retained gates; only means outside the other whole gate'))
    for endpoint in ['cd177+ all neutrophil','cd177- all neutrophils']:
        k=base+'/neutrophil/'+endpoint
        savedcd.append(dict(sample=sample,group=group,endpoint=endpoint,count=baseline[k]['canonical_count'],parent=n,percent=100*baseline[k]['canonical_count']/n,status='Legacy bounded gate, no matched CD177 FMO'))
    ids=[m['channels'].index(c) for c in m['matrix_info']['detectors']]
    raw_reconstructed=v.copy();raw_reconstructed[:,ids]=v[:,ids]@np.array(m['matrix'])
    for marker,ch in {**FUNCTIONAL,**DIAGNOSTIC}.items():
        j=m['channels'].index(ch);x=v[:,j];raw=raw_reconstructed[:,j]
        distributions['neutrophil',sample,marker]=x
        rows.append(dict(compartment='neutrophil',sample=sample,group=group,marker=marker,channel=ch,n=n,median=np.median(x),mean=x.mean(),q25=np.quantile(x,.25),q75=np.quantile(x,.75),negative_percent=100*np.mean(x<0),raw_reconstructed_median=np.median(raw),role='phenotype' if marker in FUNCTIONAL else 'lineage diagnostic'))
    for q in [.98,.99,.995]:
        masks={k:v[:,m['channels'].index(ch)]>=controls['neutrophil',k,q] for k,ch in FUNCTIONAL.items() if k!='CD177'}
        for k,mask in masks.items():
            nn=int(mask.sum());lo,hi=stats.beta.ppf([.025,.975],[nn if nn else 1,nn+1],[n-nn+1,n-nn if nn<n else 1])
            freq.append(dict(compartment='neutrophil',sample=sample,group=group,marker=k,quantile=q,positive_events=nn,parent_events=n,denominator='neutrophils',percent=100*nn/n,binomial_ci_low=100*lo if nn else 0,binomial_ci_high=100*hi if nn<n else 100))
        if q==.99:
            for a,b in itertools.combinations(masks,2):
                nn=int((masks[a]&masks[b]).sum());assert nn<=min(masks[a].sum(),masks[b].sum())
                joint.append(dict(compartment='neutrophil',sample=sample,group=group,endpoint=a+' & '+b,positive_events=nn,parent_events=n,value=100*nn/n,sparse=nn<10))
    # Temporal drift is diagnostic: split retained neutrophils by acquisition time quartiles.
    ti=m['channels'].index('Time');order=np.argsort(v[:,ti],kind='stable');blocks=np.array_split(order,4)
    drift={k:float(np.ptp([np.median(np.arcsinh(v[b,m['channels'].index(ch)]/150)) for b in blocks])) for k,ch in FUNCTIONAL.items()}
    qcs.append(dict(sample=sample,group=group,raw_events=m['raw_events'],singlets=singlets,viable_events=live,myeloid_events=myeloid,neutrophil_events=n,low_neutrophil_n=n<100,viable_percent=100*live/singlets,immune_matrix_condition=m['matrix_info']['condition_number'],**{'time_quartile_asinh_median_range_'+k:d for k,d in drift.items()}))
    qcs[-1].update(f480_events=f480,f480_neutrophil_overlap=int((mac&parent).sum()),f480_overlap_percent_neutrophils=100*(mac&parent).sum()/n)
    z.close()
    print('Immune',sample,n,flush=True)
for m in sm:
    if not m['biological']:continue
    sample=m['sample'];group=m['group'];z=np.load(WORK/'stromal_events'/(m['key']+'.npz'));v=z['events'];n=len(v)
    assert m['raw_sha256']==inventory.loc[m['file'],'sha256']
    live=m['gate_counts']['CELLS/S1/s2/live'];singlets=m['gate_counts']['CELLS/S1/s2']
    for endpoint,num,den in [('CD45-negative / stromal viable',n,live),('Stromal viable / singlets',live,singlets)]:
        pops.append(dict(sample=sample,group=group,endpoint=endpoint,numerator=num,denominator=den,value=100*num/den,source='corrected stromal workspace'))
    for k,ch in STROMAL.items():
        assert ch in m['matrix_detectors'];x=v[:,m['channels'].index(ch)]
        distributions['stromal',sample,k]=x
        rows.append(dict(compartment='stromal',sample=sample,group=group,marker=k,channel=ch,n=n,median=np.median(x),mean=x.mean(),q25=np.quantile(x,.25),q75=np.quantile(x,.75),negative_percent=100*np.mean(x<0),raw_reconstructed_median=np.nan,role='CD45-negative candidate expression'))
    for q in [.98,.99,.995]:
        masks={k:v[:,m['channels'].index(ch)]>=controls['stromal',k,q] for k,ch in STROMAL.items()}
        for k,mask in masks.items():
            nn=int(mask.sum())
            for denname,den in [('CD45-negative candidates',n),('stromal viable',live)]:
                freq.append(dict(compartment='stromal',sample=sample,group=group,marker=k,quantile=q,positive_events=nn,parent_events=den,denominator=denname,percent=100*nn/den))
        if q==.99:
            for a,b in itertools.combinations(['FAPa','PDPN','a5b1'],2):
                nn=int((masks[a]&masks[b]).sum());assert nn<=min(masks[a].sum(),masks[b].sum())
                joint.append(dict(compartment='stromal',sample=sample,group=group,endpoint=a+' & '+b,positive_events=nn,parent_events=n,value=100*nn/n,sparse=nn<10))
    # Same acquired file, two different saved viable parents: overlap is QC, not a shared denominator.
    imm=next(b for b in bio if b['sample']==sample)
    with np.load(CACHE/(imm['key']+'.npz')) as iz:immune_live=iz['event_indices'][iz['old_3']]
    overlap=len(np.intersect1d(immune_live,z['live_indices'],assume_unique=True));qc=next(r for r in qcs if r['sample']==sample)
    qc.update(stromal_viable=live,stromal_cd45_negative=n,live_parent_overlap=overlap,overlap_percent_immune=100*overlap/len(immune_live),overlap_percent_stromal=100*overlap/live,stromal_matrix_condition=m['matrix_condition_number'])
    z.close()
R=pd.DataFrame(rows);F=pd.DataFrame(freq);J=pd.DataFrame(joint);P=pd.DataFrame(pops);Q=pd.DataFrame(qcs);C=pd.DataFrame(savedcd)
assert len(R)==23*20 and R.groupby(['compartment','marker']).size().eq(23).all()
old=pd.read_csv(ROOT/'outputs/06-mouse-neutrophils/Mouse_sample_measurements.csv')
check=R[R.compartment=='neutrophil'].merge(old,on=['sample','marker']);assert len(check)==92 and np.allclose(check['median'],check.median_fi,rtol=0,atol=1e-8)
for table,name in [(R,'Marker_measurements'),(F,'FMO_frequencies'),(J,'Joint_frequencies'),(P,'Population_frequencies'),(Q,'Sample_QC'),(C,'Legacy_CD177_gate_audit'),(pd.DataFrame(hashes),'Source_hashes')]:table.to_csv(OUT/(name+'.csv'),index=False)
R.groupby(['compartment','marker','group'],sort=False).agg(n=('sample','size'),mean_median=('median','mean'),sd_median=('median','std'),mean_mean=('mean','mean'),mean_negative_percent=('negative_percent','mean')).reset_index().to_csv(OUT/'Marker_group_summary.csv',index=False)
P.groupby(['endpoint','group'],sort=False).agg(n=('sample','size'),mean=('value','mean'),sd=('value','std')).reset_index().to_csv(OUT/'Population_group_summary.csv',index=False)
F.groupby(['compartment','marker','quantile','denominator','group'],sort=False).agg(n=('sample','size'),mean=('percent','mean'),sd=('percent','std')).reset_index().to_csv(OUT/'FMO_group_summary.csv',index=False)

allocations={};tests=[];independent=0
def exact(x,y,verify=False):
    global independent
    nx,ny=len(x),len(y);pooled=np.r_[x,y];diff=x.mean()-y.mean()
    if (nx,ny) not in allocations:allocations[nx,ny]=np.array(list(itertools.combinations(range(nx+ny),nx)))
    sx=pooled[allocations[nx,ny]].sum(axis=1);d=sx/nx-(pooled.sum()-sx)/ny
    p=np.mean(abs(d)>=abs(diff)-max(1e-9,abs(diff)*1e-12))
    if verify:
        check=stats.permutation_test((x,y),lambda a,b:abs(a.mean()-b.mean()),vectorized=False,n_resamples=np.inf,alternative='greater').pvalue
        assert np.isclose(p,check,atol=1e-12,rtol=0),(p,check);independent+=1
    return float(p),len(d)
def compare(frame,endpoint,value,family,tier='primary',verify=True):
    for ref,trt in PAIRS:
        a=frame[frame.group==trt].sort_values('sample');b=frame[frame.group==ref].sort_values('sample')
        x=a[value].dropna().to_numpy();y=b[value].dropna().to_numpy();nx,ny=len(x),len(y)
        if min(nx,ny)<3:continue
        diff=x.mean()-y.mean();p,nalloc=exact(x,y,verify)
        vx=x.var(ddof=1)/nx;vy=y.var(ddof=1)/ny
        df=(vx+vy)**2/(vx*vx/(nx-1)+vy*vy/(ny-1)) if vx+vy else np.nan
        half=stats.t.ppf(.975,df)*np.sqrt(vx+vy) if vx+vy else 0
        pooledsd=np.sqrt(((nx-1)*x.var(ddof=1)+(ny-1)*y.var(ddof=1))/(nx+ny-2))
        loo=np.r_[[np.delete(x,i).mean()-y.mean() for i in range(nx)],[x.mean()-np.delete(y,i).mean() for i in range(ny)]]
        tests.append(dict(family=family,tier=tier,endpoint=endpoint,metric=value,reference=ref,comparison=trt,n_reference=ny,n_comparison=nx,mean_reference=y.mean(),mean_comparison=x.mean(),difference=diff,ci95_low=diff-half,ci95_high=diff+half,hedges_g=(1-3/(4*(nx+ny)-9))*diff/pooledsd if pooledsd else np.nan,p_exact=p,allocations=nalloc,loo_difference_min=loo.min(),loo_difference_max=loo.max(),loo_same_direction_fraction=np.mean(np.sign(loo)==np.sign(diff))))
for marker in FUNCTIONAL:
    f=R[(R.compartment=='neutrophil')&(R.marker==marker)]
    compare(f,marker,'median','Neutrophil phenotype medians (12)')
    compare(f,marker,'mean','Neutrophil phenotype means (12)','secondary')
    keep=Q.loc[~Q.low_neutrophil_n,'sample'];compare(f[f['sample'].isin(keep)],marker,'median','Neutrophil >=100-event sensitivity (12)','sensitivity')
for marker in DIAGNOSTIC:compare(R[(R.compartment=='neutrophil')&(R.marker==marker)],marker,'median','Lineage marker diagnostic (8)','diagnostic')
for marker in STROMAL:compare(R[(R.compartment=='stromal')&(R.marker==marker)],marker,'median','Stromal candidate medians (8)')
for marker in FUNCTIONAL:compare(R[(R.compartment=='F4/80 myeloid')&(R.marker==marker)],marker,'median','F4/80 myeloid phenotype medians (12)','secondary')
for endpoint in P.endpoint.unique():
    if endpoint in ['F4/80-only / viable','F4/80-only / CD45-CD11b gate']:
        compare(P[P.endpoint==endpoint],endpoint,'value','Exclusive F4/80 sensitivity (4)','secondary');continue
    if endpoint in ['Neutrophil-only / CD45-CD11b gate','Overlap / CD45-CD11b gate','Neither / CD45-CD11b gate']:continue
    if endpoint in ['Viable / singlets','Stromal viable / singlets','CD45-negative / stromal viable']:continue
    compare(P[P.endpoint==endpoint],endpoint,'value','Immune composition (10)')
for (comp,marker,q,den),f in F.groupby(['compartment','marker','quantile','denominator']):
    if comp=='stromal' and marker=='NAMPT':continue
    size=10 if comp=='neutrophil' else 12
    compare(f,marker+' / '+den,'percent',f'{comp} FMO q={q:g} ({size})','primary' if q==.99 else 'sensitivity',verify=q==.99)
for (comp,endpoint),f in J.groupby(['compartment','endpoint']):compare(f,endpoint,'value',f'{comp} joint FMO (20)' if comp=='neutrophil' else 'stromal joint FMO (6)','secondary')
T=pd.DataFrame(tests)
def holm(p):
    ix=np.argsort(p,kind='stable');v=np.minimum(1,np.maximum.accumulate(np.asarray(p)[ix]*np.arange(len(p),0,-1)));out=np.zeros(len(p));out[ix]=v;return out
for fam,ix in T.groupby('family').groups.items():T.loc[ix,'p_holm_family']=holm(T.loc[ix,'p_exact'])
ix=T[T.tier=='primary'].index;assert len(ix)==52;T.loc[ix,'p_holm_all_primary_52']=holm(T.loc[ix,'p_exact'])
T.to_csv(OUT/'Statistical_comparisons.csv',index=False)
T[T.tier=='primary'].to_csv(OUT/'Primary_comparisons.csv',index=False)

# Rank associations between biological mice; stratified permutations retain treatment group.
N=R[(R.compartment=='neutrophil')&R.marker.isin(FUNCTIONAL)].pivot(index='sample',columns='marker',values='median').reindex([m['sample'] for m in bio])
groups=np.array([m['group'] for m in bio]);rng=np.random.default_rng(20260907);perms=np.tile(np.arange(23),(9999,1))
for g in ORDER:
    ix=np.where(groups==g)[0]
    for row in perms:row[ix]=rng.permutation(ix)
def rankres(v):
    r=stats.rankdata(v)
    for g in ORDER:r[groups==g]-=r[groups==g].mean()
    return r
assoc=[]
def correlate(namea,a,nameb,b,kind):
    ra,rb=rankres(a),rankres(b);den=np.linalg.norm(ra)*np.linalg.norm(rb);rho=np.dot(ra,rb)/den if den else np.nan
    null=(rb[perms]@ra)/den if den else np.zeros(9999)
    p=(1+np.sum(abs(null)>=abs(rho)-1e-12))/10000 if den else np.nan
    assoc.append(dict(kind=kind,a=namea,b=nameb,n=23,pooled_spearman=stats.spearmanr(a,b).statistic,group_centered_rank_r=rho,p_stratified_permutation=p,resamples=9999))
for a,b in itertools.combinations(FUNCTIONAL,2):correlate(a,N[a].to_numpy(),b,N[b].to_numpy(),'neutrophil markers')
st=R[R.compartment=='stromal'].pivot(index='sample',columns='marker',values='median').reindex(N.index)
nf=P[P.endpoint=='Neutrophils / viable'].set_index('sample').value.reindex(N.index)
for a,b in [('FAPa','neutrophil_frequency'),('FAPa','OSM'),('FAPa','NAMPT'),('PDPN','neutrophil_frequency'),('PDPN','CXCR4'),('a5b1','CXCR4')]:correlate('stromal '+a,st[a].to_numpy(),'neutrophil '+b,nf.to_numpy() if b=='neutrophil_frequency' else N[b].to_numpy(),'cross compartment')
A=pd.DataFrame(assoc);ix=np.argsort(A.p_stratified_permutation);v=A.p_stratified_permutation.to_numpy()[ix]*len(A)/np.arange(1,len(A)+1);q=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);A.loc[ix,'q_BH_21']=q;A.to_csv(OUT/'Mouse_level_associations.csv',index=False)

# PCA summarizes mice, without treating events as replicates or defining cell subtypes.
mat=np.arcsinh(N[list(FUNCTIONAL)].to_numpy()/150);center=mat.mean(axis=0);scale=mat.std(axis=0,ddof=1);X=(mat-center)/scale
u,s,vt=np.linalg.svd(X,full_matrices=False);scores=u*s;variance=s*s/(s*s).sum()
pd.DataFrame(dict(sample=N.index,group=groups,PC1=scores[:,0],PC2=scores[:,1])).to_csv(OUT/'PCA_scores.csv',index=False)
pd.DataFrame(vt.T,index=list(FUNCTIONAL),columns=[f'PC{i+1}' for i in range(6)]).rename_axis('marker').reset_index().to_csv(OUT/'PCA_loadings.csv',index=False)
(OUT/'PCA_details.json').write_text(json.dumps(dict(variance_fraction=variance.tolist(),markers=list(FUNCTIONAL),asinh_cofactor=150,center=center.tolist(),scale=scale.tolist(),unit='one mouse; descriptive, no clustering inference'),indent=2))

# Compensation-direction diagnostic uses the same retained cells; raw signals are not replacement biology.
compdiag=[]
for k in FUNCTIONAL:
    f=R[(R.compartment=='neutrophil')&(R.marker==k)]
    for ref,trt in PAIRS:
        raw=f[f.group==trt].raw_reconstructed_median.mean()-f[f.group==ref].raw_reconstructed_median.mean()
        comp=f[f.group==trt]['median'].mean()-f[f.group==ref]['median'].mean()
        compdiag.append(dict(marker=k,reference=ref,comparison=trt,compensated_difference=comp,reconstructed_raw_difference=raw,direction_agrees=np.sign(raw)==np.sign(comp)))
pd.DataFrame(compdiag).to_csv(OUT/'Compensation_direction_diagnostic.csv',index=False)
shutil.copy2(WORK/'stromal_corrected_vs_saved.csv',OUT/'Stromal_corrected_vs_saved.csv')
audit=pd.read_csv(ROOT/'analysis/regating/compensation_control_audit.csv');audit[audit.experiment=='mouse'].to_csv(OUT/'Single_stain_compensation_audit.csv',index=False)
for name in ['dss_corrected_sample_key.csv','dss_preserved_exclusions.csv']:shutil.copy2(ROOT/'analysis/results'/name,OUT/name)
jsonout={'user_confirmed_model':'FAP-TK','background':'not confirmed','tissue':'not explicitly confirmed; DSS study context','pairing':'unconfirmed; unpaired assumption','groups':{'WD':'WT + DSS','FD':'FAP-TK + GCV + DSS','FC':'FAP-TK + PBS + water','FG':'FAP-TK + GCV + water'},'verification':{'raw_hashes':23,'prior_marker_medians':92,'new_marker_measurements':len(R),'independent_exact_tests':independent,'primary_tests':52,'all_passed':True},'scientific_versions':{'numpy':np.__version__,'pandas':pd.__version__,'scipy':__import__('scipy').__version__}}
(OUT/'Analysis_design_and_verification.json').write_text(json.dumps(jsonout,indent=2))
shutil.copy2(__file__,OUT/Path(__file__).name)
for fn in ['extract_mouse_stromal.py','extract_mouse_myeloid_masks.py','audit_mouse_extensive.py']:shutil.copy2(ROOT/'analysis'/fn,OUT/fn)
# Histogram arrays preserve all retained events without embedding large event files in the report package.
hist=[]
for comp,markers in [('neutrophil',FUNCTIONAL),('stromal',STROMAL)]:
    for k in markers:
        values=[np.arcsinh(distributions[comp,m['sample'],k]/150) for m in bio]
        lo=min(v.min() for v in values);hi=max(v.max() for v in values);bins=np.linspace(lo-.01,hi+.01,81)
        for m,v in zip(bio,values):
            counts,_=np.histogram(v,bins)
            hist.append(dict(compartment=comp,marker=k,sample=m['sample'],group=m['group'],edges=bins.tolist(),counts=counts.tolist()))
        cv=np.arcsinh(control_events[comp,k]/150) if (comp,k) in control_events else None
        if cv is not None:
            counts,_=np.histogram(cv,bins);hist.append(dict(compartment=comp,marker=k,sample='FMO',group='FMO',edges=bins.tolist(),counts=counts.tolist(),outside_plot_n=int(((cv<bins[0])|(cv>bins[-1])).sum())))
(WORK/'histograms.json').write_text(json.dumps(hist))
print('PRIMARY FAMILY-ADJUSTED P < 0.1',flush=True)
print(T[(T.tier=='primary')&(T.p_holm_family<.1)][['family','endpoint','reference','comparison','mean_reference','mean_comparison','difference','p_holm_family','p_holm_all_primary_52']].to_string(index=False),flush=True)
print('COMPOSITION',P.groupby(['endpoint','group']).value.mean().unstack().round(3).to_string(),flush=True)
print('ASSOCIATIONS',A.sort_values('q_BH_21').head(10).round(4).to_string(index=False),flush=True)
print('VERIFY',jsonout['verification'],flush=True)
