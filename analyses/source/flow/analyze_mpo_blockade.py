"""All-neutrophil MPO antibody intensity and focused ATN-161 contrasts."""
from pathlib import Path
import sys, json, hashlib, itertools, html, shutil
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=ROOT/'outputs/03-mpo-blockade';OUT.mkdir(exist_ok=True)
CACHE=ROOT/'analysis/regating/events'
ORDER=['N','N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK','N+FK','N+ATN']
CONTRASTS=[('N+IF+ATN','N+IF'),('N+IF+ATN+FK','N+IF+FK'),('N+ATN','N')]
inventory=pd.read_csv(ROOT/'analysis/results/raw_fcs_inventory.csv').set_index('file')
qc=pd.read_csv(ROOT/'outputs/02-gating-reanalysis/sample_qc.csv')
rows=[];distributions={g:[] for g in ORDER}
for mp in sorted(CACHE.glob('co_*.json')):
    m=json.loads(mp.read_text())
    if not m['biological']:continue
    source=ROOT/m['file'];digest=hashlib.sha256(source.read_bytes()).hexdigest()
    assert digest==inventory.loc[m['file'],'sha256'],f'Source changed: {source}'
    assert 'YG610-A' in m['matrix_info']['detectors']
    assert m['markers'][m['channels'].index('YG610-A')].lower()=='mpo'
    z=np.load(CACHE/(m['key']+'.npz'));mask=z['old_3'];v=z['events'][mask,m['channels'].index('YG610-A')]
    assert len(v)==int(qc[(qc.experiment=='co')&(qc['sample']==m['sample'])].iloc[0].neutrophil_events)
    assert len(v)>0 and np.isfinite(v).all()
    rows.append(dict(sample=m['sample'],group=m['group'],neutrophil_events=len(v),median_mpo=float(np.median(v)),mean_mpo=float(np.mean(v)),q25_mpo=float(np.quantile(v,.25)),q75_mpo=float(np.quantile(v,.75)),negative_signal_percent=float(100*np.mean(v<0)),source_file=m['file'],source_sha256=digest,compensation_sha256=m['matrix_info']['matrix_sha256']))
    distributions[m['group']].append(v)
R=pd.DataFrame(rows);assert len(R)==47 and R['sample'].is_unique
R.to_csv(OUT/'MPO_sample_measurements.csv',index=False)
summary=[]
for g in ORDER:
    x=R[R.group==g]
    summary.append(dict(group=g,n=len(x),mean_sample_median=x.median_mpo.mean(),sd_sample_median=x.median_mpo.std(),median_sample_median=x.median_mpo.median(),min_parent_events=x.neutrophil_events.min(),max_parent_events=x.neutrophil_events.max(),mean_sample_mean=x.mean_mpo.mean()))
S=pd.DataFrame(summary);S.to_csv(OUT/'MPO_group_summary.csv',index=False)
tests=[]
for metric in ['median_mpo','mean_mpo']:
    family=[]
    for a,b in CONTRASTS:
        x=R[R.group==a][metric].to_numpy();y=R[R.group==b][metric].to_numpy();diff=x.mean()-y.mean();pooled=np.r_[x,y]
        perms=[]
        for ids in itertools.combinations(range(len(pooled)),len(x)):
            total=pooled[list(ids)].sum();perms.append(total/len(x)-(pooled.sum()-total)/len(y))
        p=float(np.mean(np.abs(perms)>=abs(diff)-1e-9))
        vx=x.var(ddof=1)/len(x);vy=y.var(ddof=1)/len(y);se=np.sqrt(vx+vy);df=(vx+vy)**2/(vx*vx/(len(x)-1)+vy*vy/(len(y)-1));span=stats.t.ppf(.975,df)*se
        family.append(dict(metric=metric,treated=a,reference=b,n_treated=len(x),n_reference=len(y),mean_treated=x.mean(),mean_reference=y.mean(),difference=diff,ci95_low=diff-span,ci95_high=diff+span,relative_change_percent=100*diff/y.mean() if y.mean()>0 and x.mean()>0 else np.nan,p_exact=p,allocations=len(perms),p_welch=float(stats.ttest_ind(x,y,equal_var=False).pvalue)))
    ix=np.argsort([r['p_exact'] for r in family]);running=0
    for rank,i in enumerate(ix):
        running=max(running,min(1,family[i]['p_exact']*(len(ix)-rank)));family[i]['p_holm']=running
    tests.extend(family)
T=pd.DataFrame(tests);T.to_csv(OUT/'MPO_blockade_comparisons.csv',index=False)
verification=[]
for r in T.itertuples():
    x=R[R.group==r.treated][r.metric].to_numpy();y=R[R.group==r.reference][r.metric].to_numpy()
    independent=stats.permutation_test((x,y),lambda a,b:np.mean(a)-np.mean(b),vectorized=False,n_resamples=np.inf,alternative='two-sided').pvalue
    assert abs(independent-r.p_exact)<1e-12
    verification.append(dict(check='Independent scipy exact permutation',metric=r.metric,comparison=r.treated+' vs '+r.reference,passed=True))
# Fresh raw-file reapplication checks the new MPO endpoint in two representative conditions.
import flowkit as fk
w=fk.Workspace(str(ROOT/'analysis/regating/co_canonical.wsp'),load_missing_file_data=True)
for sid in ['P3-1','P4-1']:
    row=R[R['sample']==sid].iloc[0];sample=fk.Sample(str(ROOT/row.source_file))
    result=w.get_gating_strategy(sample.id).gate_sample(sample,cache_events=False)
    membership=result.get_gate_membership('neutrophil',gate_path='root/scatter/s1/s2')
    sample.apply_compensation(w.get_comp_matrix(sample.id))
    value=float(np.median(sample.get_events(source='comp')[membership,sample.pnn_labels.index('YG610-A')]))
    assert np.isclose(value,row.median_mpo,rtol=0,atol=1e-8)
    assert int(membership.sum())==row.neutrophil_events
    verification.append(dict(check='Fresh raw-file MPO median and parent count',sample=sid,passed=True,median_mpo=value))
(OUT/'verification.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')
plt.rcParams.update({'font.family':'Arial','font.size':11,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
colors=['#596673','#3E8196','#BE5E35','#D39252','#327E69','#735890','#959FAA','#99BFA9']
fig,ax=plt.subplots(figsize=(11,6))
for i,g in enumerate(ORDER):
    x=R[R.group==g].sort_values('sample').median_mpo.to_numpy()
    ax.scatter(i+np.linspace(-.12,.12,len(x)),x,s=48,c=colors[i],edgecolors='white',zorder=3)
    ax.errorbar(i,x.mean(),yerr=x.std(ddof=1),fmt='_',color='black',capsize=4,markersize=19)
ax.set_xticks(range(8),['N','N + NF','N + IF','N + IF\n+ FK','N + IF\n+ ATN','N + IF\n+ ATN + FK','N + FK','N + ATN'])
ax.set_ylabel('MPO median fluorescence per sample (a.u.)');ax.grid(axis='y',alpha=.2)
fig.suptitle('MPO antibody staining in all gated neutrophils',x=.1,ha='left',fontweight='bold',fontsize=17)
fig.text(.1,.03,'Points: independent biological samples; bars: mean ± SD. Original CD16/CD11b neutrophil parent; no CXCR4 restriction.\nATN = ATN-161; FK = FK-866; NF/IF = non-inflamed/inflamed fibroblasts. Antibody signal does not measure enzyme activity.',fontsize=9)
fig.subplots_adjust(bottom=.22,top=.88,left=.1,right=.97)
for ext in ['png','svg']:fig.savefig(OUT/f'MPO_all_neutrophils.{ext}',dpi=180)
plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(13,4.8))
for ax,(a,b) in zip(axes,CONTRASTS):
    transformed=[np.arcsinh(v/150) for group in [a,b] for v in distributions[group]]
    edges=np.linspace(min(v.min() for v in transformed)-.01,max(v.max() for v in transformed)+.01,151)
    for group,color in [(b,'#BE5E35'),(a,'#327E69')]:
        vals=distributions[group]
        hist=np.mean([np.histogram(np.arcsinh(v/150),bins=edges,density=True)[0] for v in vals],axis=0)
        ax.plot((edges[1:]+edges[:-1])/2,hist,color=color,label=group)
    ax.set_title(a+' vs '+b,fontsize=11);ax.set_xlabel('asinh(compensated MPO / 150)');ax.legend(fontsize=9)
axes[0].set_ylabel('Mean within-sample density')
fig.suptitle('MPO distributions: equal weight per biological sample',fontweight='bold')
fig.text(.08,.025,'All neutrophils; curves summarize the antibody-signal distribution and are not enzyme-activity measurements.',fontsize=9)
fig.subplots_adjust(bottom=.18,top=.84,wspace=.3)
for ext in ['png','svg']:fig.savefig(OUT/f'MPO_distributions.{ext}',dpi=180)
plt.close(fig)
def table(headers,rs):return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in rs)
primary=T[T.metric=='median_mpo']
report='''# Effect of ATN-161 on MPO antibody staining in co-culture neutrophils

The user confirmed that MPO was measured using an antibody. These results quantify MPO-associated fluorescence, interpreted as an intracellular abundance readout; they do not measure catalytic activity, extracellular MPO release or NET formation. The exact antibody clone, conjugate and staining conditions remain unspecified in the supplied protocol.

ATN-161 was associated with lower median MPO antibody fluorescence in inflamed fibroblast co-culture: 3,162 to 562 a.u. (82.2% lower; Holm-adjusted exact p=0.0195). In the FK-866 background, adding ATN-161 reduced the signal from 1,611 to 499 a.u. (69.0% lower; adjusted p=0.0303). Both co-culture contrasts were also supported by the secondary arithmetic-mean intensity analysis (adjusted p=0.00649 each). This supports an ATN-associated reduction in MPO staining across the whole gated neutrophil population under these conditions.

Neutrophils alone also showed lower sample-median MPO signal with ATN-161 (629 to -68 a.u.; adjusted p=0.0195), but the secondary arithmetic-mean endpoint was not significant after adjustment (p=0.0931). The negative compensated median represents signal near/below the compensated zero and is not negative protein content; a percentage reduction is omitted for this comparison. The neutrophil-only observation means the co-culture effect cannot be attributed exclusively to fibroblast mediation.

All 47 available biological samples were analyzed on YG610-A, using compensated events from the corrected workspace import. The parent is scatter/s1/s2/neutrophil (CD16/CD11b), including both CXCR4-positive and CXCR4-negative cells. There is no MPO-positive threshold or subgroup selection. Original sample identities and group assignments are preserved; N+NF has five samples because P5-5 is absent. Each of the other seven groups has six independent samples. Forty-four of 47 files lack the viability channel, so this is all gated neutrophils rather than a consistently verified live-neutrophil population.

## Primary group summaries

The primary endpoint is the median compensated MPO fluorescence within each biological sample. Group means and standard deviations summarize those sample medians; individual events are not independent replicates.

'''+table(['Condition','n','Mean sample median (a.u.)','SD (a.u.)','Neutrophil events per sample'],[[r.group,int(r.n),f'{r.mean_sample_median:.2f}',f'{r.sd_sample_median:.2f}',f'{int(r.min_parent_events)}–{int(r.max_parent_events)}'] for r in S.itertuples()])+'''

## Direct blockade comparisons

Differences are ATN-containing minus the corresponding ATN-free group. Two-sided exact label-permutation tests of sample-mean differences are Holm-adjusted across these three focused comparisons. The 95% Welch confidence intervals are unadjusted and descriptive. Pairing across conditions is unconfirmed, so tests are unpaired. These are exploratory follow-up tests selected for the current MPO question, not a prospectively registered family or a correction across every earlier manuscript analysis.

'''+table(['Comparison','Difference (a.u.)','95% CI','Relative change','Exact p','Holm p'],[[r.treated+' vs '+r.reference,f'{r.difference:.2f}',f'{r.ci95_low:.2f} to {r.ci95_high:.2f}',f'{r.relative_change_percent:.1f}%' if np.isfinite(r.relative_change_percent) else 'Not reported: nonpositive treated mean',f'{r.p_exact:.6f}',f'{r.p_holm:.6f}'] for r in primary.itertuples()])+'''

Relative changes describe the ratio of group mean sample medians and are not background-subtracted protein fold changes. Compensation can yield negative values, which were retained. A secondary analysis of arithmetic mean fluorescence uses its own three-comparison Holm family and is supplied in MPO_blockade_comparisons.csv; it is sensitive to bright tails and does not replace the primary median endpoint.

## Scope and interpretation

N+IF+ATN versus N+IF estimates the ATN-associated change in inflamed fibroblast co-culture. N+IF+ATN+FK versus N+IF+FK estimates the additional ATN-associated change in the FK-866 background. N+ATN versus N evaluates the neutrophil-only control. N+NF is shown descriptively; no N+NF+ATN condition is available. These comparisons do not demonstrate inhibitor synergy or uniquely identify which cell type mediates an effect. The study design and reagent specificity do not establish target-specific causality from this analysis alone.

The original scatter, singlet and lineage gates and compensation matrix were retained. Group assignment coincides with well/plate blocks, which can contribute technical confounding. No outcome-based sample exclusions were applied. Missing uniform viability, unverified antibody staining comparability and compensation uncertainty limit biological interpretation. No matching MPO FMO was found, so an MPO-positive percentage was not used as the primary endpoint.

## Reproducibility and files

Every analyzed FCS file was verified against the previous SHA-256 source inventory; parent-event counts matched the preceding reanalysis. The source MPO annotation and inclusion of YG610-A in the compensation matrix were checked for all 47 samples. Fresh raw-file reapplication reproduced the MPO median and parent count for P3-1 and P4-1, and an independent SciPy exact-permutation calculation matched all six primary/secondary p-values. This analysis reuses the corrected compensated-event cache under analysis/regating/events. Original FCS files and prior outputs were not changed.

MPO_sample_measurements.csv contains per-sample medians, means, quartiles, event counts and source hashes. MPO_group_summary.csv contains all eight conditions. MPO_blockade_comparisons.csv contains primary and secondary comparisons. MPO_all_neutrophils.png/svg and MPO_distributions.png/svg provide editable scientific figures. The analysis source is included for reproduction in the same project structure with the previously documented Python dependencies.
'''
(OUT/'MPO_blockade_report.md').write_text(report,encoding='utf-8')
body='<h1>ATN-161 and MPO antibody staining</h1><p>All gated neutrophils; 47 biological samples. MPO abundance readout, not enzymatic activity.</p>'
body+='<img src="MPO_all_neutrophils.png" alt="MPO antibody intensity across all conditions">'
body+=S.round(3).to_html(index=False)+primary.round(6).to_html(index=False)
body+='<img src="MPO_distributions.png" alt="MPO distributions for the three blockade comparisons">'
body+='<h2>Full results and methods</h2><pre>'+html.escape(report)+'</pre>'
(OUT/'MPO_blockade_report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><title>MPO blockade analysis</title><style>body{font:16px/1.55 Arial;max-width:1200px;margin:35px auto;padding:20px;color:#233746}img{max-width:100%}table{font-size:12px;border-collapse:collapse;display:block;overflow:auto}td,th{padding:7px;border:1px solid #ccd8df}pre{white-space:pre-wrap;font:14px/1.6 Arial}</style>'+body+'</html>',encoding='utf-8')
shutil.copy2(__file__,OUT/'analyze_mpo_blockade.py')
print(S.round(3).to_string(index=False));print(primary.round(6).to_string(index=False));print('Secondary mean intensity:',T[T.metric=='mean_mpo'][['treated','reference','difference','p_holm']].round(6).to_string(index=False))
