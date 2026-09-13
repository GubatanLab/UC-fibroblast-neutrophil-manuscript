"""Four-marker all-neutrophil ATN-161 expression analysis, with fixed contrasts."""
from pathlib import Path
import sys,json,hashlib,itertools,shutil,html
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np
import pandas as pd
from scipy import stats
OUT=ROOT/'outputs/05-expression-blockade';OUT.mkdir(exist_ok=True)
CACHE=ROOT/'analysis/regating/events'
CONFIG=ROOT/'analysis/expression_marker_assignment.json'
default={'confirmed':False,'assignment_source':'FCS labels; protocol conflicts; user confirmation pending','markers':{'OSM':'R670-A','CXCR4':'V610-A','MX1':'B710-A','PADI4':'B515-A'}}
config=json.loads(CONFIG.read_text()) if CONFIG.exists() else default
(OUT/'marker_assignments.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
MARKERS=config['markers'];ORDER=['N','N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK','N+FK','N+ATN']
PAIRS=[('N+IF','N+IF+ATN'),('N+IF+FK','N+IF+ATN+FK'),('N','N+ATN')]
inventory=pd.read_csv(ROOT/'analysis/results/raw_fcs_inventory.csv').set_index('file')
prior=pd.read_csv(ROOT/'outputs/02-gating-reanalysis/revised_marker_measurements.csv');prior=prior[(prior.experiment=='co')&(prior['quantile']==.99)]
rows=[]
for mp in sorted(CACHE.glob('co_*.json')):
    m=json.loads(mp.read_text())
    if not m['biological']:continue
    digest=hashlib.sha256((ROOT/m['file']).read_bytes()).hexdigest();assert digest==inventory.loc[m['file'],'sha256']
    z=np.load(CACHE/(m['key']+'.npz'));parent=z['old_3']
    for marker,ch in MARKERS.items():
        v=z['events'][parent,m['channels'].index(ch)];assert ch in m['matrix_info']['detectors'] and np.isfinite(v).all()
        old=prior[(prior['sample']==m['sample'])&(prior.channel==ch)].iloc[0]
        assert len(v)==old.neutrophil_events and np.isclose(np.median(v),old.median_fi,atol=1e-8,rtol=0)
        label=marker if config['confirmed'] or marker in ['OSM','CXCR4'] else f'{ch.split("-")[0]} ({marker} label; provisional)'
        rows.append(dict(marker=marker,label=label,channel=ch,sample=m['sample'],group=m['group'],neutrophil_events=len(v),median_fi=float(np.median(v)),mean_fi=float(v.mean()),q25=float(np.quantile(v,.25)),q75=float(np.quantile(v,.75)),negative_signal_percent=float(100*np.mean(v<0)),source_file=m['file'],source_sha256=digest,compensation_sha256=m['matrix_info']['matrix_sha256']))
R=pd.DataFrame(rows);assert len(R)==188 and R.groupby('marker').size().eq(47).all()
R.to_csv(OUT/'Expression_sample_measurements.csv',index=False)
S=R.groupby(['marker','label','channel','group'],sort=False).agg(n=('sample','size'),mean_sample_median=('median_fi','mean'),sd_sample_median=('median_fi','std'),median_sample_median=('median_fi','median'),mean_sample_mean=('mean_fi','mean'),min_neutrophil_events=('neutrophil_events','min'),max_neutrophil_events=('neutrophil_events','max')).reset_index()
S.to_csv(OUT/'Expression_group_summary.csv',index=False)
tests=[]
for metric in ['median_fi','mean_fi']:
    for marker,ch in MARKERS.items():
        for ref,trt in PAIRS:
            x=R[(R.marker==marker)&(R.group==trt)][metric].to_numpy();y=R[(R.marker==marker)&(R.group==ref)][metric].to_numpy();pooled=np.r_[x,y];diff=x.mean()-y.mean()
            allocations=np.array(list(itertools.combinations(range(len(pooled)),len(x))));sums=pooled[allocations].sum(axis=1);perm=sums/len(x)-(pooled.sum()-sums)/len(y)
            p=float(np.mean(abs(perm)>=abs(diff)-max(abs(diff)*1e-12,1e-9)))
            check=stats.permutation_test((x,y),lambda a,b:np.mean(a)-np.mean(b),vectorized=False,n_resamples=np.inf,alternative='two-sided').pvalue;assert np.isclose(p,check,atol=1e-12,rtol=0)
            vx=x.var(ddof=1)/len(x);vy=y.var(ddof=1)/len(y);df=(vx+vy)**2/(vx**2/(len(x)-1)+vy**2/(len(y)-1));half=stats.t.ppf(.975,df)*np.sqrt(vx+vy)
            tests.append(dict(marker=marker,label=R[R.marker==marker].label.iloc[0],channel=ch,metric=metric,reference=ref,treated=trt,n_reference=len(y),n_treated=len(x),mean_reference=y.mean(),mean_treated=x.mean(),difference=diff,ci95_low=diff-half,ci95_high=diff+half,relative_change_percent=100*diff/y.mean() if min(x.mean(),y.mean())>0 else np.nan,p_exact=p,allocations=len(perm)))
T=pd.DataFrame(tests)
def holm(ids,column):
    ix=T.loc[ids].sort_values('p_exact',kind='stable').index;running=0
    for j,i in enumerate(ix):running=max(running,min(1,T.at[i,'p_exact']*(len(ix)-j)));T.at[i,column]=running
for metric in ['median_fi','mean_fi']:
    holm(T[T.metric==metric].index,'p_holm_12')
    for marker in MARKERS:holm(T[(T.metric==metric)&(T.marker==marker)].index,'p_holm_within_marker_3')
T.to_csv(OUT/'Expression_blockade_comparisons.csv',index=False)
primary=T[T.metric=='median_fi']
def table(headers,data):return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in data)
status=('The user confirmed the detector assignments recorded in marker_assignments.json.' if config['confirmed'] else 'MX1/PADI4 identities remain provisional. FCS annotations assign B710-A to MX1 and B515-A to PADI4, whereas the protocol reverses them. The numerical findings belong to the physical detector until this conflict is resolved.')
findings=[]
for r in primary[primary.reference=='N+IF'].itertuples():
    findings.append(f'- **{r.label}:** {r.mean_reference:,.1f} to {r.mean_treated:,.1f} a.u. with ATN-161; Holm-adjusted P = {r.p_holm_12:.4f}.')
report='''# OSM, CXCR4, MX1 and PADI4: ATN-161 expression analysis

All 47 co-culture biological samples were analyzed on the original scatter/singlet/CD16/CD11b neutrophil parent, including CXCR4-negative cells. Each sample contributes one median compensated fluorescence value per marker. No marker-positive gate or CXCR4 restriction was applied. These are antibody-expression readouts, not enzyme-activity measurements.

'''+status+'''

## Main findings in inflamed co-culture

'''+ '\n'.join(findings)+'''

The response is not uniform suppression of all measured markers: OSM and CXCR4 medians fall, while B515 and B710 medians rise. All four N+IF versus N+IF+ATN contrasts are also supported by the secondary arithmetic-mean fluorescence analysis after its own 12-test correction. In the FK-866 background, median OSM and CXCR4 decline and B515 rises; B710 does not reach adjusted P < 0.05. Of these FK-background contrasts, only B515 remains significant using arithmetic-mean fluorescence. In neutrophils alone, median CXCR4 rises with ATN-161, whereas the mean-fluorescence difference is negative and nonsignificant; this endpoint is distribution-dependent and does not support a uniform CXCR4 suppression claim.

Negative treated OSM medians represent compensated signal below zero, not negative expression or a reduction exceeding 100%. Detector identities must be resolved before attaching the B515/B710 findings definitively to MX1 or PADI4.

## Primary comparisons

Two-sided exact unpaired label-permutation tests compare group means of sample medians, using 924 allocations per comparison. Holm adjustment covers all 12 primary comparisons (four markers by three ATN-161 contrasts). Arithmetic mean fluorescence is secondary and uses a separate 12-comparison adjustment. Within-marker three-comparison adjustments are included only as sensitivity results for comparability with the earlier MPO analysis; headline results and figure annotations use the 12-comparison adjustment. This focused follow-up family does not adjust jointly across all historical manuscript analyses, and was not prospectively registered. Pairing between conditions remains unconfirmed. All direct contrasts have six independent samples per condition.

'''+table(['Marker / signal','Reference','ATN-containing','Mean reference','Mean treated','Difference (95% CI), a.u.','Holm P (12)'],[[r.label,r.reference,r.treated,f'{r.mean_reference:.2f}',f'{r.mean_treated:.2f}',f'{r.difference:.2f} ({r.ci95_low:.2f}, {r.ci95_high:.2f})',f'{r.p_holm_12:.5f}'] for r in primary.itertuples()])+'''

Welch 95% confidence intervals are unadjusted and descriptive. Relative changes in the CSV describe group means of sample medians only and are omitted when either mean is nonpositive; they are not calibrated or background-subtracted protein fold changes. Negative compensated fluorescence values are preserved, and do not represent negative protein content.

## Secondary arithmetic-mean intensity check

'''+table(['Marker / signal','Reference','ATN-containing','Mean reference','Mean treated','Holm P (12)'],[[r.label,r.reference,r.treated,f'{r.mean_reference:.2f}',f'{r.mean_treated:.2f}',f'{r.p_holm_12:.5f}'] for r in T[T.metric=='mean_fi'].itertuples()])+'''

## All conditions

'''+table(['Marker / signal','Condition','n','Mean sample median','SD'],[[r.label,r.group,r.n,f'{r.mean_sample_median:.2f}',f'{r.sd_sample_median:.2f}'] for r in S.itertuples()])+'''

N = neutrophils alone; NF/IF = non-inflamed/inflamed fibroblasts; FK = FK-866; ATN = ATN-161. N+NF has five samples because P5-5 is absent. There is no N+NF+ATN arm, so blockade in non-inflamed co-culture cannot be tested. No event or biological sample was excluded on the basis of these outcomes.

## Interpretation limits and verification

The original parent gates and compensation matrices were retained, with the corrected detector-alias import. A consistent live-neutrophil population cannot be reconstructed because 44 of 47 files lack the viability channel. Matched marker FMOs are unavailable; these analyses quantify fluorescence without assigning positive/negative status. Treatment and plate/well blocks can confound comparisons. Neutrophil-only contrasts provide context but do not establish cell-specific mediation, inhibitor synergy or a target-specific mechanism.

All 47 source-file hashes matched the previous raw-file inventory. All 188 marker medians and parent counts matched the preceding corrected reanalysis. All 24 exact p-values matched an independent SciPy permutation calculation. The original source files and earlier analyses are unchanged.

Source data and complete primary/secondary comparisons are included as CSV files. Figures use an explicit median-based representative selection rule, show all biological replicates, and label provisional detector assignments. Rerun analysis/analyze_expression_blockade.py in the original project with the existing corrected event cache and documented scientific dependencies. Marker-assignment configuration, if present, is read from analysis/expression_marker_assignment.json.
'''
(OUT/'Expression_blockade_report.md').write_text(report,encoding='utf-8')
(OUT/'verification.json').write_text(json.dumps({'source_hashes_checked':47,'marker_medians_and_parent_counts_checked':188,'independent_exact_p_values_checked':24,'all_passed':True},indent=2),encoding='utf-8')
shutil.copy2(__file__,OUT/Path(__file__).name)
print(primary[['label','reference','treated','mean_reference','mean_treated','difference','p_holm_12']].round(5).to_string(index=False))
print('SECONDARY',T[T.metric=='mean_fi'][['label','reference','difference','p_holm_12']].round(5).to_string(index=False))
