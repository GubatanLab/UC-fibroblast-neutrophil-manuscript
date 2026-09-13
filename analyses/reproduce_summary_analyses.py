"""Recompute selected manuscript statistics from deposited analysis tables."""
from pathlib import Path
import itertools,json
import numpy as np
import pandas as pd
from scipy import stats

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/reproduction'
OUT.mkdir(parents=True,exist_ok=True)
checks=[]
def read(p):return pd.read_csv(ROOT/p)
def bh(p):return stats.false_discovery_control(np.asarray(p,float),method='bh')
def check(name,observed,expected):
    a=np.asarray(observed,float);b=np.asarray(expected,float)
    ok=bool(np.allclose(a,b,rtol=1e-7,atol=1e-10,equal_nan=True))
    checks.append({'check':name,'passed':ok,'comparisons':int(a.size),'max_abs_error':float(np.nanmax(np.abs(a-b)))})
    if not ok:raise AssertionError(f'{name}: observed and stored values differ')
def export(name,rows):pd.DataFrame(rows).to_csv(OUT/name,index=False)

froot='source_data/figure_01/program_panels_2026-09-12/'
fib=read(froot+'panel_D_paired_fibroblast_scores.csv')
ref=read(froot+'panel_D_all_8_program_statistics.csv')
rows=[]
for row in ref.itertuples():
    w=fib.pivot(index='PatientID',columns='condition',values=row.program).dropna()
    p=stats.wilcoxon(w['Inflamed UC'],w['Uninflamed UC'],method='approx',correction=True).pvalue
    rows.append({'program':row.program,'n_pairs':len(w),'p':p})
z=pd.DataFrame(rows);z['q']=bh(z.p)
check('Figure 1d Wilcoxon P',z.p,ref.p_value)
check('Figure 1d BH q',z.q,ref.FDR)
export('figure_1d.csv',z)

neut=read(froot+'panel_E_marker_excluded_patient_state_means.csv')
ref=read(froot+'panel_E_all_9_program_state_tests.csv')
states=['Neutrophil OSM','Neutrophil CXCR4','Neutrophil PADI4','Neutrophil MX1']
rows=[]
for row in ref.itertuples():
    w=neut.pivot(index='PatientID',columns='neut_state',values=row.program).dropna()
    test=stats.friedmanchisquare(*[w[s].to_numpy() for s in states])
    rows.append({'program':row.program,'n_patients':len(w),'statistic':test.statistic,'p':test.pvalue})
z=pd.DataFrame(rows);z['q']=bh(z.p)
check('Figure 1e Friedman statistic',z.statistic,ref.statistic)
check('Figure 1e Friedman P',z.p,ref.p_value_asymptotic)
check('Figure 1e BH q',z.q,ref.q_BH_all_9_programs)
export('figure_1e.csv',z)

flow=read('source_data/coculture_fluorescence/FACS_sample_expression.csv')
ref=read('source_data/figure_03/Canonical_FACS_RNA_2026-09-07/Figure3c_FACS_BH45.csv')
rows=[]
for row in ref.itertuples():
    a=flow.loc[(flow.marker==row.marker)&(flow.group==row.comparison),'median_fi'].to_numpy()
    b=flow.loc[(flow.marker==row.marker)&(flow.group==row.reference),'median_fi'].to_numpy()
    p=stats.permutation_test((a,b),lambda a,b:np.mean(a)-np.mean(b),n_resamples=np.inf,alternative='two-sided').pvalue
    rows.append({'marker':row.marker,'comparison':row.comparison,'reference':row.reference,'difference':a.mean()-b.mean(),'p':p})
z=pd.DataFrame(rows);z['q']=bh(z.p)
check('Figure 3c fluorescence differences',z.difference,ref.difference_fluorescence)
check('Figure 3c exact permutation P',z.p,ref.p_exact)
check('Figure 3c BH q across 45 tests',z.q,ref.q_bh_45)
export('figure_3c.csv',z)

eroot='source_data/figure_03/Canonical_FACS_RNA_2026-09-07/'
el=read(eroot+'Figure_3e_source_data.csv')
ref=read(eroot+'Figure_3e_statistics.csv')
# Historical source names use panel e; this endpoint is current main Figure 3d.
groups=list(dict.fromkeys(el.condition));rows=[]
for i,j in [(0,1),(1,2),(2,3)]:
    a=el.loc[el.condition==groups[i],'value'].to_numpy();b=el.loc[el.condition==groups[j],'value'].to_numpy()
    pooled=np.r_[a,b];difference=b.mean()-a.mean();null=[]
    for ix in itertools.combinations(range(len(pooled)),len(a)):
        mask=np.zeros(len(pooled),bool);mask[list(ix)]=True
        null.append(pooled[~mask].mean()-pooled[mask].mean())
    p=float(np.mean(np.abs(null)>=abs(difference)-1e-12))
    rows.append({'reference':groups[i],'comparison':groups[j],'difference':difference,'p':p})
z=pd.DataFrame(rows);order=np.argsort(z.p.to_numpy());adj=np.empty(len(z));adj[order]=np.minimum(1,np.maximum.accumulate(z.p.to_numpy()[order]*np.arange(len(z),0,-1)));z['p_holm']=adj
check('Figure 3d elastase difference',z.difference,ref.difference_mean)
check('Figure 3d exact permutation P',z.p,ref.p_exact)
check('Figure 3d Holm P',z.p_holm,ref.p_holm_3)
export('figure_3d.csv',z)

eroot='source_data/figure_06/epithelial_programs_2026-09-07/'
ep=read(eroot+'epithelial_program_scores.csv');ref=read(eroot+'epithelial_program_effects.csv');rows=[]
for row in ref.itertuples():
    z=ep[(ep.scope==row.scope)&(ep.state==row.state)&(ep.program==row.program)]
    a=z.loc[z.group=='DSS','score'].to_numpy();b=z.loc[z.group==row.group,'score'].to_numpy()
    rows.append({'scope':row.scope,'state':row.state,'program':row.program,'group':row.group,'effect':b.mean()-a.mean(),'p':stats.mannwhitneyu(b,a,alternative='two-sided',method='auto').pvalue})
z=pd.DataFrame(rows);z['q']=z.groupby('scope',group_keys=False).p.transform(lambda p:bh(p))
check('Figure 6 epithelial effects',z.effect,ref.effect)
check('Figure 6 epithelial Mann-Whitney P',z.p,ref.p)
check('Figure 6 epithelial scope-specific BH q',z.q,ref.q)
export('figure_6_epithelial.csv',z)

rroot='source_data/figure_07/receptor_availability_2026-09-07/'
d=read(rroot+'receptor_patient_differences.csv');ref=read(rroot+'receptor_effects.csv');rows=[]
for row in ref.itertuples():
    x=d[row.feature].dropna().to_numpy();n=len(x)
    signs=2*((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(float)-1
    p=float(np.mean(np.abs(signs@x/n)>=abs(x.mean())-1e-12))
    rows.append({'feature':row.feature,'n_patients':n,'effect':x.mean(),'p':p})
z=pd.DataFrame(rows);z['q']=bh(z.p)
check('Figure 7c receptor effects',z.effect,ref.effect)
check('Figure 7c exact sign-flip P',z.p,ref.p)
check('Figure 7c BH q',z.q,ref.q)
export('figure_7c.csv',z)
report={'status':'passed','checks':checks,'scope':'Deposited-table point estimates and P/q values; upstream raw-data processing and confidence-interval refits are not rerun.'}
(OUT/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({'status':'passed','checks':len(checks),'output':str(OUT)},indent=2))
