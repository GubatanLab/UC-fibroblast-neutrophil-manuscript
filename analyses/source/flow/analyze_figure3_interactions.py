from pathlib import Path
import sys,json,shutil
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd
from scipy import stats
OUT=ROOT/'outputs/07-figure3-preview';OUT.mkdir(exist_ok=True)
R=pd.read_csv(ROOT/'outputs/05-expression-blockade/Expression_sample_measurements.csv')
M=pd.read_csv(ROOT/'outputs/03-mpo-blockade/MPO_sample_measurements.csv').rename(columns={'median_mpo':'median_fi','mean_mpo':'mean_fi'});M['marker']='MPO';M['channel']='YG610-A'
R=pd.concat([R,M],ignore_index=True);ORDER=['OSM','CXCR4','MX1','PADI4','MPO'];GROUPS=['N','N+ATN','N+IF','N+IF+ATN'];coeff=np.array([1,-1,-1,1])
rows=[]
for metric in ['median_fi','mean_fi']:
    for marker in ORDER:
        sub=R[R.marker==marker];arrays=[sub[sub.group==g][metric].to_numpy() for g in GROUPS];assert all(len(x)==6 for x in arrays)
        means=np.array([x.mean() for x in arrays]);v=np.array([x.var(ddof=1)/len(x) for x in arrays]);estimate=float(coeff@means);se=np.sqrt(v.sum());df=v.sum()**2/sum(vv**2/(len(x)-1) for vv,x in zip(v,arrays));p=2*stats.t.sf(abs(estimate/se),df);half=stats.t.ppf(.975,df)*se
        scale=np.sqrt((arrays[0].var(ddof=1)+arrays[2].var(ddof=1))/2)
        # Independent 2x2 OLS/HC2 check of estimate and SE. Satterthwaite df remains primary.
        X=np.repeat(np.array([[1,0,0,0],[1,0,1,0],[1,1,0,0],[1,1,1,1]]),6,axis=0);y=np.concatenate(arrays);inv=np.linalg.inv(X.T@X);beta=inv@X.T@y;res=y-X@beta;h=np.einsum('ij,jk,ik->i',X,inv,X);cov=inv@X.T@np.diag(res**2/(1-h))@X@inv
        assert np.isclose(beta[3],estimate) and np.isclose(np.sqrt(cov[3,3]),se)
        rows.append(dict(marker=marker,channel=sub.channel.iloc[0],metric=metric,n_per_group=6,effect_in_neutrophils=means[1]-means[0],effect_in_coculture=means[3]-means[2],interaction=estimate,se=se,df=df,ci95_low=estimate-half,ci95_high=estimate+half,p=p,reference_sd=scale,standardized_interaction=estimate/scale,standardized_ci_low=(estimate-half)/scale,standardized_ci_high=(estimate+half)/scale))
T=pd.DataFrame(rows)
def holm(frame,pcol,outcol):
    ids=frame.sort_values(pcol,kind='stable').index;run=0
    for rank,i in enumerate(ids):run=max(run,min(1,frame.at[i,pcol]*(len(ids)-rank)));frame.at[i,outcol]=run
    return frame
T=pd.concat([holm(T[T.metric==metric].copy(),'p','p_holm_5') for metric in ['median_fi','mean_fi']]);T.to_csv(OUT/'FACS_context_interactions.csv',index=False)
D=pd.read_csv(ROOT/'outputs/05-expression-blockade/Expression_blockade_comparisons.csv');D=D[(D.metric=='median_fi')&(D.reference=='N+IF')].copy()
dm=pd.read_csv(ROOT/'outputs/03-mpo-blockade/MPO_blockade_comparisons.csv');dm=dm[(dm.metric=='median_mpo')&(dm.reference=='N+IF')].copy();dm['marker']='MPO';dm['channel']='YG610-A';dm['metric']='median_fi'
D=pd.concat([D,dm],ignore_index=True);D=holm(D,'p_exact','p_holm_5_figure3')
for i,r in D.iterrows():
    sd=T[(T.metric=='median_fi')&(T.marker==r.marker)].iloc[0].reference_sd
    D.at[i,'reference_sd']=sd
    for col in ['difference','ci95_low','ci95_high']:D.at[i,'standardized_'+col]=r[col]/sd
D.to_csv(OUT/'FACS_coculture_effects_figure3.csv',index=False)
R.to_csv(OUT/'FACS_sample_expression.csv',index=False)
def table(frame,cols):return '| '+' | '.join(cols)+' |\n|'+'|'.join(['---']*len(cols))+'|\n'+'\n'.join('| '+' | '.join(str(round(x,5)) if isinstance(x,(float,np.floating)) else str(x) for x in row)+' |' for row in frame[cols].itertuples(index=False,name=None))
report='''# Figure 3 FACS context-interaction analysis

## Design and estimand

The primary estimand is [N+IF+ATN - N+IF] - [N+ATN - N] for each of five markers. Each of the four groups contains six independent biological samples. Donor matching across FACS conditions has not been established, so this analysis is unpaired and is not assigned the RNA experiment's donor matching. No RNA/FACS cross-donor correlation is calculated.

Primary outcome: each sample's median compensated fluorescence in all gated neutrophils, without a CXCR4-positive restriction. Primary interaction inference uses a heteroskedastic four-group contrast: standard error sqrt(sum(s_g^2/n_g)), Welch-Satterthwaite degrees of freedom, a two-sided t test and unadjusted 95% confidence interval. Holm correction covers five primary marker interactions. Mean fluorescence is secondary with a separate five-test correction. The coefficient and standard error were independently verified by a saturated 2x2 linear model with HC2 robust covariance for all ten tests. This is an exploratory analysis; its assumptions and small sample size limit inference. Naive permutation of all treatment labels was not used because main effects are nuisance effects under an interaction-only null.

Positive interaction means the ATN effect is more positive (or less negative) in co-culture. Negative interaction means the ATN effect is more negative (or less positive) in co-culture. An interaction demonstrates statistical context dependence under the design assumptions, not fibroblast-specific molecular mediation. Treatment/plate confounding cannot be separated from these contrasts with this design.

## Primary median-expression interactions

'''+table(T[T.metric=='median_fi'],['marker','effect_in_neutrophils','effect_in_coculture','interaction','ci95_low','ci95_high','p_holm_5'])+'''

## Secondary mean-expression interactions

'''+table(T[T.metric=='mean_fi'],['marker','interaction','ci95_low','ci95_high','p_holm_5'])+'''

## Figure 3 within-co-culture effects

Five pre-specified marker comparisons (ATN versus no ATN in inflamed co-culture) are displayed with a separate five-marker Holm correction of the previously calculated exact two-sided permutation p-values. Their sample values and estimates are unchanged. This Figure 3 family is narrower than the previous 12-comparison expression family and different from the previous three-comparison MPO family; it is explicitly exploratory and does not replace the archived analyses. No correction was chosen by its resulting significance. All original comparisons remain available in the earlier packages.

'''+table(D,['marker','difference','ci95_low','ci95_high','p_exact','p_holm_5_figure3'])+'''

## Figure scales and retained limitations

Forest plots divide effects and confidence bounds by each marker's pooled within-group SD across the untreated N and N+IF groups. This fixed descriptive scale is used to display different detectors together; intervals are scaled raw-effect intervals and do not incorporate uncertainty in the reference SD. They are not Hedges g confidence intervals, fold changes or calibrated protein quantities. Raw-unit results are provided above.

MX1=B710 and PADI4=B515 follow FCS labels at the user's request; the protocol reverses these assignments and the source conflict remains unresolved. MPO was confirmed to be an antibody stain. The original gates and compensation are retained; viability is absent in 44/47 files. Antibody medians do not measure signaling flux, enzymatic activity, NET formation or survival. The retained elastase panel contains no blockade arm, its source statistical metadata require confirmation, and it is not a newly analyzed experiment. No contact/transwell, phospho-flow, neutralization, migration or new functional experiment was performed.
'''
(OUT/'Figure3_analysis_report.md').write_text(report,encoding='utf-8');shutil.copy2(__file__,OUT/Path(__file__).name)
print(T[T.metric=='median_fi'][['marker','interaction','ci95_low','ci95_high','p_holm_5']].round(5).to_string(index=False));print('MEAN',T[T.metric=='mean_fi'][['marker','interaction','p_holm_5']].round(5).to_string(index=False))
