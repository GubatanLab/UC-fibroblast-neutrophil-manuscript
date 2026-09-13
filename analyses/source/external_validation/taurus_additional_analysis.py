from pathlib import Path
import sys,json
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'tmp/canonical_nature_20260903/packages'))
import numpy as np,pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=R/'output/Figure7_Additional_TAURUS_2026-09-07';rng=np.random.default_rng(7092026)
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none'})
def bh(p):
 p=np.asarray(p);i=np.argsort(p);q=np.empty(len(p));q[i]=np.minimum(1,np.minimum.accumulate((p[i]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def paired(d):
 d=np.asarray(d,float);n=len(d);assert n>=3
 if n<=20:
  signs=2*((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(float)-1
  p=float(np.mean(np.abs(signs@d/n)>=abs(d.mean())-1e-12))
 else:
  null=(rng.choice([-1,1],(99999,n))*d).mean(1);p=(1+np.sum(abs(null)>=abs(d.mean())-1e-12))/100000
 ci=np.quantile(d[rng.integers(n,size=(10000,n))].mean(1),[.025,.975])
 return dict(n_patients=n,effect=d.mean(),lo=ci[0],hi=ci[1],p=p)
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(O/(name+'.'+ext),dpi=320)
 plt.close(fig)
def forest(ax,df,labels,title,xlabel,color='#8870AD'):
 for i,row in enumerate(df.itertuples()):
  ax.plot([row.lo,row.hi],[i,i],color=color,lw=1.3);ax.scatter(row.effect,i,s=24,color=color,zorder=3)
  ax.text(1.02,i,f'{row.q:.3f}' if row.q>=.001 else '<0.001',transform=ax.get_yaxis_transform(),va='center',fontsize=7)
 ax.axvline(0,color='#AAB0B7',lw=.7);ax.set_yticks(range(len(df)),labels);ax.invert_yaxis();ax.set_ylim(len(df)-.5,-.5);ax.set_title(title,loc='left',fontsize=9,pad=18);ax.set_xlabel(xlabel);ax.text(1.02,1.025,'q',transform=ax.transAxes);ax.grid(axis='x',color='#E8EBEF',lw=.5);ax.set_axisbelow(True);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
fib=pd.read_csv(O/'fib_UC_targeted_cells.csv.gz',low_memory=False)
fib=fib.loc[fib.minor.ne('Pericyte')].copy();fib['activated']=fib.final_analysis.eq('THY1pos FAPpos PDPNpos fibroblast')
genes=['ITGA5','ITGB1','OSMR','IL6ST'];features=[g+'_log' for g in genes]+['ITGA5_ITGB1_pct','OSMR_IL6ST_pct']
fib['ITGA5_ITGB1_pct']=100*(fib.ITGA5.gt(0)&fib.ITGB1.gt(0));fib['OSMR_IL6ST_pct']=100*(fib.OSMR.gt(0)&fib.IL6ST.gt(0))
keys=['sample_id','Patient','Site','Treatment','Inflammation'];means=fib.groupby(keys+['activated'])[features].mean();counts=fib.groupby(keys+['activated']).size().rename('n_cells');agg=means.join(counts).reset_index();agg.to_csv(O/'receptor_biopsy_scores.csv',index=False)
a=agg.loc[agg.activated & agg.n_cells.ge(5)];b=agg.loc[~agg.activated & agg.n_cells.ge(20)];m=a.merge(b,on=keys,suffixes=('_act','_other'));assert len(m)
diff=pd.DataFrame({'Patient':m.Patient,'sample_id':m.sample_id})
for g in features:diff[g]=m[g+'_act']-m[g+'_other']
patient=diff.groupby('Patient')[features].mean();patient.to_csv(O/'receptor_patient_differences.csv');results=pd.DataFrame([{'feature':g,**paired(patient[g])} for g in features]);results['q']=bh(results.p);results.to_csv(O/'receptor_effects.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(183/25.4,76/25.4),gridspec_kw={'width_ratios':[1.1,1]});fig.subplots_adjust(left=.12,right=.92,bottom=.24,top=.69,wspace=.98)
forest(axs[0],results.iloc[:4],genes,'a  Receptor-subunit expression','Activated FAP − other fibroblasts\n(mean log-normalized RNA)')
forest(axs[1],results.iloc[4:],['ITGA5 + ITGB1','OSMR + IL6ST'],'b  Same-cell RNA co-detection','Activated FAP − other fibroblasts\n(percentage points)')
fig.text(.04,.93,'Human UC fibroblasts: receptor availability',fontsize=11,weight='bold');fig.text(.04,.84,f'{len(patient)} paired patients; {len(m)} matched biopsies; patient-level 95% bootstrap intervals',fontsize=7)
save(fig,'Preview_Figure7_receptor_availability');print(results.to_string(index=False),flush=True)
if not (O/'epi_UC_targeted_cells.csv.gz').exists():sys.exit(0)
epi=pd.read_csv(O/'epi_UC_targeted_cells.csv.gz',low_memory=False);epi=epi.loc[epi.Treatment.eq('Pre')].copy()
mouse=json.loads((R/'output/Priority_Figure_Additions_2026-09-07/prespecified_programs.json').read_text());audit=json.loads((O/'epi_extraction_audit.json').read_text());orth=audit['ortholog_map'];programs={p:[orth[g] for g in gs if orth[g] and orth[g] in audit['present']] for p,gs in mouse.items()}
(O/'human_programs_used.json').write_text(json.dumps(programs,indent=2))
for name,gs in programs.items():epi[name]=epi[[g+'_log' for g in gs]].mean(axis=1)
ep=epi.groupby(keys)[list(programs)].mean().join(epi.groupby(keys).size().rename('n_cells')).reset_index();ep.to_csv(O/'epithelial_biopsy_scores.csv',index=False);ep=ep.loc[ep.n_cells.ge(20)]
pc=ep.groupby(['Patient','Inflammation'])[list(programs)].mean();wide=pc.unstack('Inflammation');dif=pd.DataFrame({p:wide[(p,'Inflamed')]-wide[(p,'Non_Inflamed')] for p in programs}).dropna();dif.to_csv(O/'epithelial_patient_differences.csv')
res=pd.DataFrame([{'feature':p,**paired(dif[p])} for p in programs]);res['q']=bh(res.p);res.to_csv(O/'epithelial_program_effects.csv',index=False)
fig,ax=plt.subplots(figsize=(183/25.4,94/25.4));fig.subplots_adjust(left=.38,right=.84,bottom=.19,top=.77)
forest(ax,res,list(programs),'','Inflamed − noninflamed\n(mean log-normalized RNA)',color='#CA6D69');fig.text(.04,.93,'Pretreatment UC colorectal epithelium',fontsize=11,weight='bold');fig.text(.04,.86,f'{len(dif)} paired patients; 95% patient-bootstrap intervals; six-program FDR',fontsize=7);save(fig,'Preview_Figure7_human_epithelial_programs');print(res.to_string(index=False),flush=True)
# Cross-compartment biopsy alignment; all covariates are joined explicitly, never by patient alone.
src=Path('input_data/mouse/Figure 6 TAURUS External Validation/tables/taurus_uc_fibroblast_sample_scores.tsv')
fs=pd.read_csv(src,sep='\t');print('Fib score columns',fs.columns.tolist(),flush=True)
score=[c for c in fs if 'fap' in c.lower() and 'inflamm' in c.lower()];assert len(score)==1,score
fc=fib.groupby(keys).agg(n_fib=('activated','size'),activated_FAP_fraction=('activated','mean')).reset_index();fc=fc.loc[fc.n_fib.ge(20)]
joined=ep.merge(fc,on=keys,validate='one_to_one').merge(fs[keys+[score[0]]],on=keys,validate='one_to_one');joined.to_csv(O/'matched_fibroblast_epithelial_biopsies.csv',index=False)
# Biopsy-level fixed-effects regression with patient-cluster CR1 covariance and t(G-1).
# Adjust for patient, anatomic site and binary inflammation; interpret as exploratory association.
cov=pd.get_dummies(joined[['Patient','Site','Inflammation']],drop_first=True,dtype=float);base=np.column_stack([np.ones(len(joined)),cov.to_numpy()]);out=[]
for xname in ['activated_FAP_fraction',score[0]]:
 for yname in ['Inflammatory chemokines','Barrier / junctions']:
  x=stats.zscore(joined[xname].to_numpy(float));y=stats.zscore(joined[yname].to_numpy(float));X=np.column_stack([base,x]);rank=np.linalg.matrix_rank(X);N=len(y);G=joined.Patient.nunique();assert rank>np.linalg.matrix_rank(base),'Predictor not estimable'
  inv=np.linalg.pinv(X.T@X);beta=np.linalg.lstsq(X,y,rcond=None)[0];u=y-X@beta;meat=np.zeros((X.shape[1],X.shape[1]))
  for pat in joined.Patient.unique():
   ix=joined.Patient.eq(pat).to_numpy();s=X[ix].T@u[ix];meat+=np.outer(s,s)
  assert N>rank and G>=5
  se=np.sqrt((inv@meat@inv)[-1,-1]*G/(G-1)*(N-1)/(N-rank));t=beta[-1]/se;p=2*stats.t.sf(abs(t),G-1);v=stats.t.ppf(.975,G-1)*se
  out.append(dict(predictor=xname,outcome=yname,effect=beta[-1],lo=beta[-1]-v,hi=beta[-1]+v,p=p,n_patients=G,n_biopsies=N,design_rank=rank,residual_df=N-rank))
cr=pd.DataFrame(out);cr['q']=bh(cr.p);cr.to_csv(O/'cross_compartment_adjusted_associations.csv',index=False)
fig,ax=plt.subplots(figsize=(183/25.4,89/25.4));fig.subplots_adjust(left=.43,right=.84,bottom=.22,top=.7)
labels=['FAP fraction → chemokines','FAP fraction → barrier','FAP inflammatory score → chemokines','FAP inflammatory score → barrier'];forest(ax,cr,labels,'Same-biopsy fibroblast–epithelial associations','Adjusted standardized coefficient\n(95% patient-cluster intervals)',color='#239997');fig.text(.035,.93,'Anatomically matched pretreatment UC biopsies',fontsize=11,weight='bold');fig.text(.035,.84,f'{G} patients; {N} biopsies; adjusted for patient, site and inflammation',fontsize=7);save(fig,'Preview_Figure7_matched_tissue_associations');print(cr.to_string(index=False),flush=True)
