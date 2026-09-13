from pathlib import Path
import sys,json,itertools
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'tmp/canonical_nature_20260903/packages'))
import numpy as np,pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=R/'output/Priority_Figure_Additions_2026-09-07';O.mkdir(exist_ok=True)
Q=Path('input_data/mouse')
rng=np.random.default_rng(20260907)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':.7,'svg.fonttype':'none'})
PUR='#7651A6';TEAL='#168C88';ORANGE='#D97724';INK='#30343A';GREY='#9AA4AE'
def clean(ax):
 ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3,width=.7)
def save(fig,name):
 fig.savefig(O/(name+'.png'),dpi=320);fig.savefig(O/(name+'.svg'));plt.close(fig)
def bh(p):
 p=np.asarray(p);order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
def diffboot(a,b,func=np.mean,n=10000):
 boot=func(rng.choice(a,(n,len(a))),axis=1)-func(rng.choice(b,(n,len(b))),axis=1)
 return func(a)-func(b),*np.quantile(boot,[.025,.975])

# Investigator-defined panels specified before scoring; no outcome-based gene selection.
programs={
 'Inflammatory chemokines':['Il1b','Il6','Cxcl1','Cxcl2','Cxcl5','Ccl20'],
 'Interferon response':['Mx1','Isg15','Ifit1','Ifit2','Ifit3','Oas1a','Stat1'],
 'Barrier / junctions':['Tjp1','Ocln','Cldn3','Cldn4','Cldn7','Cdh1'],
 'Mucus / secretory':['Muc2','Tff3','Fcgbp','Clca1','Zg16','Agr2'],
 'Injury-associated regeneration':['Clu','Ly6a','Anxa1','Anxa2','Krt8','Krt18','Tacstd2'],
 'Absorptive differentiation':['Krt20','Slc26a3','Aqp8','Car1','Car2','Alpi','Fabp1']}
(O/'prespecified_programs.json').write_text(json.dumps(programs,indent=2))
inv=pd.read_csv(O/'epithelial_pseudobulk_inventory.csv');units=pd.read_csv(O/'mouse_units.csv')
units['group']=np.where(units.source_batch=='M1','DSS',np.where(units.source_batch=='M3','FAP ablation','alpha5beta1 blockade'))
inv=inv.merge(units,on='sample_uid',validate='many_to_one');counts=pd.read_csv(O/'epithelial_pseudobulk_counts.csv.gz',index_col=0);expr=pd.read_csv(O/'epithelial_logcpm.csv.gz',index_col=0)
allkeys=inv.loc[(inv.scope=='All epithelial')&(inv.n_cells>=20),'key']
coverage=[];used={}
for name,genes in programs.items():
 used[name]=[]
 for g in genes:
  n=int((counts.loc[g,allkeys]>0).sum()) if g in counts.index else 0
  eligible=g in expr.index and n>=3
  coverage.append(dict(program=name,gene=g,in_assay=g in expr.index,mice_detected=n,used=eligible))
  if eligible:used[name].append(g)
 assert len(used[name])>=max(3,len(genes)/2),(name,used[name])
pd.DataFrame(coverage).to_csv(O/'epithelial_gene_coverage.csv',index=False)
scores=[]
for row in inv.itertuples():
 if row.n_cells<20:continue
 for name,genes in used.items():
  score=expr.loc[genes,row.key].mean()
  if pd.notna(score):scores.append(dict(key=row.key,sample_uid=row.sample_uid,group=row.group,scope=row.scope,state=row.state,n_cells=row.n_cells,program=name,score=score))
scores=pd.DataFrame(scores);scores.to_csv(O/'epithelial_program_scores.csv',index=False)
effects=[]
for (scope,state,name),d in scores.groupby(['scope','state','program'],sort=False):
 a=d.loc[d.group=='DSS','score'].to_numpy()
 for group in ['FAP ablation','alpha5beta1 blockade']:
  b=d.loc[d.group==group,'score'].to_numpy()
  if min(len(a),len(b))<3:continue
  eff,lo,hi=diffboot(b,a)
  p=stats.mannwhitneyu(b,a,alternative='two-sided',method='auto').pvalue
  effects.append(dict(scope=scope,state=state,program=name,group=group,n_DSS=len(a),n_treatment=len(b),effect=eff,low=lo,high=hi,p=p))
effects=pd.DataFrame(effects);effects['q']=effects.groupby('scope',group_keys=False).p.transform(lambda x:bh(x))
effects.to_csv(O/'epithelial_program_effects.csv',index=False)
main=effects[effects.scope=='All epithelial']
fig,ax=plt.subplots(figsize=(7.2,3.6));fig.subplots_adjust(left=.32,right=.91,top=.77,bottom=.19)
for j,(group,col) in enumerate([('FAP ablation',PUR),('alpha5beta1 blockade',TEAL)]):
 d=main[main.group==group].set_index('program').loc[list(programs)];y=np.arange(6)+(j-.5)*.24
 ax.errorbar(d.effect,y,xerr=[d.effect-d.low,d.high-d.effect],fmt='o',ms=4,lw=1,color=col,label=group.replace('alpha5beta1',r'$\alpha5\beta1$'))
 for yy,rr in zip(y,d.itertuples()):ax.text(1.015,yy,'<0.001' if rr.q<.001 else f'{rr.q:.3f}',transform=ax.get_yaxis_transform(),va='center',fontsize=6.5,color=col)
ax.axvline(0,c=GREY,lw=.7);ax.set_yticks(range(6),list(programs));ax.invert_yaxis();ax.set_xlabel('Intervention minus DSS (mean program log2 CPM)');clean(ax)
ax.text(1.015,1.06,'q',transform=ax.transAxes,fontsize=7);ax.legend(frameon=False,loc='lower left',bbox_to_anchor=(-.02,1.02),ncol=2,fontsize=7)
fig.text(.025,.955,'Figure 6 candidate: epithelial functional responses',weight='bold',fontsize=11)
fig.text(.025,.88,'DSS n = 6; FAP ablation n = 4; blockade n = 12. Biological-mouse means; 95% bootstrap intervals.',fontsize=8)
fig.text(.025,.035,'Exploratory gene panels; total-epithelial scores can reflect state composition. Treatment and batch are aligned.',fontsize=7)
save(fig,'Candidate_F6_Epithelial_Responses')

# State-resolved expression versus compositional changes; all eligible states, no top-effect selection.
st=effects[effects.scope=='Level 3'];eligible_states=list(st.state.unique())
if len(eligible_states):
 fig,axs=plt.subplots(1,2,figsize=(7.2,max(4.2,len(eligible_states)*.34+1.8)));fig.subplots_adjust(left=.37,right=.93,top=.79,bottom=.30,wspace=.20)
 lim=max(1,float(np.nanmax(abs(st.effect))))
 for j,group in enumerate(['FAP ablation','alpha5beta1 blockade']):
  mat=st[st.group==group].pivot(index='state',columns='program',values='effect').reindex(index=eligible_states,columns=programs)
  masked=np.ma.masked_invalid(mat.to_numpy());cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#E8EAED')
  im=axs[j].imshow(masked,cmap=cmap,vmin=-lim,vmax=lim,aspect='auto')
  axs[j].set_xticks(range(6),['Inflam.','IFN','Barrier','Mucus','Regen.','Absorp.'],rotation=55,ha='right')
  axs[j].set_yticks(range(len(eligible_states)),[s.replace('Epithelial | ','') for s in eligible_states] if j==0 else [])
  axs[j].set_title('FAP ablation' if j==0 else r'$\alpha5\beta1$ blockade')
  for ii,state in enumerate(eligible_states):
   for kk,name in enumerate(programs):
    r=st[(st.state==state)&(st.program==name)&(st.group==group)]
    if len(r) and r.iloc[0].q<.05:axs[j].text(kk,ii,'*',ha='center',va='center',fontsize=10,color='white' if abs(r.iloc[0].effect)>lim*.55 else 'black')
  axs[j].spines[:].set_visible(False);axs[j].tick_params(length=0)
 cax=fig.add_axes([.38,.12,.40,.025]);fig.colorbar(im,cax=cax,orientation='horizontal',label='Intervention minus DSS (program log2 CPM)')
 fig.text(.025,.955,'Figure 6 supporting candidate: within-state expression',weight='bold',fontsize=11)
 fig.text(.025,.875,'At least 20 cells per mouse and 3 mice per group. Gray: insufficient support; * q < 0.05.',fontsize=7.5)
 save(fig,'Candidate_F6_Within_State_Responses')

# Priority 4: joint measured programs, without a composite benefit score or selective gene removal.
ds=scores[scores.scope=='All epithelial'];fig,axs=plt.subplots(1,3,figsize=(7.2,3.1));fig.subplots_adjust(left=.09,right=.98,top=.76,bottom=.23,wspace=.38)
for ax,name in zip(axs,['Inflammatory chemokines','Barrier / junctions','Mucus / secretory']):
 for j,(group,col) in enumerate([('DSS',INK),('FAP ablation',PUR),('alpha5beta1 blockade',TEAL)]):
  vals=ds[(ds.program==name)&(ds.group==group)].score.to_numpy();ax.scatter(j+rng.uniform(-.12,.12,len(vals)),vals,s=18,c=col,edgecolor='white',lw=.3);ax.plot([j-.23,j+.23],[vals.mean()]*2,color=col,lw=2)
 ax.set_xticks(range(3),['DSS','FAP\nablation',r'$\alpha5\beta1$'+'\nblockade']);ax.set_title(name,fontsize=8);clean(ax)
axs[0].set_ylabel('Mean program log2 CPM')
fig.text(.025,.95,'Figure 6 alternative: inflammation and epithelial maintenance',weight='bold',fontsize=10)
fig.text(.025,.84,'DSS n = 6; ablation n = 4; blockade n = 12. Each dot is a mouse; bars are means.',fontsize=7.5)
fig.text(.025,.045,'Alternative view of the same epithelial scores, not an independent validation. No composite benefit score.',fontsize=7)
save(fig,'Candidate_F6_Inflammation_Maintenance')

# TAURUS: primary between-outcome delta comparison and outcome-adjusted change association.
t=pd.read_csv(R/'Nature_Style_Redrawn_Noncanonical/source_data/figure6/taurus_uc_patient_paired_deltas.tsv',sep='\t')
d=t[t.feature=='neutrophil recruitment'].copy();assert not d.Patient.duplicated().any()
a=d.loc[d.Remission_status=='Remission','delta'].to_numpy();b=d.loc[d.Remission_status=='Non_Remission','delta'].to_numpy()
eff,lo,hi=diffboot(a,b,np.median);p=stats.mannwhitneyu(a,b,method='exact').pvalue
wide=t[t.feature.isin(['activated FAP fibroblast fraction','neutrophil recruitment'])].pivot(index=['Patient','Remission_status'],columns='feature',values='delta').dropna().reset_index()
x=wide['activated FAP fibroblast fraction'].to_numpy()*100;y=wide['neutrophil recruitment'].to_numpy();grp=(wide.Remission_status=='Remission').astype(int).to_numpy()
X=np.column_stack([np.ones(len(x)),grp]);xr=stats.rankdata(x);yr=stats.rankdata(y);rx=xr-X@np.linalg.lstsq(X,xr,rcond=None)[0];ry=yr-X@np.linalg.lstsq(X,yr,rcond=None)[0]
partial=float(np.corrcoef(rx,ry)[0,1]);pooled=float(stats.spearmanr(x,y).statistic)
per=[]
for i in range(19999):
 z=ry.copy()
 for g in [0,1]:z[grp==g]=rng.permutation(z[grp==g])
 per.append(np.corrcoef(rx,z)[0,1])
pp=(1+np.sum(np.abs(per)>=abs(partial)))/20000
qq=bh([p,pp]);results=[dict(test='Recruitment delta: remission minus nonremission',effect=eff,low=lo,high=hi,p=p,q=qq[0],n_remission=len(a),n_nonremission=len(b)),dict(test='Outcome-adjusted rank association',effect=partial,p=pp,q=qq[1],n_remission=int(grp.sum()),n_nonremission=int((1-grp).sum()))]
pd.DataFrame(results).to_csv(O/'taurus_primary_results.csv',index=False);wide.to_csv(O/'taurus_common_patient_deltas.csv',index=False);d.to_csv(O/'taurus_recruitment_patient_values.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(7.2,3.4),gridspec_kw={'width_ratios':[1,1,1.2]});fig.subplots_adjust(left=.09,right=.98,top=.76,bottom=.27,wspace=.5)
for ax,status,col in zip(axs[:2],['Remission','Non_Remission'],[TEAL,ORANGE]):
 sub=d[d.Remission_status==status]
 for r in sub.itertuples():ax.plot([0,1],[r.pre,r.post],'-o',c=col,alpha=.45,lw=.8,ms=3)
 ax.plot([0,1],[sub.pre.median(),sub.post.median()],'-o',c=col,lw=2.2,ms=5);ax.set_xticks([0,1],['Pre','Post']);ax.set_title(f'{status.replace("_"," ")} (n = {len(sub)})',fontsize=8);clean(ax)
axs[0].set_ylabel('Fibroblast recruitment score');ylim=(min(d.pre.min(),d.post.min())-.04,max(d.pre.max(),d.post.max())+.05)
for ax in axs[:2]:ax.set_ylim(ylim)
ax=axs[2]
for status,col in [('Remission',TEAL),('Non_Remission',ORANGE)]:
 sel=wide.Remission_status==status;ax.scatter(x[sel],y[sel],s=24,c=col,label=status.replace('_',' '),edgecolor='white',lw=.4)
ax.axvline(0,c=GREY,lw=.7);ax.axhline(0,c=GREY,lw=.7);ax.set_xlabel('Activated FAP+ change\n(percentage points)');ax.set_ylabel('Recruitment-score change');ax.set_title(f'Common patients (n = {len(wide)})',fontsize=8);clean(ax)
fig.text(.025,.95,'Figure 7 candidate: recruitment activity after adalimumab',weight='bold',fontsize=11)
fig.text(.025,.85,f'Between-outcome change: median difference {eff:.3f}; q = {qq[0]:.3f}.',fontsize=8)
fig.text(.025,.06,f'Common-patient association: pooled Spearman r = {pooled:.2f}; outcome-adjusted r = {partial:.2f}, q = {qq[1]:.3f}.\nSite-matched patient summaries; no neutrophil measurements or response-prediction claim.',fontsize=7)
save(fig,'Candidate_F7_Longitudinal_Recruitment')

# Spatial functional extension: independent outcome genes, existing approved sections only.
D=Q/'Figure 6 Spatial External Validation/data/SCP3818'
identity=['CSF3R','CXCR1','CXCR2','FCGR3B','CEACAM8','MPO'];outcomes=['OSM','IL1B','TNF'];genes=['FAP','PDPN','THY1']+identity+outcomes
spatial_results=[];curve=[];coverage=[];spatial_frames=[]
for stem,section in [('UC1_inflamed','Inflamed'),('UC1_less_inflamed','Less inflamed')]:
 cluster=json.loads((D/f'{stem}_cluster.json').read_text())['data'];df=pd.DataFrame({'cell_id':list(map(str,cluster['cells'])),'x':cluster['x'],'y':cluster['y'],'cell_type':cluster['annotations']}).set_index('cell_id')
 for g in genes:
  gg=pd.read_csv(D/f'{stem}_{g}.tsv',sep='\t',dtype={'cell_id':str}).set_index('cell_id');assert gg.index.is_unique and df.index.isin(gg.index).all();df[g]=gg.expression.reindex(df.index)
 xy=df[['x','y']].to_numpy();tree=cKDTree(xy);spacing=np.median(tree.query(xy,k=2)[0][:,1]);neut=(df.cell_type=='Myeloid')&((df[identity]>0).sum(axis=1)>=2);active=(df.cell_type=='Stromal')&(df.FAP>0)&((df.PDPN>0)|(df.THY1>0));other=(df.cell_type=='Stromal')&~active
 n=df[neut].copy();nxy=n[['x','y']].to_numpy();n['distance_fap']=cKDTree(xy[active]).query(nxy)[0]/spacing;n['distance_other']=cKDTree(xy[other]).query(nxy)[0]/spacing
 distances,ix=tree.query(nxy,k=31);n['epithelial_fraction']=(df.cell_type.to_numpy()[ix[:,1:]]=='Epithelial').mean(axis=1);n['local_density']=30/np.maximum(distances[:,-1]/spacing,1e-6)**2
 for g in outcomes:coverage.append(dict(section=section,gene=g,n_neutrophil_like=len(n),n_detected=int((n[g]>0).sum()),fraction_detected=float((n[g]>0).mean())))
 ranked=n[outcomes].rank(pct=True);n['Feedback']=ranked.mean(axis=1)
 # Spatial tile resampling is within-section uncertainty, not a donor-level CI.
 coords=(nxy-nxy.min(axis=0))/(np.ptp(nxy,axis=0)+1e-12);tiles=np.minimum((coords*4).astype(int),3);block=tiles[:,0]*4+tiles[:,1];n['tile']=block;n['section']=section
 cov=np.column_stack([np.log1p(n.distance_fap),np.log1p(n.distance_other),n.epithelial_fraction,np.log1p(n.local_density),coords,coords[:,0]**2,coords[:,1]**2,coords[:,0]*coords[:,1]])
 cov=(cov-cov.mean(axis=0))/cov.std(axis=0);design=np.column_stack([np.ones(len(n)),cov]);blocks=np.unique(block)
 for gene in ['Feedback']+outcomes:
  yy=n.Feedback.to_numpy() if gene=='Feedback' else ranked[gene].to_numpy();sd=yy.std();yy=(yy-yy.mean())/sd if sd>0 else np.zeros(len(yy))
  beta=np.linalg.lstsq(design,yy,rcond=None)[0][1];boots=[]
  for i in range(1000):
   take=np.concatenate([np.flatnonzero(block==b) for b in rng.choice(blocks,len(blocks),replace=True)])
   if np.linalg.matrix_rank(design[take])==design.shape[1]:boots.append(np.linalg.lstsq(design[take],yy[take],rcond=None)[0][1])
  low,high=np.quantile(boots,[.025,.975]);spatial_results.append(dict(section=section,outcome=gene,n_cells=len(n),n_active=int(active.sum()),n_tiles=len(blocks),coefficient=beta,low=low,high=high,bootstrap_valid=len(boots),condition_number=float(np.linalg.cond(design))))
 edges=[0,2,4,8,16,np.inf];n['bin']=pd.cut(n.distance_fap,edges,right=False,labels=['0-2','2-4','4-8','8-16','>16'])
 for label,sub in n.groupby('bin',observed=False):curve.append(dict(section=section,bin=label,n_cells=len(sub),score=sub.Feedback.mean() if len(sub)>=10 else np.nan))
 spatial_frames.append(n)
pd.DataFrame(coverage).to_csv(O/'spatial_outcome_gene_coverage.csv',index=False);res=pd.DataFrame(spatial_results);res.to_csv(O/'spatial_adjusted_coefficients.csv',index=False);cur=pd.DataFrame(curve);cur.to_csv(O/'spatial_distance_function_curve.csv',index=False);pd.concat(spatial_frames).to_csv(O/'spatial_function_cell_values.csv.gz')
fig,axs=plt.subplots(1,2,figsize=(7.2,3.6));fig.subplots_adjust(left=.10,right=.96,top=.76,bottom=.26,wspace=.48)
covtab=pd.DataFrame(coverage)
for j,(section,col) in enumerate([('Inflamed',ORANGE),('Less inflamed',TEAL)]):
 sub=covtab[covtab.section==section].set_index('gene').loc[outcomes];y=np.arange(3)+(j-.5)*.24
 axs[0].barh(y,sub.fraction_detected*100,height=.22,color=col,label=section,alpha=.85)
 for yy,r in zip(y,sub.itertuples()):axs[0].text(r.fraction_detected*100+.35,yy,f'{r.n_detected}/{r.n_neutrophil_like}',va='center',fontsize=6)
 for k,g in enumerate(outcomes):
  yy=k+(j-.5)*.2;r=res[(res.section==section)&(res.outcome==g)].iloc[0]
  if sub.loc[g,'n_detected']>=10:axs[1].errorbar(r.coefficient,yy,xerr=[[r.coefficient-r.low],[r.high-r.coefficient]],fmt='o',c=col,ms=4,lw=1)
  else:axs[1].text(.04,yy,'Too sparse',color=col,fontsize=6.5,va='center')
axs[0].set_yticks(range(3),outcomes);axs[0].invert_yaxis();axs[0].set_xlim(0,24);axs[0].set_xlabel('Neutrophil-like cells with detection (%)');axs[0].legend(frameon=False,fontsize=7,loc='lower right');axs[0].set_title('Outcome-gene coverage',fontsize=9)
axs[1].set_yticks(range(3),outcomes);axs[1].set_ylim(2.5,-.5);axs[1].axvline(0,c=GREY,lw=.7);axs[1].set_xlim(-.21,.13);axs[1].set_xlabel('Adjusted standardized distance coefficient');axs[1].set_title('Negative: higher activity nearer FAP+',fontsize=8)
for ax in axs:clean(ax)
fig.text(.025,.95,'Figure 7 exploratory analysis: spatial functional feasibility',weight='bold',fontsize=10.5)
fig.text(.025,.85,'Two sections from one UC donor. Coefficients displayed only for genes detected in at least 10 cells.',fontsize=7.5)
fig.text(.025,.065,'Adjusted for other-stromal distance, epithelial neighborhood, density and spatial coordinates.\nIntervals: spatial-tile bootstrap, not patient-level CIs. TNF coverage is inadequate; full RNA depth unavailable.',fontsize=7)
save(fig,'Candidate_F7_Spatial_Function')
print('Epithelial effects\n',main.to_string(index=False));print('TAURUS\n',pd.DataFrame(results).to_string(index=False));print('SPATIAL\n',res.to_string(index=False));print('COVERAGE\n',pd.DataFrame(coverage).to_string(index=False))
