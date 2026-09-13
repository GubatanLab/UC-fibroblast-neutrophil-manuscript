from pathlib import Path
import sys,json,re,itertools,math,collections
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/results'
sys.path.insert(0,str(ROOT/'analysis/packages'))
sys.stdout.reconfigure(encoding='utf-8')
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
P=pd.read_csv(OUT/'primary_measurements.csv').fillna({'channel':'','qc_flag':''})
WS={
 'co_n':'all plates result\\neutrophil_plate1_3_4_5.wsp',
 'co_f':'all plates result\\fibroblast p3-p5.wsp',
 'mono_f':'all plates result\\fibroblast only_p2_p6.wsp',
 'dss_n':'DSS_mouse\\neutrophil.wsp','dss_f':'DSS_mouse\\fibroblast.wsp'}
E=[]
def add(key,name,ws,gate,stat='frequency',channel='',denom=None,status='Exploratory; saved FlowJo counts'):
    q=P[(P.workspace==WS[ws])&(P.gate==gate)&(P.statistic==stat)&(P.channel==channel)].copy()
    for _,r in q.iterrows():
        value=r.value;den=r.parent_count;unit='percent of parent' if stat=='frequency' else 'median fluorescence, a.u.'
        if denom:
            d=P[(P.workspace==WS[ws])&(P.gate==denom)&(P.statistic=='frequency')&(P['sample']==r['sample'])]
            den=d.iloc[0]['count'];value=100*r['count']/den if den else np.nan;unit='percent of stated denominator'
        E.append(dict(endpoint_id=key,endpoint=name,sample=r['sample'],group=r.group,value=value,unit=unit,count=r['count'],denominator_count=den,denominator=denom or r.parent,gate=gate,workspace=WS[ws],statistic=stat,channel=channel,qc_flag=r.qc_flag,status=status))
CN='scatter/s1/s2/neutrophil/cxcr4+'
add('co_cxcr4','CXCR4+ / neutrophil (%)','co_n',CN)
for marker,g in [('OSM','osm'),('PADI4','padi4+'),('MX1','mx1+'),('CD177','cd177'),('MPO','mpo')]:
    add(f'co_{marker.lower()}','%s+ / CXCR4+ neutrophil (%%)'%marker,'co_n',CN+'/'+g)
for marker in ['nampt','padi4','mx1','osm','cd177','mpo']:
    add(f'co_mfi_{marker}',marker.upper()+' median FI in CXCR4+ neutrophils','co_n',CN,'Median','Comp-'+marker)
CF='scatter/s1/s2/fibroblast/pdpn+'
for marker,g in [('FAPa','fapa+'),('a5b1','a5b1+'),('FAPa a5b1 joint gate','fapa+ a5b1+')]:
    key={'FAPa':'co_fapa','a5b1':'co_a5b1','FAPa a5b1 joint gate':'co_joint'}[marker]
    add(key,marker+' gate / PDPN+ fibroblast (%)','co_f',CF+'/'+g,status='Gate definitions inconsistent; descriptive only' if key=='co_joint' else 'Exploratory; fixed co-culture panel')
for marker in ['a5b1','fapa']:add('co_f_mfi_'+marker,marker.upper()+' median FI in PDPN+ fibroblasts','co_f',CF,'Median','Comp-'+marker)
DN='CELLS/s1/s2/cell1/cd11b+/neutrophil'
add('dss_n_parent','Neutrophils / CD45/CD11b parent gate (%)','dss_n',DN)
add('dss_n_live','Neutrophils / viable scatter-selected singlets (%)','dss_n',DN,denom='CELLS/s1/s2/cell1')
for marker,g in [('CXCR4','cxcr4+'),('OSM','osm+'),('PADI4','padi4+'),('MX1','mx1+'),('CD177','cd177+ all neutrophil')]:
    add('dss_'+marker.lower(),marker+'+ / neutrophil (%)','dss_n',DN+'/'+g)
add('dss_mfi_nampt','NAMPT median FI in neutrophils','dss_n',DN,'Median','Comp-NAMPT')
DF='CELLS/S1/s2/live'
add('dss_fapa_live','FAPa gate / viable singlets (%)','dss_f',DF+'/fapa+')
add('dss_pdpn_live','PDPN+ gate / viable singlets (%)','dss_f',DF+'/pdpn+')
add('dss_fapa_pdpn','FAPa+ / PDPN+ gate (%)','dss_f',DF+'/pdpn+/fapa+')
add('dss_a5b1_pdpn','a5b1+ / PDPN+ gate (%)','dss_f',DF+'/pdpn+/a5b1+')
ep=pd.DataFrame(E)
# Per plate map, P5 wells 7-12 have no fibroblasts. Preserve in source table, exclude fibroblast inference.
ep=ep[~((ep.endpoint_id.str.startswith('co_f')|ep.endpoint_id.isin(['co_a5b1','co_joint']))&(ep.group=='N+ATN'))]
ep.to_csv(OUT/'endpoint_values.csv',index=False)
summary=ep.groupby(['endpoint_id','endpoint','group'],sort=False).agg(n=('value','count'),mean=('value','mean'),sd=('value','std'),median=('value','median'),minimum=('value','min'),maximum=('value','max'),min_gated_events=('count','min'),min_denominator=('denominator_count','min')).reset_index()
summary.to_csv(OUT/'group_summary.csv',index=False)
CO=[('N+IF','N'),('N+NF','N'),('N+IF','N+NF'),('N+IF+FK','N+IF'),('N+IF+ATN','N+IF'),('N+IF+ATN+FK','N+IF'),('N+IF+ATN+FK','N+IF+ATN'),('N+FK','N'),('N+ATN','N')]
comparisons=[]
def compare(key,a,b,family):
    q=ep[ep.endpoint_id==key];x=q[q.group==a].value.dropna().to_numpy();y=q[q.group==b].value.dropna().to_numpy()
    if len(x)<2 or len(y)<2:return
    nx,ny=len(x),len(y);vx=x.var(ddof=1);vy=y.var(ddof=1);se=math.sqrt(vx/nx+vy/ny)
    df=(vx/nx+vy/ny)**2/((vx/nx)**2/(nx-1)+(vy/ny)**2/(ny-1)) if se else np.nan
    diff=x.mean()-y.mean();t=stats.ttest_ind(x,y,equal_var=False);margin=stats.t.ppf(.975,df)*se
    z=np.r_[x,y]; comb=np.array(list(itertools.combinations(range(len(z)),nx)));sx=z[comb].sum(1);d=sx/nx-(z.sum()-sx)/ny
    pp=np.mean(np.abs(d)>=abs(diff)-1e-10)
    comparisons.append(dict(family=family,endpoint_id=key,endpoint=q.iloc[0].endpoint,group_a=a,group_b=b,n_a=nx,n_b=ny,mean_a=x.mean(),sd_a=x.std(ddof=1),mean_b=y.mean(),sd_b=y.std(ddof=1),difference_a_minus_b=diff,ci95_low=diff-margin,ci95_high=diff+margin,p_welch=t.pvalue,p_exact_permutation=pp,df_welch=df,status='Exploratory; verify saved gates before manuscript use'))
for a,b in CO:compare('co_cxcr4',a,b,'Co-culture CXCR4, 9 contrasts')
for key in ['co_osm','co_padi4','co_mx1']:
    for a,b in [CO[i] for i in [0,2,3,4,5,6]]:compare(key,a,b,'Co-culture subset markers, 18 contrasts')
for key in ['co_fapa','co_a5b1']:
    for a,b in [CO[i] for i in [2,3,4,5,6]]:compare(key,a,b,'Co-culture fibroblast markers, 10 contrasts')
for key in ['dss_n_live','dss_cxcr4','dss_osm','dss_padi4','dss_mx1','dss_fapa_live','dss_pdpn_live','dss_fapa_pdpn','dss_a5b1_pdpn']:
    for a,b in [('FD','WD'),('FG','FC')]:compare(key,a,b,'DSS, 18 contrasts')
statsdf=pd.DataFrame(comparisons)
for family,ix in statsdf.groupby('family').groups.items():
    for column in ['p_welch','p_exact_permutation']:
        vals=statsdf.loc[ix,column];order=vals.sort_values().index
        adjusted=np.minimum(1,np.maximum.accumulate(vals.loc[order].to_numpy()*np.arange(len(order),0,-1)))
        statsdf.loc[order,column+'_holm']=adjusted
statsdf.to_csv(OUT/'statistical_comparisons.csv',index=False)
print('CXCR4 comparisons')
print(statsdf[statsdf.endpoint_id=='co_cxcr4'][['group_a','group_b','difference_a_minus_b','ci95_low','ci95_high','p_welch','p_welch_holm','p_exact_permutation_holm']].round(5).to_string(index=False))
print('DSS comparisons')
print(statsdf[statsdf.family.str.startswith('DSS')][['endpoint_id','group_a','group_b','mean_a','mean_b','p_welch','p_welch_holm']].round(5).to_string(index=False))

# Quantify impossible subset relationships without assuming identical positivity thresholds.
q=P[(P.workspace==WS['co_f'])&(P.statistic=='frequency')&(P.gate.isin([CF+'/fapa+',CF+'/a5b1+',CF+'/fapa+ a5b1+']))]
piv=q.pivot(index='sample',columns='gate',values='value')
piv['joint_exceeds_either_single_gate']=piv[CF+'/fapa+ a5b1+']>piv[[CF+'/fapa+',CF+'/a5b1+']].min(axis=1)+1e-9
piv.to_csv(OUT/'fibroblast_gate_consistency.csv')
print('Joint exceeds single:',int(piv.joint_exceeds_either_single_gate.sum()),'/',len(piv))

# Detect DSS Prism misplaced rows by a multivariate fingerprint, keeping all original numbers.
src=json.loads((ROOT/'analysis/source_extract.json').read_text(encoding='utf-8'))
pn=next(p for p in src['prism'] if p['file']=='DSS_mouse\\neutrophil.pzfx')
prism_lookup={t['id']:t for t in pn['tables']}
spec=[('Table0',0,DN+'/cd177+ all neutrophil'),('Table0',1,DN+'/cd177- all neutrophils'),('Table4',0,DN+'/cxcr4+'),('Table4',1,DN+'/cxcr4-'),('Table8',0,DN+'/mx1+'),('Table8',1,DN+'/osm+'),('Table8',2,DN+'/padi4+')]
saved=P[(P.workspace==WS['dss_n'])&(P.statistic=='frequency')].pivot(index='sample',columns='gate',values='value')
matchrows=[]
for group in ['WD','FD','FC','FG']:
    for pos in range(6):
        vals=[];names=[]
        for tid,ri,g in spec:
            c=next(c for c in prism_lookup[tid]['columns'] if c['title']==group)
            try:v=float(c['subcolumns'][pos][ri]['value'])
            except (ValueError,IndexError):continue
            vals.append(v);names.append(g)
        if not vals:continue
        dist=(saved[names]-vals).abs().max(axis=1).sort_values()
        matchrows.append(dict(prism_group=group,replicate_position=pos+1,matched_sample=dist.index[0],max_difference_percentage_points=dist.iloc[0],metrics_compared=len(vals),group_mismatch=not dist.index[0].startswith(group)))
pd.DataFrame(matchrows).to_csv(OUT/'dss_prism_sample_identity_audit.csv',index=False)
print('DSS Prism mismatches:',[x for x in matchrows if x['group_mismatch']])

FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
colors={'N':'#596673','N+NF':'#3E8196','N+IF':'#BE5E35','N+IF+FK':'#D39252','N+IF+ATN':'#327E69','N+IF+ATN+FK':'#735890','N+FK':'#959FAA','N+ATN':'#99BFA9','WD':'#BD633E','FD':'#267A78','FC':'#667A91','FG':'#8D699E'}
coorder=['N','N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK','N+FK','N+ATN']
labels={'N':'N','N+NF':'N + NF','N+IF':'N + IF','N+IF+FK':'N + IF\n+ FK','N+IF+ATN':'N + IF\n+ ATN','N+IF+ATN+FK':'N + IF\n+ ATN + FK','N+FK':'N + FK','N+ATN':'N + ATN'}
def panel(ax,key,order,ylabel):
    q=ep[ep.endpoint_id==key]
    for i,g in enumerate(order):
        vals=q[q.group==g].sort_values('sample').value.dropna().to_numpy()
        jit=np.linspace(-.12,.12,len(vals))
        ax.scatter(i+jit,vals,s=31,color=colors[g],edgecolors='white',linewidths=.45,zorder=3)
        if len(vals)>1:ax.errorbar(i,vals.mean(),yerr=vals.std(ddof=1),fmt='_',markersize=19,capsize=4,color='black',lw=1.1,zorder=4)
    ax.set_xticks(range(len(order)),[labels.get(g,g) for g in order],fontsize=8)
    ax.set_ylabel(ylabel,fontsize=9);ax.grid(axis='y',color='#E4E8EC',lw=.6);ax.set_axisbelow(True)
    ax.set_ylim(bottom=0)
def save(fig,name):
    for ext in ['png','svg','pdf']:fig.savefig(FIG/f'{name}.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)
fig,axes=plt.subplots(2,2,figsize=(12,8.2))
for ax,key,title in zip(axes.flat,['co_cxcr4','co_osm','co_padi4','co_mx1'],['A  CXCR4','B  OSM within CXCR4+','C  PADI4 within CXCR4+','D  MX1 within CXCR4+']):
    panel(ax,key,coorder,'Percent of neutrophils' if key=='co_cxcr4' else 'Percent of CXCR4+ neutrophils');ax.set_title(title,loc='left',fontweight='bold')
fig.suptitle('Co-culture: saved FlowJo results',x=.055,ha='left',fontsize=16,fontweight='bold')
fig.text(.055,.02,'Independent replicates (n=6; N+NF n=5); mean ± SD. N=neutrophils; IF/NF=inflamed/non-inflamed fibroblasts; FK=FK-866; ATN=ATN-161.\nPlate 5 labels corrected from the schematic. Exploratory: validate raw gate reproduction and MX1/PADI4 channel assignments.',fontsize=9,color='#515B66')
fig.subplots_adjust(top=.91,bottom=.14,hspace=.4,wspace=.26);save(fig,'01_coculture_neutrophils')
fig,axes=plt.subplots(2,3,figsize=(12,7.8))
for ax,key,title in zip(axes.flat,['dss_n_live','dss_cxcr4','dss_osm','dss_padi4','dss_fapa_live','dss_pdpn_live'],['A  Neutrophil frequency','B  CXCR4+ neutrophils','C  OSM+ neutrophils','D  PADI4+ neutrophils','E  FAPa gate','F  PDPN+ gate']):
    panel(ax,key,['WD','FD','FC','FG'],'Percent of viable singlets' if key in ['dss_n_live','dss_fapa_live','dss_pdpn_live'] else 'Percent of neutrophils');ax.set_title(title,loc='left',fontweight='bold',fontsize=11)
fig.suptitle('DSS mouse study: sample IDs recovered from FlowJo',x=.055,ha='left',fontsize=16,fontweight='bold')
fig.text(.055,.02,'Points: individual mice; bars: mean ± SD. WD n=6, FD n=5 (FD3 absent), FC/FG n=6.\nWD=WT+DSS; FD=FAP+GCV+DSS; FC=FAP+PBS+water; FG=FAP+GCV+water. Frequencies are not absolute tissue counts.\nExploratory: saved gates, low event counts and Prism/workspace discrepancies require review.',fontsize=9,color='#515B66')
fig.subplots_adjust(top=.9,bottom=.15,hspace=.4,wspace=.32);save(fig,'02_dss_mouse')
fig,axes=plt.subplots(1,3,figsize=(12,4.5))
order=['N+NF','N+IF','N+IF+FK','N+IF+ATN','N+IF+ATN+FK']
for ax,key,title in zip(axes,['co_fapa','co_a5b1','co_joint'],['A  FAPa gate','B  α5β1 gate','C  Separate joint gate']):
    panel(ax,key,order,'Percent of PDPN+ gate');ax.set_title(title,loc='left',fontweight='bold')
fig.suptitle('Fibroblast gates in fixed co-cultures',x=.055,ha='left',fontsize=16,fontweight='bold')
fig.text(.055,.015,'Points: independent biological replicates; bars: mean ± SD. Joint gate is not a validated intersection of the single-marker gates.\nP5 wells 7–12 excluded from fibroblast interpretation because the schematic specifies neutrophils + ATN-161 only.',fontsize=9,color='#515B66')
fig.subplots_adjust(top=.85,bottom=.23,wspace=.3);save(fig,'03_coculture_fibroblast_gates')

# Small workbook payload: full source tables remain available as auditable CSV files.
payload={'sample_map':pd.read_csv(OUT/'sample_map.csv').fillna('').to_dict('records'),'endpoints':ep.fillna('').to_dict('records'),'summary':summary.fillna('').to_dict('records'),'comparisons':statsdf.fillna('').to_dict('records'),'dss_identity_audit':matchrows}
(OUT/'workbook_payload.json').write_text(json.dumps(payload,indent=2,allow_nan=False),encoding='utf-8')
