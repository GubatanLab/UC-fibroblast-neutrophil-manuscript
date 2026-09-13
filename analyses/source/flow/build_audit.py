from pathlib import Path
import sys, json, re, hashlib, csv, collections
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis'/'results'
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/'analysis'/'packages'))
import numpy as np
import pandas as pd
from lxml import etree as ET
D=json.loads((ROOT/'analysis/source_extract.json').read_text(encoding='utf-8'))
def norm(s):
    m=re.search(r'p([1-6])[-_](\d+)',s.lower())
    if m:return f'P{m[1]}-{int(m[2])}'
    m=re.search(r'(WD|FD|FC|FG)-?(\d+)',s.upper())
    if m:return f'{m[1]}{int(m[2])}'
    return s
MAP={
 (1,0):('N','Neutrophils'),(1,1):('N+FK','Neutrophils + FK-866'),
 (2,0):('NF','Non-inflamed fibroblasts'),(2,1):('IF','Inflamed fibroblasts'),
 (3,0):('N+IF','Neutrophils + inflamed fibroblasts'),(3,1):('N+IF+FK','Neutrophils + inflamed fibroblasts + FK-866'),
 (4,0):('N+IF+ATN','Neutrophils + inflamed fibroblasts + ATN-161'),(4,1):('N+IF+ATN+FK','Neutrophils + inflamed fibroblasts + ATN-161 + FK-866'),
 (5,0):('N+NF','Neutrophils + non-inflamed fibroblasts'),(5,1):('N+ATN','Neutrophils + ATN-161'),
 (6,0):('IF+FK','Inflamed fibroblasts + FK-866'),(6,1):('IF+ATN','Inflamed fibroblasts + ATN-161')}
PRIMARY=['all plates result\\neutrophil_plate1_3_4_5.wsp','all plates result\\fibroblast p3-p5.wsp','all plates result\\fibroblast only_p2_p6.wsp','DSS_mouse\\neutrophil.wsp','DSS_mouse\\fibroblast.wsp']
mapping=[]
for (pl,half),(code,label) in MAP.items():
    for well in range(1+half*6,7+half*6):
        mapping.append(dict(sample=f'P{pl}-{well}',group=code,condition=label,plate=pl,well=well,replicate_position=(well-1)%6+1,independence='six independent biological replicates, user confirmed; pairing unconfirmed',mapping_source='Co-culture Experiment Schematic, embedded plate map',note='Prism copies corrected to schematic: no inhibitors; P5-5 absent' if pl==5 and half==0 else ''))
mp={x['sample']:x for x in mapping}
for group,n in [('WD',6),('FD',6),('FC',6),('FG',6)]:
    for i in range(1,n+1):
        sid=f'{group}{i}'
        mp[sid]=dict(sample=sid,group=group,condition={'WD':'WT + DSS','FD':'FAP + GCV + DSS','FC':'FAP + PBS + water','FG':'FAP + GCV + water'}[group],plate='',well='',replicate_position=i,independence='mouse ID; presumed independent animal',mapping_source='DSS_mouse/read me.docx',note='Explicitly skipped in source; no FCS present' if sid=='FD3' else '')
mapping=list(mp.values())
rows=[]; samples=[]; gates=[]
ns={'g':'http://www.isac-net.org/std/Gating-ML/v2.0/gating','d':'http://www.isac-net.org/std/Gating-ML/v2.0/datatypes'}
for w in D['flowjo']:
    for s in w['samples']:
        sid=norm(s['node'].get('name',''))
        if sid not in mp:continue
        count=float(s['node'].get('count','nan'))
        counts={p['path']:float(p['attrs']['count']) for p in s['populations']}
        meta=dict(workspace=w['file'],primary=w['file'] in PRIMARY,sample=sid,group=mp[sid]['group'],sample_name=s['node']['name'],acquisition_date=s['keywords'].get('$DATE'),acquired_events=count,workspace_modified=w['attrs'].get('modDate'))
        samples.append(meta)
        for p in s['populations']:
            path=p['path']; par=path.rsplit('/',1)[0] if '/' in path else 'root'
            n=float(p['attrs']['count']); den=counts.get(par,count)
            flags=[]
            if n<20:flags.append('fewer than 20 gated events')
            if den<100:flags.append('parent fewer than 100 events')
            if den==0:flags.append('empty parent')
            row={**meta,'gate':path,'parent':par,'count':int(n),'parent_count':int(den),'percent_parent':100*n/den if den>0 else None,'percent_acquired':100*n/count if count>0 else None,'qc_flag':'; '.join(flags)}
            rows.append({**row,'statistic':'frequency','channel':'','value':row['percent_parent'],'unit':'percent of parent'})
            for st in p['stats']:
                try: val=float(st['value'])
                except (ValueError,KeyError):continue
                if n==0:val=None
                rows.append({**row,'statistic':st.get('name'),'channel':st.get('id'),'value':val,'unit':'arbitrary fluorescence units'})
            if p['gate_xml']:
                r=ET.fromstring(p['gate_xml'].encode()); ge=next(iter(r))
                dims=[{'channel':dim.find('d:fcs-dimension',ns).get('{'+ns['d']+'}name'),'min':dim.get('{'+ns['g']+'}min'),'max':dim.get('{'+ns['g']+'}max')} for dim in ge.findall('g:dimension',ns)]
                gates.append({**meta,'gate':path,'gate_type':ET.QName(ge).localname,'dimensions':json.dumps(dims),'geometry_sha256':hashlib.sha256(ET.tostring(ge,method='c14n')).hexdigest()})
allrows=pd.DataFrame(rows); allrows.to_csv(OUT/'all_workspace_measurements.csv',index=False)
pd.DataFrame(mapping).to_csv(OUT/'sample_map.csv',index=False)
pd.DataFrame(samples).to_csv(OUT/'workspace_sample_inventory.csv',index=False)
pd.DataFrame(gates).to_csv(OUT/'gate_definitions.csv',index=False)
primary=allrows[allrows.primary].copy(); primary.to_csv(OUT/'primary_measurements.csv',index=False)
prismrows=[]
for p in D['prism']:
    for t in p['tables']:
        for ci,c in enumerate(t['columns']):
            for si,sc in enumerate(c['subcolumns']):
                for ri,v in enumerate(sc):
                    try:val=float(v['value'])
                    except ValueError:val=None
                    prismrows.append(dict(file=p['file'],table_id=t['id'],table=t['title'],column=c['title'],column_position=ci+1,row=t['rows'][ri] if ri<len(t['rows']) else '',row_position=ri+1,subcolumn_position=si+1,value=val,excluded=v['attrs'].get('Excluded','0'),attributes=json.dumps(v['attrs'])))
pd.DataFrame(prismrows).to_csv(OUT/'prism_source_values.csv',index=False)
xrows=[]
for b in D['workbooks']:
    for sh in b['sheets']:
        headers=[]
        for ri,row in enumerate(sh['rows']):
            if row and row[0]=='Sample:':headers=row;continue
            for ci,v in enumerate(row[1:],1):
                if ci<len(headers) and headers[ci] and isinstance(v,(int,float)):
                    xrows.append(dict(file=b['file'],sheet=sh['name'],row_number=ri+1,column_number=ci+1,sample=norm(headers[ci]),metric=row[0],value=v))
pd.DataFrame(xrows).to_csv(OUT/'excel_source_values.csv',index=False)

# Compare each Excel value to the matching saved FlowJo statistic, retaining source/version.
reconcile=[]
for x in xrows:
    if ' | ' not in str(x['metric']):continue
    gate,metric=x['metric'].split(' | ',1)
    channel=''; stat='frequency'
    if metric.startswith('Median'):
        stat='Median';channel=metric[metric.find('(')+1:metric.rfind(')')].split(' :: ')[0]
        if not channel.startswith('Comp-'):channel='Comp-'+channel
    cand=allrows[(allrows['sample']==x['sample'])&(allrows.gate==gate)&(allrows.statistic==stat)&(allrows.channel==channel)]
    for _,c in cand.iterrows():
        tol=max(0.5*10**(np.floor(np.log10(abs(x['value'])))-2),0.000001) if x['value'] else 0.000001
        reconcile.append({**x,'workspace':c.workspace,'workspace_value':c.value,'difference':x['value']-c.value,'matches_3_significant_figures':bool(abs(x['value']-c.value)<=tol)})
pd.DataFrame(reconcile).to_csv(OUT/'excel_workspace_reconciliation.csv',index=False)

# Finite, positive event counts are descriptive. No automatic sample exclusions.
pri_counts=primary[primary.statistic=='frequency']
print('Primary biological sample counts:',pri_counts.groupby('workspace')['sample'].nunique().to_dict())
print('Key group means from saved workspace')
for ws,suffix in [(PRIMARY[0],'/neutrophil/cxcr4+'),(PRIMARY[1],'/pdpn+/fapa+'),(PRIMARY[3],'/neutrophil')]:
    q=pri_counts[(pri_counts.workspace==ws)&pri_counts.gate.str.endswith(suffix)]
    print(ws,suffix,q.groupby('group')['value'].agg(['count','mean','std']).round(3).to_dict('index'))
print('Excel reconciliation summary:')
rr=pd.DataFrame(reconcile)
print(rr.groupby(['file','workspace']).matches_3_significant_figures.agg(['size','mean']).round(3).to_string())
