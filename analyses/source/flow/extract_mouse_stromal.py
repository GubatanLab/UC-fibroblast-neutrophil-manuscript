"""Corrected-alias import of the existing mouse stromal workspace, without new gates."""
from pathlib import Path
import sys,json,hashlib,re,warnings
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,pandas as pd,flowkit as fk
warnings.filterwarnings('ignore',message='WSP references')
OUT=ROOT/'analysis/mouse_extensive';OUT.mkdir(exist_ok=True)
CACHE=OUT/'stromal_events';CACHE.mkdir(exist_ok=True)
w=fk.Workspace(str(ROOT/'analysis/regating/mouse_fibroblast_canonical.wsp'),load_missing_file_data=True)
known=w.get_sample_ids(loaded_only=False)
P=pd.read_csv(ROOT/'analysis/results/primary_measurements.csv')
index=[];compares=[]
for fp in sorted((ROOT/'DSS_mouse').glob('*.fcs')):
    match=re.search(r'samples_(WD|FD|FC|FG)-?(\d+)',fp.name)
    biological=bool(match)
    if not biological and fp.name not in ['samples_FMO FAPA.fcs','samples_FMO PDPN.fcs','samples_FMO A5B1.fcs','samples_FMO NAMPT.fcs']:continue
    name=match[1]+match[2] if match else fp.stem
    key=re.sub('[^A-Za-z0-9_-]','_',name)
    dest=CACHE/(key+'.npz');meta=CACHE/(key+'.json')
    if dest.exists() and meta.exists():
        m=json.loads(meta.read_text());index.append(m);compares.extend(m['comparison']);continue
    s=fk.Sample(str(fp));ref=s.id if s.id in known else 'samples_FG1.fcs'
    gs=w.get_gating_strategy(ref);result=gs.gate_sample(s,cache_events=False)
    paths={'/'.join(list(gp)[1:]+[n]):(n,gp) for n,gp in gs.get_gate_ids()}
    desired=['CELLS','CELLS/S1','CELLS/S1/s2','CELLS/S1/s2/live','CELLS/S1/s2/live/cd45 neg']
    assert all(p in paths for p in desired),list(paths)
    masks={}
    for i,p in enumerate(desired):
        n,gp=paths[p];masks['gate_'+str(i)]=result.get_gate_membership(n,gate_path='/'.join(gp))
    comp=w.get_comp_matrix(ref);s.apply_compensation(comp);ev=s.get_events(source='comp')
    gate_counts={}
    for p,(n,gp) in paths.items():gate_counts[p]=int(result.get_gate_membership(n,gate_path='/'.join(gp)).sum())
    comparison=[]
    if biological:
        pp=P[(P.workspace=='DSS_mouse\\fibroblast.wsp')&(P['sample']==name)&(P.statistic=='frequency')]
        for r in pp.itertuples():
            if r.gate not in gate_counts:continue
            n=gate_counts[r.gate];den=s.event_count if r.parent=='root' else gate_counts[r.parent]
            comparison.append(dict(sample=name,gate=r.gate,saved_count=r.count,corrected_count=n,count_difference=n-r.count,saved_percent=r.value,corrected_percent=100*n/den if den else np.nan,difference_pp=100*n/den-r.value if den else np.nan))
    parent=masks['gate_4']
    m=dict(sample=name,group=match[1] if match else '',biological=biological,key=key,file=str(fp.relative_to(ROOT)),raw_sha256=hashlib.sha256(fp.read_bytes()).hexdigest(),channels=s.pnn_labels,markers=s.pns_labels,raw_events=s.event_count,gate_counts=gate_counts,matrix_detectors=comp.detectors,matrix_sha256=hashlib.sha256(comp.matrix.tobytes()).hexdigest(),matrix_condition_number=float(np.linalg.cond(comp.matrix)),comparison=comparison,reference=ref)
    # Preserve all CD45-negative candidates and gate indices for cross-workspace overlap audit.
    np.savez_compressed(dest,events=ev[parent],indices=np.flatnonzero(parent),live_indices=np.flatnonzero(masks['gate_3']),scatter_indices=np.flatnonzero(masks['gate_0']))
    meta.write_text(json.dumps(m,indent=2));index.append(m);compares.extend(comparison)
    print(name,'CD45-negative candidates',int(parent.sum()),flush=True)
(OUT/'stromal_index.json').write_text(json.dumps(index,indent=2))
pd.DataFrame(compares).to_csv(OUT/'stromal_corrected_vs_saved.csv',index=False)
print('Completed',len(index),'stromal samples/controls',flush=True)
