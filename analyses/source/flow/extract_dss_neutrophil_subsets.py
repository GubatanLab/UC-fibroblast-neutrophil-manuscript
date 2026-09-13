from pathlib import Path
import sys,json,warnings
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,flowkit as fk
warnings.filterwarnings('ignore',message='WSP references')
OUT=ROOT/'analysis/dss_neutrophil_subsets';OUT.mkdir(exist_ok=True)
w=fk.Workspace(str(ROOT/'analysis/regating/mouse_canonical.wsp'),load_missing_file_data=True)
base='CELLS/s1/s2/cell1/cd11b+/neutrophil'
meta=[m for m in json.loads((ROOT/'analysis/regating/event_cache_index.json').read_text()) if m['experiment']=='mouse' and m['biological'] and m['group'] in ['WD','FD']]
for m in meta:
    s=fk.Sample(str(ROOT/m['file']));gs=w.get_gating_strategy(m['reference_workspace_sample']);result=gs.gate_sample(s,cache_events=False)
    with np.load(ROOT/'analysis/regating/events'/(m['key']+'.npz')) as z:
        inds=z['event_indices'];parent=z['old_5'];masks={}
        for row in m['baseline_comparison']:
            g=row['gate']
            if g==base or g.startswith(base+'/'):
                path,name=g.rsplit('/',1);mask=result.get_gate_membership(name,gate_path='root/'+path)[inds]
                assert int(mask.sum())==row['canonical_count'],(m['sample'],g)
                assert not np.any(mask&~parent)
                masks[g]=mask[parent]
    np.savez_compressed(OUT/(m['sample']+'.npz'),**masks)
    print(m['sample'],len(masks),'verified gates',flush=True)
