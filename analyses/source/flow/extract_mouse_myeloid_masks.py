from pathlib import Path
import sys,json,warnings
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/packages'))
import numpy as np,flowkit as fk
warnings.filterwarnings('ignore',message='WSP references')
CACHE=ROOT/'analysis/regating/events';OUT=ROOT/'analysis/mouse_extensive'
w=fk.Workspace(str(ROOT/'analysis/regating/mouse_canonical.wsp'),load_missing_file_data=True);masks={}
for m in json.loads((ROOT/'analysis/regating/event_cache_index.json').read_text()):
    if m['experiment']!='mouse' or not m['biological']:continue
    s=fk.Sample(str(ROOT/m['file']));gs=w.get_gating_strategy(m['reference_workspace_sample']);result=gs.gate_sample(s,cache_events=False)
    full=result.get_gate_membership('f480+',gate_path='root/CELLS/s1/s2/cell1/cd11b+')
    with np.load(CACHE/(m['key']+'.npz')) as z:mask=full[z['event_indices']]
    count=next(r['canonical_count'] for r in m['baseline_comparison'] if r['gate']=='CELLS/s1/s2/cell1/cd11b+/f480+')
    assert int(mask.sum())==count,(m['sample'],mask.sum(),count)
    masks[m['sample']]=mask;print(m['sample'],count,flush=True)
np.savez_compressed(OUT/'f480_masks.npz',**masks)
