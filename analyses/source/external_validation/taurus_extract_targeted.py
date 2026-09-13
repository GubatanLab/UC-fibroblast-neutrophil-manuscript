from pathlib import Path
import sys,json
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'tmp/taurus_analysis_packages'));sys.path.insert(0,str(R/'tmp/canonical_nature_20260903/packages'))
import h5py,numpy as np,pandas as pd
from scipy.sparse import csr_matrix
O=R/'output/Figure7_Additional_TAURUS_2026-09-07';O.mkdir(exist_ok=True)
def decode(ds):
 a=ds[:]
 return np.asarray([v.decode() if isinstance(v,bytes) else v for v in a]) if a.dtype.kind in 'OS' else a
def readobs(f):
 d={}
 for k in f['obs']:
  v=f['obs'][k]
  if isinstance(v,h5py.Dataset):d[k]=decode(v)
  elif 'codes' in v:
   cats=decode(v['categories']);codes=v['codes'][:];d[k]=np.asarray([cats[i] if i>=0 else None for i in codes])
 return pd.DataFrame(d)
mode=sys.argv[1] if len(sys.argv)>1 else 'fib'
path=Path('input_data/mouse/Figure 6 TAURUS External Validation/data/fibperi_final.h5ad') if mode=='fib' else O/'data/epicolonic_final.h5ad'
mouse=json.loads((R/'output/Priority_Figure_Additions_2026-09-07/prespecified_programs.json').read_text())
# Ortholog mapping is explicit: Ly6a has no direct human ortholog and is omitted, not replaced.
special={'Ly6a':None,'Oas1a':'OAS1','Car1':'CA1','Car2':'CA2'}
orth={g:special.get(g,g.upper()) for gs in mouse.values() for g in gs}
genes=['ITGA5','ITGB1','OSMR','IL6ST'] if mode=='fib' else list(dict.fromkeys(g for g in orth.values() if g))
with h5py.File(path,'r') as f:
 obs=readobs(f);obs.to_csv(O/f'{mode}_metadata.csv',index=False);gn=decode(f['var']['_index']);print(mode,'shape',f['X'].attrs['shape'],'annotations',obs.final_analysis.value_counts().to_dict(),flush=True)
 present=[g for g in genes if g in gn];missing=[g for g in genes if g not in gn];ix=np.asarray([np.flatnonzero(gn==g)[0] for g in present]);uc=obs.Disease.eq('UC').to_numpy()
 if mode!='fib' and 'Ileum_vs_Colon' in obs:uc=uc & obs.Ileum_vs_Colon.isin(['Colon','Rectum']).to_numpy()
 n=len(obs);values=np.zeros((n,len(present)),dtype=np.float32);totals=np.zeros(n)
 x=f['X'];ip=x['indptr'][:];shape=x.attrs['shape']
 for start in range(0,n,5000):
  end=min(start+5000,n)
  if not uc[start:end].any():continue
  lo,hi=ip[start],ip[end];m=csr_matrix((x['data'][lo:hi],x['indices'][lo:hi],ip[start:end+1]-lo),shape=(end-start,shape[1]))
  values[start:end]=m[:,ix].toarray();totals[start:end]=np.asarray(m.sum(axis=1)).ravel()
  if start%50000==0:print(mode,start,'/',n,flush=True)
 assert np.all(values>=0);assert np.allclose(totals[uc],obs.loc[uc,'total_counts'].astype(float),rtol=1e-5,atol=1)
 out=obs.loc[uc].reset_index(drop=True);vv=values[uc];scale=10000/out.total_counts.to_numpy(float)
 for j,g in enumerate(present):out[g]=vv[:,j];out[g+'_log']=np.log1p(vv[:,j]*scale)
 out.to_csv(O/f'{mode}_UC_targeted_cells.csv.gz',index=False)
 (O/f'{mode}_extraction_audit.json').write_text(json.dumps({'source':str(path),'all_cells':n,'UC_cells':len(out),'present':present,'missing':missing,'total_counts_verified':True,'ortholog_map':orth if mode!='fib' else None},indent=2))
print('Done',mode,flush=True)
