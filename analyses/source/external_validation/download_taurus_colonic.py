from pathlib import Path
import urllib.request,json,hashlib,time
R=Path(__file__).resolve().parents[1];O=R/'output/Figure7_Additional_TAURUS_2026-09-07';O.mkdir(exist_ok=True);D=O/'data';D.mkdir(exist_ok=True)
meta=json.load(urllib.request.urlopen('https://zenodo.org/api/records/13768607',timeout=40));(D/'zenodo_record.json').write_text(json.dumps(meta,indent=2))
entry=next(x for x in meta['files'] if x['key']=='epicolonic_final.h5ad');dest=D/entry['key'];part=D/(entry['key']+'.part')
print('Downloading approved TAURUS colonic epithelial file:',entry['size'],'bytes',flush=True)
if not dest.exists():
 with urllib.request.urlopen(entry['links']['self'],timeout=60) as resp,part.open('wb') as out:
  n=0;last=time.time()
  while True:
   chunk=resp.read(8*1024*1024)
   if not chunk:break
   out.write(chunk);n+=len(chunk)
   if time.time()-last>20:print(round(n/1e9,2),'GB',round(100*n/entry['size'],1),'%',flush=True);last=time.time()
 assert part.stat().st_size==entry['size'];part.rename(dest)
h=hashlib.md5()
with dest.open('rb') as inp:
 for chunk in iter(lambda:inp.read(8*1024*1024),b''):h.update(chunk)
assert 'md5:'+h.hexdigest()==entry['checksum'],(h.hexdigest(),entry['checksum'])
(D/'download_verification.json').write_text(json.dumps({'file':dest.name,'size':dest.stat().st_size,'checksum':entry['checksum'],'verified':True},indent=2));print('Verified download',flush=True)
