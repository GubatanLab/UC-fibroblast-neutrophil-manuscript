"""Validate deposited file integrity, script syntax and figure completeness."""
from pathlib import Path
import ast,csv,hashlib,json
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((ROOT/'metadata/file_manifest.csv').open(encoding='utf-8')))
errors=[]
for row in rows:
    p=ROOT/row['path']
    if not p.is_file():errors.append('Missing '+row['path']);continue
    if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:errors.append('Changed '+row['path'])
for p in (ROOT/'analyses').rglob('*.py'):
    try:ast.parse(p.read_text(encoding='utf-8-sig'))
    except SyntaxError as e:errors.append(str(e))
for p in (ROOT/'figures/main').glob('*.pdf'):
    if len(PdfReader(p).pages)!=1:errors.append('Expected one page '+p.name)
for p in (ROOT/'figures/extended_data').glob('*.pdf'):
    if len(PdfReader(p).pages)!=1:errors.append('Expected one page '+p.name)
for p in (ROOT/'figures/supplementary').glob('*.pdf'):
    if len(PdfReader(p).pages)!=1:errors.append('Expected one page '+p.name)
if len(list((ROOT/'figures/main').glob('*.pdf')))!=7:errors.append('Expected seven main figures')
if len(list((ROOT/'figures/extended_data').glob('*.pdf')))!=10:errors.append('Expected ten extended data figures')
if len(list((ROOT/'figures/supplementary').glob('*.pdf')))!=12:errors.append('Expected twelve supplementary figures')
for filename,n in [('Main_Figures_1_to_7.pdf',7),('Supplementary_Figures_1_to_12.pdf',12),('Extended_Data_Figures_1_to_10.pdf',10)]:
    if len(PdfReader(ROOT/'figures'/filename).pages)!=n:errors.append('Wrong page count '+filename)
for row in rows:
    p=Path(row['path'])
    if p.suffix.lower() in {'.rds','.h5ad','.fcs','.qptiff','.bam','.cram','.pem','.key'}:errors.append('Unsupported deposit type '+str(p))
    if int(row['bytes'])>=100*1024*1024:errors.append('Oversized Git file '+str(p))
report={'status':'passed' if not errors else 'failed','files_checked':len(rows),'main_figures':7,'supplementary_figures':12,'errors':errors}
(ROOT/'metadata/validation_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
raise SystemExit(bool(errors))
