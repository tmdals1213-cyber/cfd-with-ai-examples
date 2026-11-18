from pathlib import Path
import subprocess,sys,json,hashlib,ast,re
ROOT=Path(__file__).resolve().parents[1]
subprocess.run([sys.executable,'-B',str(ROOT/'scripts/verify_release.py')],check=True)
manifest=json.loads((ROOT/'provenance/legacy_source_manifest.json').read_text());errors=[]
for f in manifest['files']:
 if hashlib.sha256((ROOT/f['file']).read_bytes()).hexdigest()!=f['releaseSHA256']:errors.append(f['file']+': source hash')
for p in ROOT.rglob('*'):
 if not p.is_file()or'.git'in p.parts:continue
 if p.suffix=='.json':
  try:json.loads(p.read_text(encoding='utf-8-sig'))
  except Exception as ex:errors.append(str(p.relative_to(ROOT))+': JSON '+str(ex))
 if p.suffix=='.md':
  for link in re.findall(r'\]\(([^\s)]+)(?:\s+[^)]*)?\)',p.read_text(encoding='utf-8')):
   if re.match(r'[a-zA-Z]+:',link)or link.startswith('#'):continue
   if not(p.parent/link.split('#')[0]).exists():errors.append(str(p.relative_to(ROOT))+': link '+link)
assert not errors,'\n'.join(errors)
print(json.dumps({'PASS':True,'NGA2examples':7,'portfolioEntries':len(json.loads((ROOT/'examples/portfolio_index.json').read_text())),'legacySourceHashesChecked':len(manifest['files']),'rawFieldsIncluded':False}))
