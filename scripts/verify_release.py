#!/usr/bin/env python3
"""Static release checks; no network, no solver, no secret recovery."""
import ast,hashlib,json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ignore={'.git','vendor','runs','build','.cache','__pycache__'}
files=[p for p in ROOT.rglob('*')if p.is_file()and not(ignore&set(p.relative_to(ROOT).parts))]
problems=[]
for p in files:
 rel=p.relative_to(ROOT).as_posix()
 if p.is_symlink():problems.append(rel+': symlink')
 if p.stat().st_size>50*1024*1024:problems.append(rel+': >50 MiB')
 if p.suffix in {'.bin','.exe','.pem','.key','.o','.mod'}:problems.append(rel+': runtime/private artifact')
 if p.suffix=='.py':
  try:ast.parse(p.read_text(encoding='utf-8-sig'))
  except SyntaxError as ex:problems.append(rel+': '+str(ex))
 if p.suffix in {'.mp4','.png','.jpg','.gif','.wav'}:continue
 try:text=p.read_text(encoding='utf-8-sig')
 except UnicodeDecodeError:problems.append(rel+': unexpected binary');continue
 # Exclude this scanner's own pattern literals.
 if p.name=='verify_release.py':continue
 patterns={
  'private Windows path':r'[A-Za-z]:[\\/](?:Users|Workspace|CFD_with_AI_agents)[\\/]',
  'private Linux account path':r'/(?:home/nga2|mnt/[cd]/Users)/',
  'private IPv4':r'\b(?:166\.104\.\d+\.\d+|192\.168\.\d+\.\d+|10\.\d+\.\d+\.\d+)\b',
  'private key':r'-----BEGIN (?:OPENSSH|RSA|EC|DSA)?\s*PRIVATE KEY-----',
  'token-shaped literal':r'\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{30,}|sk-[A-Za-z0-9]{32,})\b',
 }
 for label,pattern in patterns.items():
  if re.search(pattern,text):problems.append(rel+': '+label)
m=json.loads((ROOT/'provenance/source_manifest.json').read_text())
for f in m['files']:
 if hashlib.sha256((ROOT/f['release']).read_bytes()).hexdigest()!=f['releaseSHA256']:problems.append(f['release']+': source hash mismatch')
c=json.loads((ROOT/'examples/catalog.json').read_text())
for name,entry in c.items():
 for mode in ['smoke','demo']:
  p=ROOT/entry['directory']/('input.'+mode);s=p.read_text()
  part=re.search(r'(?m)^Partition\s*:\s*([^\n]+)',s);v=[int(x)for x in part[1].split()]
  if v!=entry['partition']:problems.append(str(p.relative_to(ROOT))+': partition mismatch')
  if any(line.partition(':')[2].strip()!='none'for line in s.splitlines()if line.startswith('Restart prefix')):problems.append(str(p.relative_to(ROOT))+': non-cold start')
assert not problems,'\n'.join(problems)
print(json.dumps({'PASS':True,'files':len(files),'bytes':sum(p.stat().st_size for p in files),'sourceHashesChecked':len(m['files']),'entries':len(c),'note':'Pattern scan is scoped verification, not a guarantee against every possible secret or licensing issue.'}))
