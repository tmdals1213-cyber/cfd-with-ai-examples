#!/usr/bin/env python3
"""Small Linux controller: isolated builds, cold-start inputs and one owned MPI job."""
import argparse,csv,fcntl,hashlib,json,math,os,re,shlex,shutil,signal,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PIN='28a58d937af89dd506456048be1cad58f6d383a9'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def params(p):return {m[1].strip():m[2].strip() for m in re.finditer(r'(?m)^([^#\n:]+)\s*:\s*([^\n]*)',p.read_text())}
def catalog():return json.loads((ROOT/'examples/catalog.json').read_text())
def source_check(case):
 manifest=json.loads((ROOT/'provenance/source_manifest.json').read_text())
 for f in manifest['files']:
  if f['release'].startswith(catalog()[case]['directory']+'/'):
   assert digest(ROOT/f['release'])==f['releaseSHA256'],'Source changed; review and update provenance deliberately: '+f['release']
def upstream():
 p=Path(os.environ.get('NGA_HOME',str(ROOT/'vendor/NGA2'))).expanduser().resolve()
 assert p.is_dir(),'Set NGA_HOME or run bash scripts/fetch_nga2.sh'
 assert not any(c.isspace()for c in str(p)),'NGA make expects a path without spaces'
 actual=subprocess.check_output(['git','-C',str(p),'rev-parse','HEAD'],text=True).strip()
 assert actual==PIN,'Upstream commit differs from the pinned dependency'
 assert not subprocess.check_output(['git','-C',str(p),'status','--porcelain','--untracked-files=no'],text=True).strip(),'Upstream tracked files changed'
 return p
def case_build(c,jobs):
 up=upstream();source_check(c);src=ROOT/catalog()[c]['directory'];dst=ROOT/'build'/c;dst.mkdir(parents=True,exist_ok=True)
 shutil.copytree(src/'src',dst/'src',dirs_exist_ok=True);shutil.copyfile(src/'GNUmakefile',dst/'GNUmakefile')
 suffix='dp.gnu.opt.mpi';flags=f'-g -O2 -ftree-vectorize -ffree-line-length-none -fno-range-check -fno-second-underscore -Jbuild/o/{suffix} -Ibuild/o/{suffix} -fimplicit-none -fcheck=bounds -fbacktrace -ffpe-trap=zero,overflow'
 args=['make','-j',str(jobs),'NGA_HOME='+str(up),'DEBUG=FALSE','F90FLAGS='+flags,'FFLAGS='+flags.replace('-ffree-line-length-none','-ffixed-line-length-none')]
 started=time.time()
 with(dst/'build.log').open('w')as log:ret=subprocess.run(args,cwd=dst,stdout=log,stderr=subprocess.STDOUT)
 assert ret.returncode==0,'Build failed; inspect '+str(dst/'build.log')
 exe=dst/'nga.dp.gnu.opt.mpi.exe';assert exe.is_file()
 record={'case':c,'upstreamCommit':PIN,'source':{p.name:digest(p)for p in(src/'src').iterdir()if p.is_file()},'binarySHA256':digest(exe),'seconds':time.time()-started,'compiler':subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],'command':args}
 save(dst/'BUILD.json',record);print(json.dumps({'built':c,'seconds':record['seconds'],'binarySHA256':record['binarySHA256']}))
def audit(run,expected):
 log=(run/'run.log').read_text(errors='replace')
 assert not re.search(r'MPI_ABORT|ERROR STOP|SIGFPE|Floating.point exception|Fortran runtime error',log,re.I),'Solver signalled failure, even if return code was zero'
 p=run/'history.csv'
 if p.exists():
  with p.open()as f:rows=list(csv.DictReader(f))
  assert rows,'Empty custom history'
  values=[{k:float(v.replace('D','E'))for k,v in r.items()}for r in rows]
  assert all(math.isfinite(v)for r in values for v in r.values()),'Nonfinite history'
  assert abs(values[-1]['time_s']-expected)<max(1e-10,expected*1e-8),'Target time not reached'
  checks={'total_energy_error_J':1e-5,'div_minus_source_s_1':1e-4,'pressure_relative_error':1.01e-8}
  maxima={k:max(abs(r[k])for r in values if k in r)for k in checks if k in values[0]}
  assert all(maxima[k]<=checks[k]for k in maxima),'History residual/ledger limit exceeded'
  if 'mass_balance_relative' in values[0]:
   significant=[r for r in values if abs(r.get('phase_power_W',1))>1e-8]
   if significant:
    maxima['mass_balance_relative']=max(abs(r['mass_balance_relative'])for r in significant)
    assert maxima['mass_balance_relative']<=1e-4,'Mass residual exceeded (nonzero phase transfer)'
  return {'rows':len(values),'final':values[-1],'maxima':maxima,'scope':'Execution/history checks only; inspect saved phase/field positivity separately.'}
 # Upstream / early water examples use a native simulation monitor.
 monitor=run/'monitor/simulation';assert monitor.exists(),'Expected monitor/simulation missing'
 rows=[]
 for line in monitor.read_text().splitlines():
  try:row=[float(t.replace('D','E'))for t in line.split()]
  except ValueError:continue
  if row:assert all(math.isfinite(v)for v in row);rows.append(row)
 assert len(rows)>1,'No accepted monitor steps'
 assert abs(rows[-1][1]-expected)<max(1e-8,expected*1e-6),'Target time not reached'
 return {'monitorRows':len(rows),'lastMonitorRow':rows[-1],'scope':'Startup only. Inspect pressure/IBM/thermodynamic extrema before longer work.'}
def owned(exe,sid):
 out={}
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   raw=(p/'stat').read_text();fields=raw[raw.rfind(')')+2:].split()
   if int(fields[3])==sid and(p/'exe').resolve()==exe.resolve():out[int(p.name)]=int(fields[19])
  except(OSError,ValueError):pass
 return out
def close_owned(proc,exe):
 if proc.poll()is None:
  proc.terminate()
  try:proc.wait(timeout=25)
  except subprocess.TimeoutExpired:pass
 # MPI ranks may each have a separate process group. Match session/executable/start tick before any fallback signal.
 children=owned(exe,proc.pid)
 for pid,start in children.items():
  if owned(exe,proc.pid).get(pid)==start:
   try:os.kill(pid,signal.SIGTERM)
   except ProcessLookupError:pass
 deadline=time.monotonic()+5
 while owned(exe,proc.pid)and time.monotonic()<deadline:time.sleep(.2)
 for pid,start in children.items():
  if owned(exe,proc.pid).get(pid)==start:
   try:os.kill(pid,signal.SIGKILL)
   except ProcessLookupError:pass
 if proc.poll()is None:proc.kill();proc.wait()
 assert not owned(exe,proc.pid),'Owned ranks remain; inspect before another launch'
def run(c,choice):
 assert sys.platform.startswith('linux'),'Run inside Linux/WSL'
 source_check(c);b=ROOT/'build'/c;record=json.loads((b/'BUILD.json').read_text());exe=b/'nga.dp.gnu.opt.mpi.exe'
 assert digest(exe)==record['binarySHA256'],'Binary changed'
 assert {p.name:digest(p)for p in(ROOT/catalog()[c]['directory']/'src').iterdir()if p.is_file()}==record['source'],'Rebuild after source changes'
 original=ROOT/catalog()[c]['directory']/('input.'+choice);pr=params(original);partition=[int(x)for x in pr['Partition'].split()];ranks=math.prod(partition)
 assert ranks==catalog()[c]['ranks'],'Partition/rank count changed; review the case'
 assert pr.get('Restart prefix','none')=='none','Portable entries are cold-start only'
 parent=ROOT/'runs';parent.mkdir(exist_ok=True)
 with(parent/'.solver.lock').open('w')as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  stamp=time.strftime('%Y%m%d_%H%M%S');out=parent/(stamp+'_'+c);out.mkdir()
  shutil.copyfile(original,out/'input')
  data=ROOT/catalog()[c]['directory']/'data'
  if data.exists():shutil.copytree(data,out/'data')
  flags=['--bind-to','none']
  if os.environ.get('CFD_WSL1')=='1':
   (out/'private_shm').mkdir();(out/'private_tmp').mkdir()
   flags+=['--mca','plm','isolated','--mca','pml','ob1','--mca','btl','self,vader','--mca','btl_vader_single_copy_mechanism','none','--mca','btl_vader_backing_directory',str(out/'private_shm'),'--mca','orte_tmpdir_base',str(out/'private_tmp'),'--host',f'localhost:{ranks}']
  flags+=shlex.split(os.environ.get('CFD_MPI_EXTRA',''))
  args=['mpirun',*flags,'-np',str(ranks),str(exe),'-i','input'];env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
  d={'state':'RUNNING','case':c,'inputSHA256':digest(original),'binarySHA256':digest(exe),'partition':partition,'command':args};save(out/'EXECUTION.json',d)
  with(out/'run.log').open('w')as f:
   proc=subprocess.Popen(args,cwd=out,stdout=f,stderr=subprocess.STDOUT,env=env,start_new_session=True);d['launcherPID']=proc.pid;save(out/'EXECUTION.json',d)
   print(str(out),flush=True)
   def interrupt_owned(signum,frame):raise KeyboardInterrupt('Owned solver stop requested')
   old_term=signal.signal(signal.SIGTERM,interrupt_owned)
   try:
    ret=proc.wait();assert ret==0,'Nonzero solver return code';d['audit']=audit(out,float(pr['Max time']));d['state']='COMPLETED_EXECUTION_CHECKS'
   except BaseException as ex:
    d['state']='FAILED_OR_INTERRUPTED_PRESERVED';d['error']=str(ex);raise
   finally:
    signal.signal(signal.SIGTERM,old_term)
    close_owned(proc,exe);d['allOwnedRanksClosed']=True;d['exitCode']=proc.returncode;save(out/'EXECUTION.json',d)
  print(json.dumps(d,indent=2))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['list','build','run']);p.add_argument('case',nargs='?',choices=sorted(catalog()));p.add_argument('--jobs',type=int,default=8);p.add_argument('--input',choices=['smoke','demo','rest'],default='smoke');a=p.parse_args()
 if a.action=='list':
  for k,v in catalog().items():print(k,'|',v['title'],'|',v['ranks'],'MPI')
 elif not a.case:p.error('Select a case')
 elif a.action=='build':assert a.jobs>0;case_build(a.case,a.jobs)
 else:run(a.case,a.input)
if __name__=='__main__':main()
