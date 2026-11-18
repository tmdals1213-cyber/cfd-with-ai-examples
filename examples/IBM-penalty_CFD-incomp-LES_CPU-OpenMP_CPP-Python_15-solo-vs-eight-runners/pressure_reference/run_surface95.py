from pathlib import Path
from fractions import Fraction as F
from concurrent.futures import ThreadPoolExecutor
import sys,json,time,hashlib,importlib.util,threading
R=Path(__file__).resolve().parent.parent;WS=R.parents[1];OLD=WS/'work/pacemaker_penalty_les_20261112';REC=WS/'work/pacemaker_penalty_les_recovery_20261113';sys.dont_write_bytecode=True
def pin(p):
 b=p.read_bytes();return dict(path=str(p.resolve()),bytes=len(b),SHA=hashlib.sha256(b).hexdigest())
def save(p,v):
 with p.open('x',encoding='utf8')as f:json.dump(v,f,indent=2)
spec=importlib.util.spec_from_file_location('surface95_owned',REC/'group8/local_owned_runtime.py');W=importlib.util.module_from_spec(spec);spec.loader.exec_module(W);W.H=R;W.THREADS=1
tasks=[];closures=[]
for case in('solo','group8'):
 closure=REC/case/'LIVE_CONTINUATION_FINAL_CLOSURE.json';end=json.loads(closure.read_text());assert end['PASS']and F(end['finalTimeExact'])>=F(5,2);closures.append(pin(closure));entries=[]
 for root in(OLD,REC):
  for q in sorted((root/'data'/case).glob('STEP_*/STEP_COMMIT.json')):
   x=json.loads(q.read_text());entries.append(dict(time_s=float(F(x['timeExact'])),result=x['result'],acceptedFluidStep=x['acceptedFluidStep']))
 entries.sort(key=lambda x:x['time_s']);targets=[.5+j*.025 for j in range(81)];chosen={min(range(len(entries)),key=lambda i:abs(entries[i]['time_s']-t))for t in targets}
 for t in(1.25,2.5):
  chosen.add(max(i for i,e in enumerate(entries)if e['time_s']<=t));chosen.add(min(i for i,e in enumerate(entries)if e['time_s']>=t))
 for ordinal,i in enumerate(sorted(chosen)):
  O=R/'data'/('%s_%03d'%(case,ordinal));O.mkdir(exist_ok=False);req=dict(case=case,entry=entries[i],config=pin(OLD/case/'BASE_LOCAL_INPUT.json'),worker=pin(R/'source/surface95_frame.py'),kernel=pin(R/'source/pressure_kernel.py'),geometryWorker=pin(R/'source/exposed_surface.py'))
  save(O/'REQUEST.json',req);tasks.append(O)
save(R/'SURFACE95_RUNTIME_PLAN.json',dict(admitted=True,reason='Saved47 quadrature needed excessive moment correction for a pose; use genuinely finer union quadrature, not looser tolerances',
 originalFailedHistoryPreserved=pin(WS/'work/runner_pressure_load_history_20261113/PRESSURE_HISTORY_ATTENTION.json'),completedCFDClosures=closures,
 inputOutputComposition=dict(frameRoles=len(tasks),filesPerFrame=5,controlAndSourceFileUpper=32,physicalFilesUpper=1024),
 newEnvelope=dict(readUpper_B=4*1073741824,writeRetainUpper_B=1073741824,physicalFilesUpper=1024),
 noOldBudgetReuse=True,geometryQuadrature=95,maxRelativeAreaCorrectionUnchanged=.25,workers=16,CPUPerWorker=1,workerWall_s=30,memoryCap=None,
 pressureLoadModelEstimate=True,viscousForceIncluded=False,physicalPressureExteriorLimitQualified=False,physicalTotalDragQualified=False,newCFDSteps=0))
lock=threading.Lock();nextcore=0;local=threading.local()
def init():
 global nextcore
 with lock:local.core=nextcore;nextcore+=1
def run(O):
 rec=W.owned([sys.executable,'-B',str(R/'source/surface95_frame.py'),str(O/'REQUEST.json')],mask=1<<(2*local.core),stdout_path=O/'STDOUT.log',wall=30);save(O/'JOB.json',rec)
 if not(rec.get('returncode')==0 and rec.get('allOwnedHandlesClosed')and not rec.get('timedOut')):return dict(PASS=False,folder=str(O),job=rec)
 out=json.loads((O/'PRESSURE_LOAD.json').read_text());assert out['PASS'];return out
with ThreadPoolExecutor(max_workers=16,initializer=init)as pool:
 futures=[pool.submit(run,O)for O in tasks];outputs=[f.result()for f in futures]
good=all(x['PASS']for x in outputs);save(R/'SURFACE95_CLOSURE.json',dict(PASS=good,frames=outputs,newCFDSteps=0,physicalTotalDragQualified=False))
assert good,[(x.get('folder'),x.get('job',{}).get('returncode'))for x in outputs if not x['PASS']]
print(json.dumps(dict(PASS=True,frames=len(outputs),cases=2)),flush=True)
