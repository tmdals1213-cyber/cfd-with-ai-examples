"""Continuous actual inverted8 simulation to a full gait cycle; no packet-confirmation stop."""
import ast,hashlib,json,mmap,os,secrets,shutil,sys,time
from fractions import Fraction as F
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
sys.dont_write_bytecode=True
H=Path(__file__).resolve().parent
import local_owned_runtime as W

def digest(b):return hashlib.sha256(b).hexdigest()
def pin(q):
 b=q.read_bytes();assert len(b)<=134217728;return dict(path=str(q),bytes=len(b),SHA=digest(b))
def write(q,v):
 with q.open('x',encoding='utf8')as f:json.dump(v,f,indent=2)
 return pin(q)
def readmeta(q):
 b=q.read_bytes();assert len(b)<=1048576;return json.loads(b)
ad=readmeta(H/'LIVE_CONTINUATION_ADMISSION.json');assert ad['admitted']and ad['masterSHA']==digest(Path(__file__).read_bytes())
for p in ad['pins']:
 b=Path(p['path']).read_bytes();assert len(b)==p['bytes']and digest(b)==p['SHA']
ast.parse((H/'continuous_worker.py').read_bytes());cfg=readmeta(H/'BASE_LOCAL_INPUT.json');control=readmeta(H/'CONTINUOUS_MEMORY_RUNTIME_CONTROL_PASS.json');assert control['PASS']
O=Path(cfg['dataRoot']).resolve();assert str(O).lower()==str(Path(ad['dataRoot']).resolve()).lower()and O.is_dir()
assert shutil.disk_usage(O).free>=16*1073741824
write(H/'LIVE_CONTINUATION_STARTED.json',dict(pid=os.getpid(),startedEpoch=time.time(),startResult=cfg['initialAcceptedResult'],startTimeExact=cfg['initialTimeExact'],targetEndExact=cfg['targetEndExact'],campaignTotalDurationExact='1281/512',additionalTargetDurationExact=str(F(cfg['targetEndExact'])-F(cfg['initialTimeExact'])),maximumAdditionalSteps=10000,noPacketConfirmationStop=True,dataRoot=str(O)))
k=W.job_api();assert k.SetPriorityClass(k.GetCurrentProcess(),16384)==1;t0=time.monotonic();accepted=readmeta(Path(cfg['initialAcceptedResult']['path']));acceptedPin=cfg['initialAcceptedResult'];seen=list(cfg['priorC1AlreadyUsed']);history=[];currentStep=0;usedWrite=0;maps={};S=None
staticPins=[ad['worker'],ad['common'],ad['base']]
def role(mode,S,extra,mask=21845,threads=8,suffix=''):
 assert time.monotonic()-t0<172800
 f=dict(mode=mode,pins=staticPins+extra.pop('pins',[]),**extra);name=mode.upper()+suffix;fp=write(S/('FRAME_'+name+'.json'),f);W.THREADS=threads
 r=W.owned([sys.executable,'-B',str(H/'continuous_worker.py'),mode,fp['path'],fp['SHA']],mask=mask,wall=30)
 write(S/('JOB_'+name+'.json'),r)
 assert r.get('returncode')==0 and not r.get('timedOut')and not r.get('error')and r.get('allOwnedHandlesClosed'),('Owned stage failure',mode,suffix,r)
 return r

def preserve_failed_maps():
 if S is None:return
 for name in('state','trace','stats','output','rhs'):
  if name in maps:
   q=S/('FAILURE_BUFFER_'+name+'.bin')
   if not q.exists():
    b=maps[name][:];assert len(b)<=134217728
    with q.open('xb')as f:f.write(b)
 if 'prepared'in maps:
  mm=maps['prepared'];size=len(mm)
  for n,start in enumerate(range(0,size,134217728)):
   q=S/('FAILURE_PREPARED_PART_%02d.bin'%n)
   if not q.exists():
    with q.open('xb')as f:f.write(mm[start:min(size,start+134217728)])

def close_maps():
 for mm in reversed(list(maps.values())):mm.close()
 maps.clear()

try:
 for step in range(1,10001):
  if F(accepted['timeExact'])>=F(cfg['targetEndExact']):break
  assert usedWrite<128*1073741824 and shutil.disk_usage(O).free>=16*1073741824 and time.monotonic()-t0<172800
  if(H/'STOP_AT_ACCEPTED_BOUNDARY.json').exists():raise RuntimeError('User stop at accepted-state boundary')
  currentStep=step;S=O/('STEP_%06d'%step);S.mkdir();motionOrdinal=cfg['initialMotionOrdinal']+step
  binding=write(S/'STEP_SOURCE_BINDING.json',dict(stepNumber=motionOrdinal,acceptedFluidStep=cfg['initialAcceptedStepCount']+step,initialMomentum=accepted['newP'],initialTimeExact=accepted['timeExact'],previousResult=acceptedPin,previousC1Bits=seen))
  role('cfl',S,dict(binding=binding,pins=[binding]),threads=1);tp=pin(S/'NEXT_TIME_PLAN.json');timeplan=readmeta(Path(tp['path']))
  role('motion',S,dict(timePlan=tp,pins=[tp,cfg['motionCache']]));mo=pin(S/'MOTION_PASS.json');motion=readmeta(Path(mo['path']));seen.extend(motion['newC1Bits'])
  def geometry(a):return role('geometry',S,dict(axis=a,timePlan=tp,motion=mo,pins=[tp,mo]),mask=[4,16,64][a],threads=1,suffix='_%d'%a)
  with ThreadPoolExecutor(max_workers=3)as pool:list(pool.map(geometry,range(3)))
  role('sgs_mask',S,dict(timePlan=tp,motion=mo,pins=[tp,mo]),threads=1)
  gp=[pin(S/('EXPOSED47_AXIS%d_PASS.json'%a))for a in range(3)];role('training',S,dict(timePlan=tp,motion=mo,pins=[tp,mo]+gp));tr=pin(S/'TRAINING_AND_PC_PASS.json');tm=readmeta(Path(tr['path']));D=tm['D'];N=tm['N'];K=tm['K'];L=N+K
  nonce=secrets.token_hex(16);prefix='Local\\LICA_I8LIVE_'+nonce;counts=dict(prepared=10*D+N,rhs=L,state=3*L+32,trace=385*9,output=3*D+2*L,stats=12,predictorStats=6);desc={key:dict(name=prefix+'_'+key,bytes=count*8)for key,count in counts.items()};desc['guard']=dict(name=prefix+'_guard',bytes=64);desc['guardHex']=('I8LIVE:'+nonce).encode().ljust(64,b'\0').hex();maps={}
  try:
   for key in list(counts)+['guard']:maps[key]=mmap.mmap(-1,desc[key]['bytes'],tagname=desc[key]['name'],access=mmap.ACCESS_WRITE)
   maps['guard'][:64]=bytes.fromhex(desc['guardHex']);descPin=write(S/'LIVE_SHARED_BUFFER_DESCRIPTOR.json',desc)
   role('prepare',S,dict(timePlan=tp,motion=mo,maps=desc,mapsPin=descPin,pins=[tp,mo,tr,descPin]));prep=pin(S/'ACTUAL_PREPARED_PASS.json');pm=readmeta(Path(prep['path']));assert pm['L']==L
   previous=None;slices=[];beg=time.monotonic()
   for chunk in range(6):
    rp=role('pressure',S,dict(chunk=chunk,previous=previous,maps=desc,mapsPin=descPin,timePlan=tp,motion=mo,pins=[tp,mo,prep,tr,descPin]+([previous]if previous else[])),mask=21845,threads=8,suffix='_%03d'%chunk)
    previous=pin(S/('SLICE_%03d.json'%chunk));sl=readmeta(Path(previous['path']));assert sl['PASS']and sl['transactionId']==pm['transactionId']and sl['rawBSHA']==pm['rawBSHA'];slices.append(dict(chunk=chunk,iteration=sl['iteration'],status=sl['status'],sliceSHA=previous['SHA'],pid=rp['pid'],allOwnedHandlesClosed=True))
    print(json.dumps(dict(stage='LIVE_INVERTED8_PRESSURE',acceptedStep=cfg['initialAcceptedStepCount']+step,newStep=step,iteration=sl['iteration'],status=sl['status'],full=sl['trueFullRaw'],tau=sl['tauFull'],CFL=pm['CFL'],dt_s=timeplan['dt_s'],timeExact=timeplan['timeExact'],seconds=time.monotonic()-beg)),flush=True)
    if sl['status']==0:
     write(S/'PRESSURE_TERMINAL_COMPLETE.json',dict(PASS=True,lastSlice=previous,iterations=sl['iteration'],slices=slices,wall_s=time.monotonic()-beg,transactionId=pm['transactionId'],allOwnedHandlesClosed=True));break
    assert sl['status']==101
   else:raise RuntimeError('Global384 exhausted without restart')
   term=pin(S/'PRESSURE_TERMINAL_COMPLETE.json');role('audit',S,dict(timePlan=tp,motion=mo,maps=desc,mapsPin=descPin,pins=[tp,mo,prep,tr,term,descPin]));acceptedPin=pin(S/'INVERTED8_EXPLORATORY_ACTUAL_RESULT.json');accepted=readmeta(Path(acceptedPin['path']));assert accepted['PASS']
  except BaseException:
   preserve_failed_maps();raise
  finally:close_maps()
  hist=dict(newStep=step,acceptedFluidStep=cfg['initialAcceptedStepCount']+step,timeExact=accepted['timeExact'],dtExact=accepted['dtExact'],CFL=pm['CFL'],diffusion=pm['diffusion'],iterations=accepted['iterations'],normalizedContinuity=accepted['normalizedContinuity'],result=acceptedPin,field=accepted['newP'])
  history.append(hist);write(S/'STEP_COMMIT.json',hist);usedWrite+=sum(q.stat().st_size for q in S.iterdir()if q.is_file());print(json.dumps(dict(stage='LIVE_ACTUAL_TIME_STEP_COMMITTED',**hist)),flush=True)
 targetReached=F(accepted['timeExact'])>=F(cfg['targetEndExact']);result=dict(PASS=targetReached,status='TARGET_DURATION_COMPLETED'if targetReached else'STEP_LIMIT_UNFINISHED',newCFDSteps=len(history),history=history,finalResult=acceptedPin,campaignTotalDurationExact='1281/512',additionalTargetDurationExact=str(F(cfg['targetEndExact'])-F(cfg['initialTimeExact'])),initialTimeExact=cfg['initialTimeExact'],finalTimeExact=accepted['timeExact'],allOwnedActorsClosed=True,allSharedBuffersClosed=True,newMotionQueries=3*len(history),physicalTotalDragQualified=False,writeRetained_B=usedWrite,wall_s=time.monotonic()-t0,CFDCurrentlyRunning=False)
 write(H/'LIVE_CONTINUATION_FINAL_CLOSURE.json',result);print(json.dumps(dict(stage=result['status'],newSteps=len(history),finalTimeExact=result['finalTimeExact'])),flush=True)
except BaseException as e:
 if maps:
  try:preserve_failed_maps()
  finally:close_maps()
 write(H/'LIVE_CONTINUATION_CRITICAL_FAILURE.json',dict(PASS=False,error=repr(e),newStep=currentStep,successfulHistory=history,finalAcceptedResult=acceptedPin,failedRHSReplayForbidden=True,allOwnedHandleCleanupBySource=True,wall_s=time.monotonic()-t0));raise
