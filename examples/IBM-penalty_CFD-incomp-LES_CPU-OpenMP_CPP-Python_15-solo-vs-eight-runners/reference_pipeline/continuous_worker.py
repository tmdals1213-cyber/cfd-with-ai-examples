"""Source-bound adaptive-CFL continuation of the actual inverted8 field."""
import sys,io,json,time,ctypes as C,struct,math,traceback
from pathlib import Path
from fractions import Fraction as F
sys.dont_write_bytecode=True
import runtime_common_local as W
np,p,t0=W.start()
from divfree_mac_reconstruction import rows,coordinates
from tensor_pressure_pc import make,save
from surface_block_preconditioner import make_blocks,save_blocks,make_action
framePath=Path(sys.argv[2]);body=framePath.read_bytes();assert len(body)<=1048576 and W.sha(body)==sys.argv[3];frame=json.loads(body);S=framePath.parent.resolve();assert S.parent==Path(p['dataRoot']).resolve()and S.name.startswith('STEP_');W.H=S
for pin in frame['pins']:W.frozen(pin)
if frame.get('mapsPin'):assert frame['maps']==json.loads(W.frozen(frame['mapsPin']))
mode=sys.argv[1];DT=None
if frame.get('timePlan'):
 tp=json.loads(W.frozen(frame['timePlan']));assert tp['PASS'];DT=float(F(tp['dtExact']));p['dtExact']=tp['dtExact'];p['timeExact']=tp['timeExact'];p['initialMomentum']=tp['initialMomentum']
if frame.get('motion'):
 mo=json.loads(W.frozen(frame['motion']));assert mo['PASS'];p['pose']=mo['pose'];p['poseEpoch']=mo['poseEpoch']
def cfl_plan():
 b=json.loads(W.frozen(frame['binding']));previous=json.loads(W.frozen(b['previousResult']));assert previous['PASS']and previous['newP']['SHA']==b['initialMomentum']['SHA']and previous['timeExact']==b['initialTimeExact']
 faces,n,N,h,shapes,sizes,off,D,mass,areas,md=W.grid(np,p['grid']);old=np.load(io.BytesIO(W.frozen(b['initialMomentum'])),allow_pickle=False);assert old.shape==(D,)and W.sha(old.tobytes())==b['initialMomentum']['momentumPayloadSHA']
 lib=W.dll(p['fineDLL']);fn=lib.pm_flow_rates_les;fn.argtypes=[W.DP]*3+[C.c_double]*3+[C.c_int,W.DP,W.DP];fn.restype=C.c_int;rates=np.zeros(11);assert fn(*[W.ptr(q)for q in faces],1.225,1.5e-5,.035,1,W.ptr(old),W.ptr(rates))==0
 assert np.isfinite(rates).all()and rates[0]>0 and max(rates[5:7])*rates[4]/(2813/480)<=1e-4
 den=1048576;initial=F(b['initialTimeExact']);remaining=F(p['targetEndExact'])-initial;m=min(int(np.floor(.45/(rates[0]+rates[1])*den)),int(remaining*den),int(p['bodyDtUpper_s']*den));assert m>0;eps=F(3,268435456);forbidden=set(json.loads(W.frozen(p['priorC1Registry']))['bits'])|set(b['previousC1Bits']);bits=None
 for trial in range(4096):
  assert m>0;dt=F(m,den);target=initial+dt;candidate=[struct.pack('<d',float(target+d)).hex()for d in(F(0),eps,-eps)]
  if len(set(candidate))==3 and not any(q in forbidden for q in candidate):bits=candidate;break
  m-=1
 assert bits is not None and float(dt)*(rates[0]+rates[1])<=.5
 r=dict(PASS=True,initialMomentum=b['initialMomentum'],previousResult=b['previousResult'],initialTimeExact=str(initial),timeExact=str(target),dtExact=str(dt),dt_s=float(dt),maxAdvectiveRate_s_1=float(rates[0]),maxDiffusiveRate_s_1=float(rates[1]),maxCombinedLocalRate_s_1=float(rates[2]),CFL=float(dt*rates[0]),diffusion=float(dt*rates[1]),combinedConservative=float(dt*(rates[0]+rates[1])),CFLTarget=.45,CFLUpper=.5,previousMaximumFaceVelocity_m_s=float(rates[3]),newC1BitTriplet=bits,poseEpoch=600000+b['stepNumber'],stepNumber=b['stepNumber'],read_B=W.READ)
 W.save('NEXT_TIME_PLAN.json',r);print(json.dumps(dict(mode='CFL_PLAN',step=b['stepNumber'],CFL=r['CFL'],dt_s=r['dt_s'],timeExact=r['timeExact'])),flush=True)

def motion():
 tp=json.loads(W.frozen(frame['timePlan']));assert tp['PASS'];raw=W.frozen(p['motionCache']);top=W.frozen(p['topology']);from motion_model import Motion
 model=Motion(io.BytesIO(raw),192.,mode='periodic_c1');tri=np.frombuffer(top,'<i4',49112*3,84).reshape(-1,3);canon=np.frombuffer(top,'<i4',28374,84+49112*24);_,rep=np.unique(canon,return_index=True);assert np.array_equal(model.faces,tri)
 target=F(tp['timeExact']);eps=F(3,268435456)
 seen=[]
 W.save('C1_INTENT.json',dict(beforeC1=True,once=True,timeExact=str(target),bits=tp['newC1BitTriplet'],sourceSHA=p['motionCache']['SHA']))
 def sample(t):
  bit=struct.pack('<d',float(t)).hex();assert bit in tp['newC1BitTriplet']and bit not in seen;seen.append(bit);return model.sample(float(t%F(5,8)))
 xyz,vel=sample(target);xp,_=sample(target+eps);xm,_=sample(target-eps);ids=np.array([107,1201,6007,13421,18013,24011,27971]);fd=float(abs((xp[ids]-xm[ids])/(2*float(eps))-vel[ids]).max());assert fd<=2e-5 and np.array_equal(xyz,xyz[rep[canon]])and np.array_equal(vel,vel[rep[canon]])
 lo=xyz.min(0);hi=xyz.max(0);off=np.array(p['offsets_m']);bl=lo+off;bh=hi+off;domain=np.array([[-8,6],[0,4],[-6,6]]);assert np.all(bl>domain[:,0])and np.all(bh<domain[:,1]);gaps=[]
 for i in range(len(p['offsets_m'])):
  for j in range(i):
   gap=float(np.maximum(bl[i]-bh[j],bl[j]-bh[i]).max());assert gap>0;gaps.append(gap)
 body=struct.pack('<8sidQ',b'RDXPOS01',28374,float(target),tp['poseEpoch'])+np.ascontiguousarray(xyz,dtype='<f8').tobytes()+np.ascontiguousarray(vel,dtype='<f8').tobytes()
 with(W.H/'pose.bin').open('xb')as f:f.write(body)
 r=dict(PASS=True,pose=W.pinfile('pose.bin'),poseEpoch=tp['poseEpoch'],timeExact=str(target),dtExact=tp['dtExact'],FDmax=fd,canonicalWeldingByteExact=True,minimumAABBSeparatingGap_m=min(gaps)if gaps else None,allBodiesInsideSameDomain=True,C1Queries=3,newC1Bits=seen,sourceMotionSHA=p['motionCache']['SHA'],read_B=W.READ)
 W.save('MOTION_PASS.json',r);print(json.dumps(dict(mode='MOTION',step=tp['stepNumber'],FDmax=fd,timeExact=str(target))),flush=True)

def geometry():
 axis=frame['axis'];assert 0<=axis<3
 pose=W.frozen(p['pose']);top=W.frozen(p['topology']);magic,nv,t,epoch=struct.unpack_from('<8sidQ',pose);assert magic==b'RDXPOS01'and nv==28374 and t==float(F(p['timeExact']))and epoch==p['poseEpoch']
 xyz=np.frombuffer(pose,dtype='<f8',count=nv*3,offset=28).copy().reshape(-1,3);vel=np.frombuffer(pose,dtype='<f8',count=nv*3,offset=28+nv*24).copy().reshape(-1,3);lo=xyz.min(0);hi=xyz.max(0);tri=np.frombuffer(top,'<i4',49112*3,84).reshape(-1,3);component=np.frombuffer(top,'<i4',49112,84+49112*12);context=top[20:84]
 class Event(C.Structure):_fields_=[('info',C.c_int*4),('data',C.c_double*13)]
 assert C.sizeof(Event)==120
 lib=W.dll(p['geometryDLL']);lib.rnu1_rd_load.argtypes=[C.c_char_p,C.c_char_p,C.POINTER(C.c_void_p)];lib.rnu1_rd_load.restype=C.c_int;lib.rnu1_rd_update.argtypes=[C.c_void_p,C.c_int,W.DP,W.DP,C.c_double,C.c_uint64,C.c_char_p];lib.rnu1_rd_update.restype=C.c_int;lib.rnu1_rd_segment.argtypes=[C.c_void_p,W.DP,W.DP,C.c_double,C.c_uint64,C.c_int,C.c_int,W.IP,W.IP,W.IP,W.IP,C.POINTER(Event)];lib.rnu1_rd_segment.restype=C.c_int;lib.rnu1_rd_release.argtypes=[C.c_void_p]
 handle=C.c_void_p();assert lib.rnu1_rd_load(p['topology']['path'].encode('utf8'),context,C.byref(handle))==0 and handle.value
 assert lib.rnu1_rd_update(handle,nv,W.ptr(xyz),W.ptr(vel),t,1,context)==0
 NQ=47;other=[a for a in range(3)if a!=axis];planeArea=(hi[other[0]]-lo[other[0]])*(hi[other[1]]-lo[other[1]])/(NQ*NQ);records=[];volume=0.;eventsTotal=0;maxNormal=maxBary=maxWall=0.;result=None
 try:
  for i in range(NQ):
   for j in range(NQ):
    a=np.zeros(3);b=np.zeros(3)
    for q,ax in enumerate(other):a[ax]=b[ax]=lo[ax]+((i if q==0 else j)+.5)*(hi[ax]-lo[ax])/NQ
    a[axis]=lo[axis]-.125;b[axis]=hi[axis]+.125;raw=C.c_int();num=C.c_int();si=C.c_int();ei=C.c_int();ev=(Event*256)();code=lib.rnu1_rd_segment(handle,W.ptr(a),W.ptr(b),t,1,1,256,C.byref(raw),C.byref(num),C.byref(si),C.byref(ei),ev)
    assert code==0,('Original provider refusal',axis,i,j,code)
    assert si.value==ei.value==0 and 0<=num.value<=raw.value<=256 and num.value%2==0
    previous=0.;entry=None
    for k in range(num.value):
     info=list(ev[k].info);data=list(ev[k].data);sign=-1 if k%2==0 else 1;assert info[3]==sign and component[info[0]]==info[1]and previous<data[0]<1 and data[7+axis]*sign>0;previous=data[0]
     normal=np.array(data[7:10]);bary=np.array(data[10:13]);ids=tri[info[0]];point=bary@xyz[ids];wall=bary@vel[ids]
     ne=abs(float(normal@normal)-1);be=abs(float(bary.sum())-1);we=float(np.max(np.abs(wall-data[4:7])));pe=float(np.max(np.abs(point-data[1:4])))
     assert np.isfinite(data).all()and ne<=1e-10 and be<=1e-12 and we<=2e-8 and pe<=2e-10
     maxNormal=max(maxNormal,ne);maxBary=max(maxBary,be);maxWall=max(maxWall,we);eventsTotal+=1
     if sign<0:entry=data[1+axis]
     else:assert entry is not None and data[1+axis]>entry;volume+=(data[1+axis]-entry)*planeArea;entry=None
     if int(np.argmax(abs(normal)))==axis:records.append((data[1:4],data[4:7],data[7:10],planeArea/abs(normal[axis]),info[1],axis,info[3]))
    assert entry is None
  result=dict(PASS=True,axis=axis,N=NQ,queries=NQ*NQ,events=eventsTotal,retainedDominantNormalPoints=len(records),volumeMidpointEstimate_m3=float(volume),maxNormalError=maxNormal,maxBaryError=maxBary,maxWallError=maxWall,sourcePoseSHA=p['pose']['SHA'],newC1Calls=0,newMotionCalls=0,originalProviderNumericsUnmodified=True)
 finally:lib.rnu1_rd_release(handle)
 assert result is not None
 name='EXPOSED47_AXIS%d.npz'%axis
 with(W.H/name).open('xb')as f:np.savez(f,points=np.array([r[0]for r in records]),wall=np.array([r[1]for r in records]),normal=np.array([r[2]for r in records]),area=np.array([r[3]for r in records]),component=np.array([r[4]for r in records],dtype='i4'),axis=np.array([r[5]for r in records],dtype='i4'),sign=np.array([r[6]for r in records],dtype='i4'))
 result['asset']=W.pinfile(name);result['wall_s']=time.monotonic()-t0;result['read_B']=W.READ;W.save('EXPOSED47_AXIS%d_PASS.json'%axis,result);print(json.dumps(result),flush=True)
 
def sgs_mask():
 pose=W.frozen(p['pose']);top=W.frozen(p['topology']);magic,nv,t,epoch=struct.unpack_from('<8sidQ',pose)
 assert magic==b'RDXPOS01'and epoch==p['poseEpoch']and t==float(F(p['timeExact']))
 xyz=np.frombuffer(pose,dtype='<f8',count=nv*3,offset=28).copy().reshape(-1,3)
 vel=np.frombuffer(pose,dtype='<f8',count=nv*3,offset=28+nv*24).copy().reshape(-1,3)
 faces,n,N,h,shapes,sizes,off,D,mass,areas,md=W.grid(np,p['grid']);cent=[.5*(f[:-1]+f[1:])for f in faces]
 class Event(C.Structure):_fields_=[('info',C.c_int*4),('data',C.c_double*13)]
 lib=W.dll(p['geometryDLL']);context=top[20:84]
 lib.rnu1_rd_load.argtypes=[C.c_char_p,C.c_char_p,C.POINTER(C.c_void_p)]
 lib.rnu1_rd_update.argtypes=[C.c_void_p,C.c_int,W.DP,W.DP,C.c_double,C.c_uint64,C.c_char_p]
 lib.rnu1_rd_segment.argtypes=[C.c_void_p,W.DP,W.DP,C.c_double,C.c_uint64,C.c_int,C.c_int,W.IP,W.IP,W.IP,W.IP,C.POINTER(Event)]
 lib.rnu1_rd_release.argtypes=[C.c_void_p]
 handle=C.c_void_p();assert lib.rnu1_rd_load(p['topology']['path'].encode('utf8'),context,C.byref(handle))==0
 assert lib.rnu1_rd_update(handle,nv,W.ptr(xyz),W.ptr(vel),t,1,context)==0
 fluid=np.ones(n,dtype='u1');lo=xyz.min(0);hi=xyz.max(0);queries=0
 try:
  for offset in np.array(p['offsets_m']):
   js=np.flatnonzero((cent[1]-offset[1]>lo[1])&(cent[1]-offset[1]<hi[1]))
   ks=np.flatnonzero((cent[2]-offset[2]>lo[2])&(cent[2]-offset[2]<hi[2]))
   for j in js:
    for k in ks:
     a=np.ascontiguousarray([lo[0]-.125,cent[1][j]-offset[1],cent[2][k]-offset[2]])
     b=a.copy();b[0]=hi[0]+.125;raw=C.c_int();num=C.c_int();si=C.c_int();ei=C.c_int();ev=(Event*256)()
     code=lib.rnu1_rd_segment(handle,W.ptr(a),W.ptr(b),t,1,1,256,C.byref(raw),C.byref(num),C.byref(si),C.byref(ei),ev)
     assert code==0 and si.value==ei.value==0 and num.value%2==0,('SGS exact union row refusal',j,k,code)
     queries+=1
     for q in range(0,num.value,2):
      assert ev[q].info[3]==-1 and ev[q+1].info[3]==1
      x0=ev[q].data[1]+offset[0];x1=ev[q+1].data[1]+offset[0]
      inside=(cent[0]>=x0-1e-10)&(cent[0]<=x1+1e-10);fluid[inside,j,k]=0
 finally:lib.rnu1_rd_release(handle)
 asset=W.np_save('SGS_FLUID_MASK.npy',fluid.ravel(),np)
 W.save('SGS_FLUID_MASK_PASS.json',dict(PASS=True,asset=asset,fluidCells=int(fluid.sum()),solidCells=int((fluid==0).sum()),
   originalExactUnionQueries=queries,poseSHA=p['pose']['SHA'],newMotionQueries=0))

def training():
 faces,n,N,h,shapes,sizes,off,D,mass,areas,md=W.grid(np,p['grid']);pieces=[]
 for a in range(3):
  meta=json.loads((W.H/('EXPOSED47_AXIS%d_PASS.json'%a)).read_text('utf8'));assert meta['PASS']and meta['N']==47 and meta['sourcePoseSHA']==p['pose']['SHA']
  with np.load(io.BytesIO(W.frozen(meta['asset'])),allow_pickle=False)as z:pieces.append({k:z[k].copy()for k in z.files})
 raw={k:np.concatenate([d[k]for d in pieces])for k in pieces[0]};assert np.all(raw['area']>0);selected=[];counts=[]
 for person,offset in enumerate(np.array(p['offsets_m'])):
  groups={}
  for j,point0 in enumerate(raw['points']):
   point=point0+offset;cell=tuple(int(np.searchsorted(faces[a],point[a],side='right')-1)for a in range(3));assert all(0<=cell[a]<n[a] for a in range(3));key=cell+(int(raw['component'][j]),int(raw['axis'][j]),int(raw['sign'][j]));groups.setdefault(key,[]).append(j)
  start=len(selected)
  for key,jj in sorted(groups.items()):
   aa=raw['area'][jj];area=float(aa.sum());centre=np.sum((raw['points'][jj]+offset)*aa[:,None],axis=0)/area;j=jj[int(np.argmin(np.sum((raw['points'][jj]+offset-centre)**2,axis=1)))];selected.append((raw['points'][j]+offset,raw['wall'][j],raw['normal'][j],area,person,raw['component'][j]))
  counts.append(len(selected)-start)
 points=np.ascontiguousarray([q[0]for q in selected]);wall=np.ascontiguousarray(np.array([q[1]for q in selected]).ravel());normal=np.array([q[2]for q in selected]);area=np.repeat([q[3]for q in selected],3);person=np.array([q[4]for q in selected],dtype='i4');component=np.array([q[5]for q in selected],dtype='i4');K=3*len(points);assert K<=49152
 ops=[rows(faces,points,a)for a in range(3)];ids=np.ascontiguousarray(np.stack([v[0]for v in ops],axis=1).reshape(-1,24),dtype='i8');w=np.ascontiguousarray(np.stack([v[1]for v in ops],axis=1).reshape(-1,24));comp,pos=coordinates(faces,ids);constantError=0.
 for a in range(3):
  mask=comp==a;constantError=max(constantError,float(abs(np.sum(w*mask,axis=1)-(np.arange(K)%3==a)).max()));q=ids[mask]-off[a];shape=shapes[a];j=q//(shape[1]*shape[2])if a==0 else q//shape[2]%shape[1]if a==1 else q%shape[2];assert np.all((j>=1)&(j<n[a]))
 U=np.array([.23,-.11,.17]);G=np.array([[.12,-.3,.2],[.4,-.07,-.1],[.23,.14,-.05]]);known=U[comp]+np.sum(G[comp]*pos,axis=2);got=np.sum(w*known,axis=1).reshape(-1,3);affineError=float(abs(got-(U+points@G.T)).max());assert constantError<=4e-15 and affineError<=1e-12
 rawArea=float(raw['area'].sum());personArea=np.array([area.reshape(-1,3)[person==j,0].sum()for j in range(len(p['offsets_m']))]);assert np.allclose(personArea,rawArea,rtol=2e-14,atol=2e-14)
 with(W.H/'INVERTED8_TRAINING.npz').open('xb')as f:np.savez(f,points=points,wall=wall,normal=normal,area=area,person=person,component=component,indices=ids,weights=w)
 training=W.pinfile('INVERTED8_TRAINING.npz');pc=make(np,faces,1.225);save(np,pc,W.H/'TENSOR_PC.npz');blocks=make_blocks(np,faces,points,person,ids,w,md,DT,10000.,area);save_blocks(np,blocks,W.H/'SURFACE_BLOCK_PC.npz');action,C,Ct=make_action(np,pc,blocks,n,shapes,off,areas,md,ids,w)
 x=1e-4*np.sin(np.arange(N+K)*.001731);y=7e-5*np.cos(np.arange(N+K)*.002391);Mx=action(x);My=action(y);energy=float(x@Mx);sym=float(abs(x@My-y@Mx));adj=float(abs(x[:N]@C(y[N:])-Ct(x[:N])@y[N:]));assert energy>0 and sym<=1e-12 and adj<=1e-12
 r=dict(PASS=True,formation='inverted triangle front 2+2+1',training=training,tensorPC=W.pinfile('TENSOR_PC.npz'),blockPC=W.pinfile('SURFACE_BLOCK_PC.npz'),K=K,N=N,D=D,counts=counts,rawExposedPoints=len(raw['points']),quadratureN=47,areaPerPerson_m2=personArea.tolist(),constantComponentMaximumError=constantError,affineMaximumError=affineError,PCPositiveEnergy=energy,PCSymmetryError=sym,couplingAdjointError=adj,poseSHA=p['pose']['SHA'],sourceOffsets_m=p['offsets_m'],allPositiveOriginalExposedAreaPreserved=True,oldN64QuadratureNotClaimedEqual=True,newMotionQueries=0,newCFDSteps=0,read_B=W.READ,wall_s=time.monotonic()-t0)
 W.save('TRAINING_AND_PC_PASS.json',r);print(json.dumps(r),flush=True)
 
def prepare():
 meta=json.loads((W.H/'TRAINING_AND_PC_PASS.json').read_text('utf8'));assert meta['PASS']and meta['sourceOffsets_m']==p['offsets_m']
 faces,n,N,h,shapes,sizes,off,D,mass,areas,md=W.grid(np,p['grid'])
 with np.load(io.BytesIO(W.frozen(meta['training'])),allow_pickle=False)as z:data={k:z[k].copy()for k in z.files}
 ids=np.ascontiguousarray(data['indices'],dtype='i8');w=np.ascontiguousarray(data['weights']);area=data['area'];wall=data['wall'];K=len(ids);L=N+K
 old=np.load(io.BytesIO(W.frozen(p['initialMomentum'])),allow_pickle=False);assert old.shape==(D,)and old.dtype==np.dtype('f8')and W.sha(old.tobytes())==p['initialMomentum']['momentumPayloadSHA'];zero=np.zeros(D);V=W.shared_arrays(np,D,N,K,frame['maps'],'prepare');prepared=V['prepared'];stats=V['predictorStats'];prepared.fill(-37.);stats.fill(-37.);V['output'].fill(-37.);V['stats'].fill(-37.);V['state'].fill(0.);V['trace'].fill(np.nan)
 W.save('ACTUAL_PREDICTOR_ENTRY.json',dict(beforeActualPredictor=True,once=True,formation='inverted',training=meta['training'],oldP=p['initialMomentum'],dtExact=p['dtExact'],pose=p['pose'],operator=p['fineDLL']))
 maskmeta=json.loads((W.H/'SGS_FLUID_MASK_PASS.json').read_text('utf8'));assert maskmeta['poseSHA']==p['pose']['SHA']
 fluid=np.load(io.BytesIO(W.frozen(maskmeta['asset'])),allow_pickle=False);assert fluid.dtype==np.dtype('u1')and fluid.shape==(N,)
 nut=np.zeros(N);sgsstats=np.zeros(9)
 lib=W.dll(p['fineDLL']);fn=lib.pm_prepare_les;fn.argtypes=[W.DP]*3+[C.c_double]*5+[C.c_int,W.DP,C.POINTER(C.c_ubyte),W.DP,W.DP,W.DP,W.DP];fn.restype=C.c_int
 code=fn(*[W.ptr(q)for q in faces],1.225,DT,1.5e-5,-2813/480,.035,1,W.ptr(old),fluid.ctypes.data_as(C.POINTER(C.c_ubyte)),W.ptr(prepared),W.ptr(nut),W.ptr(stats),W.ptr(sgsstats));assert code==0,('Actual LES predictor refused',code)
 assert np.isfinite(nut).all()and np.all(nut>=0)and np.all(nut[fluid==0]==0)
 W.save('SGS_PREDICTOR_PASS.json',dict(PASS=True,coefficient=.035,stats=sgsstats.tolist(),nut=W.np_save('SGS_NUT.npy',nut,np),fluidMask=maskmeta['asset'],poseSHA=p['pose']['SHA']))
 report=(C.c_ulonglong*4)();lib.pm_parallel_report(report)
 assert list(report)[0]==255,('Actual transport OpenMP mask',list(report))
 W.save('TRANSPORT_PARALLEL_EVIDENCE.json',dict(actualMasks=list(report),solver=p['fineDLL']))
 assert np.array_equal(md,prepared[8*D:9*D])and float(stats[0]+stats[1])<=.5
 star=prepared[5*D:6*D];canonical=np.zeros(n)
 for a in range(3):canonical+=np.diff(areas[a]*(star[off[a]:off[a+1]].reshape(shapes[a])/mass[a]),axis=a)
 assert np.array_equal((-canonical/(DT)).ravel(),prepared[10*D:]);rhs=V['rhs'];rhs.fill(0.);rhs[:N]=prepared[10*D:];ju=np.zeros(K)
 for k in range(24):ju+=w[:,k]*(star[ids[:,k]]/md[ids[:,k]])
 rhs[N:]=(ju-wall)/(DT);assert np.isfinite(rhs).all()and np.linalg.norm(rhs)>0
 parts=[];rp=W.np_save('COUPLED_RHS.npy',rhs,np)
 raw=W.sha(rhs.tobytes());signature=W.sha(json.dumps(dict(rhs=raw,training=meta['training']['SHA'],pose=p['pose']['SHA'],old=p['initialMomentum']['SHA'],offsets=p['offsets_m'],dt=p['dtExact'],dll=p['fineDLL']['SHA']),sort_keys=True).encode())
 r=dict(PASS=True,formation='inverted 2+2+1',K=K,N=N,D=D,L=L,parts=parts,preparedSHA=W.sha(prepared.tobytes()),rhs=rp,rawBSHA=raw,initialPressureRHS_SHA=W.sha(rhs[:N].tobytes()),transactionId=signature,training=meta['training'],tensorPC=meta['tensorPC'],blockPC=meta['blockPC'],oldMomentumPayloadSHA=W.sha(old.tobytes()),MdSHA=W.sha(md.tobytes()),CFL=float(stats[0]),diffusion=float(stats[1]),timeExact=p['timeExact'],dtExact=p['dtExact'],poseSHA=p['pose']['SHA'],nativePredictorReturnCode=code,physicalArraysComputed=True,newCFDPredictorCalls=1,newCoupledPressureCalls=0,read_B=W.READ,wall_s=time.monotonic()-t0)
 W.save('ACTUAL_PREPARED_PASS.json',r);print(json.dumps({k:v for k,v in r.items()if k!='parts'}),flush=True)
 
def pressure():
 chunk=frame['chunk'];previousSHA=frame['previous']['SHA']if frame.get('previous')else'NONE'
 meta=json.loads((W.H/'ACTUAL_PREPARED_PASS.json').read_text('utf8'));assert meta['PASS']and meta['nativePredictorReturnCode']==0
 faces,n,N,h,shapes,sizes,off,D,mass,areas,md=W.grid(np,p['grid'])
 with np.load(io.BytesIO(W.frozen(meta['training'])),allow_pickle=False)as z:data={k:z[k].copy()for k in z.files}
 ids=np.ascontiguousarray(data['indices'],dtype='i8');w=np.ascontiguousarray(data['weights']);area=np.ascontiguousarray(data['area']);wall=np.ascontiguousarray(data['wall']);K=len(ids);L=N+K
 V=W.shared_arrays(np,D,N,K,frame['maps'],'pressure');prepared=V['prepared'];rhs=V['rhs'];assert W.sha(prepared.tobytes())==meta['preparedSHA']and W.sha(rhs.tobytes())==meta['rawBSHA'];prepared.flags.writeable=False;rhs.flags.writeable=False
 state=V['state'];trace=V['trace']
 if chunk==0:
  assert previousSHA=='NONE'and np.count_nonzero(state)==0;trace.fill(np.nan)
 else:
  body=W.frozen(frame['previous']);assert W.sha(body)==previousSHA;previous=json.loads(body);assert previous['status']==101 and previous['chunk']==chunk-1 and previous['rawBSHA']==meta['rawBSHA']and previous['transactionId']==meta['transactionId']
  assert W.sha(state.tobytes())==previous['statePayloadSHA']and int(state[1])==previous['iteration']
 assert state.dtype==np.dtype('f8')and state.shape==(3*L+32,);initial=int(state[1]);output=V['output'];stats=V['stats']
 from tensor_pressure_pc import load
 from surface_block_preconditioner import load_blocks,make_action
 with np.load(io.BytesIO(W.frozen(meta['tensorPC'])),allow_pickle=False)as z:pc=load(np,z)
 with np.load(io.BytesIO(W.frozen(meta['blockPC'])),allow_pickle=False)as z:blocks=load_blocks(np,z)
 assert blocks['K']==K;action,_,_=make_action(np,pc,blocks,n,shapes,off,areas,md,ids,w);PC=C.CFUNCTYPE(C.c_int,C.c_void_p,C.c_void_p,C.c_void_p);errors=[];calls=[]
 @PC
 def callback(ip,op,ctx):
  try:
   src=np.ctypeslib.as_array((C.c_double*L).from_address(ip));dst=np.ctypeslib.as_array((C.c_double*L).from_address(op));dst[:]=action(src);calls.append(1);return 0
  except BaseException:errors.append(traceback.format_exc());return 1
 lib=W.dll(p['fineDLL']);abi=(C.c_int*4)();assert lib.pm_quad_precision(abi)==0 and list(abi)==[8,53,1794426759,7]
 fn=lib.pm_surface_quad_project_slice;fn.argtypes=[W.DP]*3+[C.c_double]*3+[W.DP,C.c_int,W.LP,W.DP,W.DP,W.DP,W.DP,W.DP,W.DP,W.DP,W.DP,C.c_int,PC,C.c_void_p];fn.restype=C.c_int
 W.save('PRESSURE_ENTRY_%03d.json'%chunk,dict(beforeCAPI=True,once=True,chunk=chunk,iterationStart=initial,transactionId=meta['transactionId'],rawBSHA=meta['rawBSHA'],previousSuccessfulSliceSHA=previousSHA,maxAdvance=64,globalCG384=True,ABI=list(abi),newRHS=False))
 for advance in range(1):
  code=fn(*[W.ptr(q)for q in faces],1.225,DT,10000.,W.ptr(prepared),K,ids.ctypes.data_as(W.LP),W.ptr(w),W.ptr(area),W.ptr(wall),W.ptr(rhs),W.ptr(output),W.ptr(trace),W.ptr(stats),W.ptr(state),64,callback,None)
  if code!=101:break
  assert state[5]==0 and 0<int(state[1])<=384

 report=(C.c_ulonglong*4)();lib.pm_parallel_report(report)
 assert list(report)[1:3]==[255,255],('Actual CG/operator OpenMP mask',list(report))
 W.save('PRESSURE_PARALLEL_EVIDENCE_%03d.json'%chunk,dict(actualMasks=list(report),solver=p['fineDLL']))
 iteration=int(state[1]);good=code in(0,101)and not errors;rec=dict(PASS=good,status=code,chunk=chunk,iterationStart=initial,iteration=iteration,rawBSHA=meta['rawBSHA'],transactionId=meta['transactionId'],callbackCalls=len(calls),callbackErrors=errors,wall_s=time.monotonic()-t0,read_B=W.READ,checkpoint=None,statePayloadSHA=W.sha(state.tobytes()),trueFullRaw=float(trace[iteration,1])if np.isfinite(trace[iteration,1])else None,truePressureRaw=float(trace[iteration,2])if np.isfinite(trace[iteration,2])else None,tauFull=float(trace[iteration,4])if np.isfinite(trace[iteration,4])else None,tauPressure=float(trace[iteration,5])if np.isfinite(trace[iteration,5])else None)
 if good:
  assert state[0]==1794426759 and list(state[6:11])==list(n)+[K,7]and initial<iteration<=initial+64 and iteration<=384
  assert W.sha(prepared.tobytes())==meta['preparedSHA']and W.sha(rhs.tobytes())==meta['rawBSHA']
  if code==101:assert state[5]==0 and np.all(output==-37.)and np.all(stats==-37.)
  else:assert state[5]==1;rec['candidateMemorySHA']=W.sha(output.tobytes());rec['stats']=stats.tolist()
  if code==0:
   name='CHECKPOINT_TERMINAL.npz'
   with(W.H/name).open('xb')as f:np.savez(f,state=state,trace=trace)
   rec['checkpoint']=W.pinfile(name)
 if not good:
  with(W.H/('FAILED_KRYLOV_%03d.npz'%chunk)).open('xb')as f:np.savez(f,state=state,trace=trace)
 W.save('SLICE_%03d.json'%chunk,rec)
 assert good,('Actual coupled pressure refused',code,errors)
 
def audit():
 done=json.loads((W.H/'PRESSURE_TERMINAL_COMPLETE.json').read_text('utf8'));assert done['PASS']and done['allOwnedHandlesClosed'];last=json.loads(W.frozen(done['lastSlice']));assert last['status']==0 and last['PASS'];meta=json.loads((W.H/'ACTUAL_PREPARED_PASS.json').read_text('utf8'))
 faces,n,N,h,shapes,sizes,off,D,mass,areas,md=W.grid(np,p['grid'])
 with np.load(io.BytesIO(W.frozen(meta['training'])),allow_pickle=False)as z:data={k:z[k].copy()for k in z.files}
 ids=data['indices'];w=data['weights'];a=data['area'];wall=data['wall'];K=len(ids);L=N+K
 V=W.shared_arrays(np,D,N,K,frame['maps'],'audit');out=V['output'];rhs=V['rhs'];prepared=V['prepared'];assert W.sha(out.tobytes())==last['candidateMemorySHA']and W.sha(prepared.tobytes())==meta['preparedSHA'];assert out.shape==(3*D+2*L,)and np.isfinite(out).all();assert W.sha(rhs.tobytes())==meta['rawBSHA'];x=out[3*D:3*D+L]+out[3*D+L:];reaction=x[N:];pr=x[:N].reshape(n);spread=np.zeros(D);np.add.at(spread,ids.ravel(),(w*reaction[:,None]).ravel());gradient=np.zeros(D)
 for ax in range(3):
  mid=[slice(None)]*3;mid[ax]=slice(1,-1);left=[slice(None)]*3;right=left.copy();left[ax]=slice(None,-1);right[ax]=slice(1,None);gradient[off[ax]:off[ax+1]].reshape(shapes[ax])[tuple(mid)]=areas[ax][tuple(mid)]*(pr[tuple(right)]-pr[tuple(left)])
 gradient[:sizes[0]].reshape(shapes[0])[0]=areas[0][0]*pr[0]
 def B(value):
  r=np.zeros(n)
  for ax in range(3):r-=np.diff(areas[ax]*value[off[ax]:off[ax+1]].reshape(shapes[ax]),axis=ax)
  return r.ravel()
 def J(value):
  r=np.zeros(K)
  for k in range(24):r+=w[:,k]*value[ids[:,k]]
  return r
 force=gradient+spread;velocity=out[:D]/md;action=np.r_[B(force/md),J(force/md)+reaction/((DT)*10000*a)];res=rhs-action;norm=lambda x:np.sqrt(np.sum(x*x));tau=max(1e-10,1e-5*norm(rhs));effective=rhs[:N]-B(spread/md);taup=max(1e-10,1e-5*norm(effective));raw=float(norm(res));praw=float(norm(res[:N]));assert raw<=tau and praw<=taup
 volume=h[0][:,None,None]*h[1][None,:,None]*h[2][None,None,:];dc=np.zeros(n);dv=np.zeros(n)
 for ax in range(3):
  v=out[off[ax]:off[ax+1]].reshape(shapes[ax])/mass[ax];q=[1]*3;q[ax]=n[ax];dc+=np.diff(v,axis=ax)/h[ax].reshape(q);dv+=np.diff(areas[ax]*v,axis=ax)
 both=[float(abs(dc).max()),float(abs(dv/volume).max())];U=2813/480;normalized=max(both)*min(float(q.min())for q in h)/U;dtdiv=max(both)*DT;vmax=float(abs(velocity).max());assert normalized<=1e-4 and dtdiv<=1e-4 and vmax<=10*U
 star=prepared[5*D:6*D];pt=prepared[3*D:4*D];js=prepared[7*D:8*D];momentum=float(abs(out[:D]-pt+(DT)*force-js).max());projection=float(abs(out[:D]-star+(DT)*force).max());energy=float(np.sum(out[:D]*out[:D]/(2*md)));old=np.load(io.BytesIO(W.frozen(p['initialMomentum'])),allow_pickle=False);oldEnergy=float(np.sum(old*old/(2*md)));slip=(J(velocity)-wall).reshape(-1,3);pointArea=a.reshape(-1,3)[:,0];normalSlip=np.sum(slip*data['normal'],axis=1);slipRMS=float(np.sqrt(np.sum(pointArea*normalSlip**2)/pointArea.sum()));assert np.isfinite([energy,oldEnergy,momentum,projection,slipRMS]).all()
 newpin=W.np_save('INVERTED8_newP.npy',out[:D],np);newpin['momentumPayloadSHA']=W.sha(out[:D].tobytes())
 with(W.H/'INVERTED8_pressure_reaction_HiLo.npz').open('xb')as f:np.savez(f,hi=out[3*D:3*D+L],lo=out[3*D+L:])
 r=dict(PASS=True,acceptance='EXPLORATORY_STABILITY_ONLY',formation='inverted triangle front: 2+2+1; protected P0; two rear',offsets_m=p['offsets_m'],speed_km_h=21.0975,cadence_spm=192,gridShape=list(n),N=N,D=D,K=K,timeExact=p['timeExact'],dtExact=p['dtExact'],poseSHA=p['pose']['SHA'],newP=newpin,pressureHiLo=W.pinfile('INVERTED8_pressure_reaction_HiLo.npz'),iterations=done['iterations'],trueFullRaw=raw,tauFull=float(tau),truePressureRaw=praw,tauPressure=float(taup),normalizedContinuity=normalized,bothStoredDiv=both,dtDiv=dtdiv,maximumFaceVelocity_m_s=vmax,kineticEnergy_J=energy,previousKineticEnergy_J=oldEnergy,sourceMomentumDiagnostic=momentum,projectionDiagnostic=projection,normalSlipRMS_m_s=slipRMS,bodyReactionIncludingFictitiousFluid_N={str(j):reaction.reshape(-1,3)[data['person']==j].sum(0).tolist()for j in range(len(p['offsets_m']))},bodyReactionIsNotPhysicalDrag=True,transientInitializationFromSolo=False,quadratureN=47,oldN64DirectForceComparisonQualified=False,physicalTotalDragQualified=False,physicalImpermeabilityQualified=False,steadyOrPeriodicQualified=False,parallelSolver=p['fineDLL'],explicitSGSClosureIncluded=True,SGSModel='Vreman C=0.035',movingGround=True,newCFDTimeSteps=1,allPressureOwnedActorsClosed=True,rawBSHA=meta['rawBSHA'],transactionId=meta['transactionId'],read_B=W.READ,wall_s=time.monotonic()-t0)
 W.save('INVERTED8_EXPLORATORY_ACTUAL_RESULT.json',r);print(json.dumps(r),flush=True)
 
try:
 if mode=='cfl':cfl_plan()
 elif mode=='motion':motion()
 elif mode=='geometry':geometry()
 elif mode=='sgs_mask':sgs_mask()
 elif mode=='training':training()
 elif mode=='prepare':prepare()
 elif mode=='pressure':pressure()
 elif mode=='audit':audit()
 else:raise RuntimeError('Unadmitted role')
finally:W.close_shared()
