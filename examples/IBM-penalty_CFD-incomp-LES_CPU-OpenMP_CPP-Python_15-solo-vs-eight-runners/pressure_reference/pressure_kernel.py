"""Stored-pressure surface-load estimate on saved exposed-union quadrature.

No IBM reaction, fictitious-fluid load, viscous force, or old raw N is used.
Pressure reconstruction can use solid-side cells: this remains a model estimate,
not qualified physical exterior traction or total aerodynamic drag.
"""
from pathlib import Path
import sys,io,json,hashlib
R=Path(__file__).resolve().parent.parent;WS=R.parents[1];sys.dont_write_bytecode=True
sys.path.insert(0,str(WS/'work/pacemaker_promo_20261113/source'))
from media_common import np,frozen,pin,grid_info

def sample_pressure(p,points,cent):
 ii=[np.clip(np.searchsorted(cent[a],points[:,a],side='right')-1,0,len(cent[a])-2)for a in range(3)]
 s=[(points[:,a]-cent[a][ii[a]])/(cent[a][ii[a]+1]-cent[a][ii[a]])for a in range(3)]
 assert min(float(q.min())for q in s)>=-.500001 and max(float(q.max())for q in s)<=1.500001
 out=np.zeros(len(points))
 for a in(0,1):
  for b in(0,1):
   for c in(0,1):out+=p[ii[0]+a,ii[1]+b,ii[2]+c]*np.where(a,s[0],1-s[0])*np.where(b,s[1],1-s[1])*np.where(c,s[2],1-s[2])
 return out

def moment_weights(points,normal,area,volume):
 # Minimum relative correction of positive saved quadrature weights to enforce
 # closed-surface constant/affine-pressure identities. No source area is edited.
 center=.5*(points.min(0)+points.max(0));length=float(np.ptp(points,axis=0).max());x=(points-center)/length
 C=np.vstack([np.ones(len(area)),normal.T,*[x[:,i]*normal[:,j]for i in range(3)for j in range(3)]])
 b=np.r_[area.sum(),np.zeros(3),(np.eye(3)*volume/length).ravel()]
 G=(C*area)@C.T;assert np.linalg.cond(G)<1e10
 correction=C.T@np.linalg.solve(G,b-C@area);assert abs(correction).max()<.25
 weight=area*(1+correction);assert np.all(weight>0)and np.max(abs(C@weight-b))<1e-10
 uniform=np.sum(normal*weight[:,None],axis=0)
 affine=points.T@(normal*weight[:,None]);expected=np.eye(3)*volume
 assert abs(uniform).max()<1e-10 and abs(affine-expected).max()<1e-10
 return weight,dict(maximumRelativeAreaCorrection=float(abs(correction).max()),normalClosure_m2=uniform.tolist(),affinePressureIdentityError_m3=float(abs(affine-expected).max()),originalArea_m2=float(area.sum()),fittedArea_m2=float(weight.sum()),volume_m3=volume)

def evaluate(entry,cent,n,offset):
 source=json.loads(frozen(entry['result']));assert source['PASS']and source['explicitSGSClosureIncluded']
 folder=Path(entry['result']['path']).parent;motion=json.loads((folder/'MOTION_PASS.json').read_text('utf8'));assert motion['pose']['SHA']==source['poseSHA']
 if len(source['offsets_m'])>1:assert motion['minimumAABBSeparatingGap_m']>0
 parts=[];volumes=[];pins=[]
 for a in range(3):
  meta=json.loads((folder/('EXPOSED47_AXIS%d_PASS.json'%a)).read_text());assert meta['PASS']and meta['sourcePoseSHA']==source['poseSHA']and meta['N']==47
  with np.load(io.BytesIO(frozen(meta['asset'])),allow_pickle=False)as z:parts.append({k:z[k].copy()for k in z.files})
  volumes.append(meta['volumeMidpointEstimate_m3']);pins.append(meta['asset'])
 data={k:np.concatenate([q[k]for q in parts])for k in parts[0]};points=data['points']+offset;normal=data['normal'];area=data['area'];assert np.all(area>0)and abs(np.linalg.norm(normal,axis=1)-1).max()<1e-10
 volume=float(np.mean(volumes));weights,control=moment_weights(points,normal,area,volume)
 N=int(np.prod(n))
 with np.load(io.BytesIO(frozen(source['pressureHiLo'])),allow_pickle=False)as z:p=(z['hi'][:N]+z['lo'][:N]).reshape(n)
 pressure=sample_pressure(p,points,cent);quantity=pressure*normal[:,0]
 force=float(np.sum(quantity*weights));gauge=float(np.sum((pressure+137)*normal[:,0]*weights));assert abs(force-gauge)<1e-9
 return dict(timeExact=source['timeExact'],time_s=entry['time_s'],acceptedFluidStep=entry['acceptedFluidStep'],pressureLoadEstimate_N=force,
  control=control,sourcePressure=source['pressureHiLo'],sourceResult=entry['result'],sourcePose=motion['pose'],sourceExposedSurface=pins,
  originalUnfittedPressureLoad_N=float(np.sum(quantity*area)),gaugeOffsetCheckError_N=abs(force-gauge),
  quantityMinimum_Pa=float(quantity.min()),quantityMaximum_Pa=float(quantity.max()),
  duplicateInteriorGeometryExcludedBySavedUnion=True,pressureExteriorLimitQualified=False,physicalTotalDragQualified=False,
  pressureCanUseSolidSideCells=True,viscousForceIncluded=False),dict(points=points,normal=normal,area=weights,pnx=quantity)
