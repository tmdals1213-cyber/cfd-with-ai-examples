from exposed_surface import *
from media_common import save
import time
t0=time.monotonic();request=Path(sys.argv[1]);req=json.loads(request.read_text('utf8'));O=request.parent
assert pin(Path(__file__))==req['worker']and pin(R/'source/pressure_kernel.py')==req['kernel']and pin(R/'source/exposed_surface.py')==req['geometryWorker']
cfg=json.loads(frozen(req['config']));result=json.loads(frozen(req['entry']['result']));assert result['PASS']and result['explicitSGSClosureIncluded'];folder=Path(req['entry']['result']['path']).parent
motion=json.loads((folder/'MOTION_PASS.json').read_text());assert motion['pose']['SHA']==result['poseSHA']
if req['case']=='group8':assert motion['minimumAABBSeparatingGap_m']>0
faces,cent,n,shapes,off,md=grid_info(cfg['grid']);points,normal,area,volume,gcontrol,xyz=exposed(cfg,motion['pose'],95);points+=np.array(cfg['offsets_m'][0]);weights,wcontrol=moment_weights(points,normal,area,volume)
N=int(np.prod(n))
with np.load(io.BytesIO(frozen(result['pressureHiLo'])),allow_pickle=False)as z:p=(z['hi'][:N]+z['lo'][:N]).reshape(n)
pv=sample_pressure(p,points,cent);pnx=pv*normal[:,0];load=float(np.sum(pnx*weights));gauge=float(np.sum((pv+137)*normal[:,0]*weights));assert abs(gauge-load)<1e-9
with(O/'SURFACE.npz').open('xb')as f:np.savez_compressed(f,points=points,normal=normal,area=weights,pnx=pnx,body=xyz)
out=dict(PASS=True,case=req['case'],time_s=req['entry']['time_s'],timeExact=result['timeExact'],acceptedFluidStep=req['entry']['acceptedFluidStep'],pressureLoadEstimate_N=load,
 sourceResult=req['entry']['result'],sourcePressure=result['pressureHiLo'],sourcePose=motion['pose'],surface=pin(O/'SURFACE.npz'),geometry=gcontrol,quadrature=wcontrol,
 originalUnfittedPressureLoad_N=float(np.sum(pnx*area)),gaugeOffsetCheckError_N=abs(gauge-load),
 interiorAndOverlapFacesExcluded=True,physicalPressureExteriorLimitQualified=False,pressureCanUseSolidSideCells=True,physicalTotalDragQualified=False,viscousForceIncluded=False,
 modelEstimate=True,newCFDSteps=0,newMotionQueries=0,seconds=time.monotonic()-t0)
save(O/'PRESSURE_LOAD.json',out);print(json.dumps(dict(PASS=True,case=req['case'],time_s=out['time_s'],load_N=load,maxAreaCorrection=wcontrol['maximumRelativeAreaCorrection'],seconds=out['seconds'])))
