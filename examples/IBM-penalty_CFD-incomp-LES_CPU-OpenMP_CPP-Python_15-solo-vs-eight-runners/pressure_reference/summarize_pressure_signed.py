"""Keep signed model pressure components; never force a positive drag result."""
from pressure_kernel import *
from media_common import save
state=json.loads((R/'SURFACE95_CLOSURE.json').read_text());assert state['PASS'];curves={};means={};checks={}
for case in('solo','group8'):
 rows=sorted([q for q in state['frames']if q['case']==case],key=lambda q:q['time_s']);t=np.array([q['time_s']for q in rows]);v=np.array([q['pressureLoadEstimate_N']for q in rows]);assert np.all(np.diff(t)>0)and t[0]<.625 and t[-1]>2.5
 a,b=1.25,2.5;inside=(t>a)&(t<b);tt=np.r_[a,t[inside],b];vv=np.interp(tt,t,v);means[case]=float(np.trapezoid(vv,tt)/(b-a));curves[case]=(t,v)
 checks[case]=dict(frames=len(rows),maximumRelativeAreaCorrection=max(q['quadrature']['maximumRelativeAreaCorrection']for q in rows),maximumGaugeError_N=max(q['gaugeOffsetCheckError_N']for q in rows),maximumAffinePressureIdentityError_m3=max(q['quadrature']['affinePressureIdentityError_m3']for q in rows),meanUnfitted_N=float(np.trapezoid(np.interp(tt,t,[q['originalUnfittedPressureLoad_N']for q in rows]),tt)/(b-a)))
assert means['solo']>0
delta=means['solo']-means['group8'];percent=100*delta/means['solo'];common=np.linspace(.625,2.5,751)
with(R/'PRESSURE_SERIES_SIGNED.npz').open('xb')as f:np.savez(f,soloTime=curves['solo'][0],soloLoad=curves['solo'][1],pacedTime=curves['group8'][0],pacedLoad=curves['group8'][1],commonTime=common,soloCommon=np.interp(common,*curves['solo']),pacedCommon=np.interp(common,*curves['group8']))
summary=dict(PASS=True,quantity='Signed pressure-only stored-field load estimate on exposed-union surface',meanWindow_s=[1.25,2.5],gaitCyclesInMean=2.,
 meanSolo_N=means['solo'],meanPaced_N=means['group8'],pressureLoadDifference_N=delta,algebraicPressureComponentChange_percent=percent,
 publishableDragReduction_percent=None,percentagePublicationHeld=True,
 percentageHoldReason='Paced model pressure component is near zero and slightly negative; fluid-side pressure extraction and missing viscous contribution prevent a total-drag reduction claim.',
 checks=checks,physicalTotalDragQualified=False,physicalPressureExteriorLimitQualified=False,steadyOrPeriodicQualified=False,viscousForceIncluded=False,
 pressureCanUseSolidSideCells=True,sourceNumericalFieldsUnchanged=True,negativeValuesClipped=False,oldRawPressureDiagnosticNotPromoted=True,
 limitation='Numerical pressure contribution at the regularized interface on the 0.10 m grid; not validated exterior pressure drag or pressure-plus-viscous total drag.',
 sourceHistory=pin(R/'SURFACE95_CLOSURE.json'),sourceSeries=pin(R/'PRESSURE_SERIES_SIGNED.npz'))
save(R/'PRESSURE_COMPARISON_SIGNED_SUMMARY.json',summary);print(json.dumps(summary))
