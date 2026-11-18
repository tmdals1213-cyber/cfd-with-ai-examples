#!/usr/bin/env python3
"""Build/prepare legacy numerical sources. This tool never launches a flow solver."""
from pathlib import Path
import argparse,json,subprocess,shutil,hashlib,os,time,ctypes,math,itertools
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def catalog():return json.loads((ROOT/'examples/legacy_catalog.json').read_text())
def kernel_test(libpath):
 lib=ctypes.CDLL(str(libpath));DP=ctypes.POINTER(ctypes.c_double);IP=ctypes.POINTER(ctypes.c_int)
 dims=(ctypes.c_int*3)();lib.pm_dimensions.argtypes=[IP];assert lib.pm_dimensions(dims)==0 and list(dims)==[5,7,9]
 precision=(ctypes.c_int*4)();lib.pm_quad_precision.argtypes=[IP];assert lib.pm_quad_precision(precision)==0 and precision[0]==8 and precision[1]==53
 grad=lib.pm_vreman_gradient;grad.argtypes=[DP,DP,ctypes.c_double];grad.restype=ctypes.c_double
 h=(ctypes.c_double*3)(.1,.15,.2);a=(ctypes.c_double*9)(1,0,0,0,2,0,0,0,3)
 beta=[.1**2,(2*.15)**2,(3*.2)**2];truth=.035*math.sqrt((beta[0]*beta[1]+beta[0]*beta[2]+beta[1]*beta[2])/14)
 assert abs(grad(a,h,.035)-truth)<1e-15
 zero=(ctypes.c_double*9)();assert grad(zero,h,.035)==0
 shear=(ctypes.c_double*9)(0,1,0,0,0,0,0,0,0);assert grad(shear,h,.035)==0
 return {'dimensions':list(dims),'doublePrecisionABI':True,'VremanAnalyticChecks':3,'fluidStepsAdvanced':0}
def build(id,compiler):
 e=catalog()[id];source=ROOT/e['directory'];b=ROOT/'build/legacy'/id;b.mkdir(parents=True,exist_ok=False);start=time.monotonic()
 if e['backend']=='LICA-Fortran':
  shutil.copytree(source/'src',b/'src');files=[str(Path('src')/n)for n in e['sourceOrder']]
  flags=['-std=legacy','-fdefault-integer-8','-fdefault-real-8','-fdefault-double-8','-ffree-line-length-none','-O2','-fopenmp']
  if os.name=='nt':flags+=['-Wl,--stack,134217728']
  targets=[('lica.exe'if os.name=='nt'else'lica',flags+files)]
 elif id=='15_runner_final':
  shutil.copytree(source/'src',b/'src');flags=['-O2','-std=c++17','-shared','-fopenmp','-fno-fast-math','-ffp-contract=off']
  flags+=['-static']if os.name=='nt'else['-fPIC']
  ext='.dll'if os.name=='nt'else'.so';targets=[]
  for label,shape in [('cold',[5,7,9]),('preview',[76,29,72])]:targets.append(('runner_'+label+ext,flags+[f'-D{k}={v}'for k,v in zip(['NX','NY','NZ'],shape)]+['src/penalty_les.cpp']))
 else:raise ValueError('Reference-only entry has no asserted complete build recipe')
 records=[]
 for name,args in targets:
  cmd=[compiler,*args,'-o',name]
  with(b/(name+'.build.log')).open('w')as log:result=subprocess.run(cmd,cwd=b,stdout=log,stderr=subprocess.STDOUT)
  assert result.returncode==0,'Compile/link failed: '+str(b/(name+'.build.log'))
  records.append({'target':name,'binarySHA256':sha(b/name),'argv':cmd,'exitCode':result.returncode})
 d={'state':'BUILD_PASSED_NO_CFD_LAUNCHED','case':id,'source':{p.name:sha(p)for p in(source/'src').iterdir()if p.is_file()},'targets':records,'seconds':time.monotonic()-start}
 if id=='15_runner_final':d['kernelChecks']=kernel_test(b/targets[0][0])
 (b/'BUILD.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'case':id,'state':d['state'],'seconds':d['seconds'],'kernelChecks':d.get('kernelChecks')}))
def prepare(id,variant):
 e=catalog()[id];source=ROOT/e['directory']/'inputs'/variant;assert source.is_dir(),'Select a retained input variant'
 out=ROOT/'runs/legacy'/(time.strftime('%Y%m%d_%H%M%S')+'_'+id+'_'+variant);out.parent.mkdir(parents=True,exist_ok=True);shutil.copytree(source,out)
 p=out/'lica.in';text=p.read_text();rows=text.splitlines();found=False
 for i,line in enumerate(rows):
  if line.strip().startswith('RE IRESET IREAD '):
   values=rows[i+1].split();assert len(values)>=6;values[1]='1';values[2]='0';rows[i+1]=' '.join(values);found=True
 assert found,'Unknown LICA input schema; manual review required'
 p.write_text('\n'.join(rows)+'\n')
 missing=['matching IBM preprocessor output for this body/grid and kind width']
 if 'mod_stl_sdf.f90'in e.get('sourceOrder',[]):missing.append('body.sdf from the correctly oriented/scaled licensed geometry')
 (out/'PREPARATION.json').write_text(json.dumps({'state':'SETUP_ONLY_NOT_RUN','case':id,'variant':variant,'coldStartIREAD':0,'IRESET':1,'missingExternalInputs':missing,'noRawFieldBorrowed':True,'noFlowSolverLaunched':True},indent=2)+'\n')
 print(json.dumps({'prepared':str(out),'missing':missing,'CFDlaunched':False}))
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['list','build','prepare']);ap.add_argument('case',nargs='?');ap.add_argument('--compiler');ap.add_argument('--jobs',type=int,default=1,help='Single compiler invocation; Fortran module order is preserved.');ap.add_argument('--variant',default='default');a=ap.parse_args()
 if a.action=='list':
  for k,v in catalog().items():print(k,v.get('backend',v.get('state')),v['directory'])
 elif a.action=='build':build(a.case,a.compiler or('g++'if a.case=='15_runner_final'else'gfortran'))
 else:prepare(a.case,a.variant)
if __name__=='__main__':main()
