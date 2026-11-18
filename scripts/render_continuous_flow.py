"""Actual accepted XYZ PLIC + temperature/velocity sequence. No interpolated CFD states."""
from pathlib import Path
import argparse,csv,io,json,hashlib,os,subprocess,time,sys,importlib.util
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from PIL import Image,ImageDraw,ImageFont

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',required=True,type=Path);ap.add_argument('--helper',type=Path,default=Path(__file__).with_name('render_showcase.py'));ap.add_argument('--output',required=True,type=Path);ap.add_argument('--ffmpeg',required=True);ap.add_argument('--limit',type=float,required=True);ap.add_argument('--seconds',type=float,default=12);ap.add_argument('--form',choices=['landscape','portrait'],default='landscape');a=ap.parse_args()
 spec=importlib.util.spec_from_file_location('showcase',a.helper);h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
 source=a.run;out=a.output;out.mkdir(parents=True,exist_ok=False);assets=out/'assets';assets.mkdir()
 def load(p):return h.stable(p)
 raw,hpin=load(source/'history.csv');history=[{k:float(v)for k,v in x.items()}for x in csv.DictReader(raw.decode().splitlines())];assert history
 assert all(all(np.isfinite(list(x.values())))for x in history)
 accepted=min(a.limit,history[-1]['time_s']);base=source/'ensight/VaporizingDrop'
 raw,cpin=load(base/'nga.case');times=[float(x)for x in raw.decode().split('time values:')[1].split()]
 raw,ppin=load(base/'plic.case');ptimes=[float(x)for x in raw.decode().split('time values:')[1].split()];assert times==ptimes,'PLIC and field schedules differ'
 params={x.partition(':')[0].strip():x.partition(':')[2].strip()for x in(source/'input').read_text().splitlines()if ':'in x}
 assert [int(params[k])for k in['nx','ny','nz']]==[80]*3 and all(abs(float(params[k])-.02)<1e-14 for k in['Lx','Ly','Lz'])
 cells=80**3;pins=[];states=[]
 def field(name,index,components=1):
  p=base/name/f'{name}.{index:06d}';raw,pin=load(p);assert len(raw)==244+4*cells*components
  v=np.frombuffer(raw[244:],'<f4').copy();assert np.isfinite(v).all(),name
  return (v.reshape(components,80,80,80)if components>1 else v.reshape(80,80,80)),pin
 entries=[]
 for index,t in enumerate(times,1):
  row=min(history,key=lambda x:abs(x['time_s']-t))
  if t>accepted+1e-9 or abs(row['time_s']-t)>1e-7:continue
  for key,lim in [('total_energy_error_J',1e-5),('pressure_relative_error',1.01e-8),('div_minus_source_s_1',1e-4)]:assert abs(row[key])<=lim
  if abs(row['phase_power_W'])>1e-8:assert abs(row['mass_balance_relative'])<=1e-4
  entries.append((index,row))
 assert len(entries)>=3,'At least three actual matched saved times needed'
 assert np.all(np.diff([row['time_s']for _,row in entries])>0),'Saved times must increase strictly'
 # Verify one exported complete field against all eight independent physical CSV packets.
 matches=[(i,x)for i,x in entries if all((source/f'field_step{int(x["step"]):08d}_rank{rank:04d}.csv').exists()for rank in range(8))];assert matches,'No complete CSV/EnSight cross-check pair'
 ci,crow=matches[-1];arrays=[];csvpins=[]
 for rank in range(8):
  raw,pin=load(source/f'field_step{int(crow["step"]):08d}_rank{rank:04d}.csv');v=np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1);assert v.shape==(64000,12)and np.isfinite(v).all();arrays.append(v);csvpins.append(pin)
 v=np.concatenate(arrays);ijk=np.rint(np.column_stack(((v[:,0]+.01)/.00025-.5,v[:,1]/.00025-.5,(v[:,2]+.01)/.00025-.5))).astype(int);assert np.all((ijk>=0)&(ijk<80));off=ijk[:,0]+80*(ijk[:,1]+80*ijk[:,2]);assert len(np.unique(off))==cells
 expected=np.empty((cells,12));expected[off]=v
 verified={}
 for name,col in [('temperature',3),('latent_sink',9),('VOF',4)]:
  actual,_=field(name,ci);want=expected[:,col].astype('<f4').reshape(80,80,80);assert np.array_equal(actual,want),name+' differs from actual doubleCSV';verified[name]=True
 actual,_=field('velocity',ci,3);want=expected[:,5:8].astype('<f4').T.reshape(3,80,80,80);assert np.array_equal(actual,want),'Velocity differs from actual doubleCSV';verified['velocity']=True
 solid=expected[:,8].reshape(80,80,80)>.5;gas=np.sum((1-expected[:,4])*(expected[:,8]<.5))*.00025**3;assert abs(gas-crow['gas_volume_m3'])<1e-13
 X,Y=np.meshgrid(np.linspace(-.009875,.009875,80),np.linspace(.000125,.019875,80));sl=np.zeros((80,80,12));sl[:,:,0]=X;sl[:,:,1]=Y;sl[:,:,2]=.000125;sl[:,:,8]=solid[40]
 import matplotlib
 fontdir=Path(matplotlib.get_data_path())/'fonts/ttf';font=lambda n,b=False:ImageFont.truetype(str(fontdir/('DejaVuSans-Bold.ttf'if b else'DejaVuSans.ttf')),n)
 portrait=a.form=='portrait';W,H=(1080,1920)if portrait else(1920,1080)
 def fit(c,im,rect):
  x,y,w,z=map(int,rect);q=im.copy();q.thumbnail((w,z),Image.Resampling.LANCZOS);c.paste(q,(x+(w-q.width)//2,y+(z-q.height)//2))
 saved=[];previous=None;maxchange=0.
 for j,(index,row)in enumerate(entries):
  T,tp=field('temperature',index);vel,vp=field('velocity',index,3);vf,fp=field('VOF',index);assert np.all((vf>=0)&(vf<=1));poly,pp=h.plic(base/'plic'/f'plic.{index:06d}')
  sl[:,:,3]=T[40];sl[:,:,4]=vf[40];sl[:,:,5:8]=np.moveaxis(vel[:,40],0,-1)
  im3=h.render3d(poly);imt=h.slice_image(sl,poly);imv=h.slice_image(sl,poly,flow=True)
  c=Image.new('RGB',(W,H),h.BG);d=ImageDraw.Draw(c);d.text((60,42),'3D IMMERSION COOLING',font=font(42,True),fill=h.FG);d.text((60,103),f"t = {row['time_s']*1000:.3f} ms",font=font(30,True),fill=h.TEAL)
  if portrait:
   fit(c,im3,(45,190,990,870));fit(c,imt,(12,1120,530,680));fit(c,imv,(538,1120,530,680));labels=[(60,180,'LIQUID / VAPOR'),(60,1080,'TEMPERATURE'),(570,1080,'VELOCITY')]
  else:
   fit(c,im3,(0,195,910,795));fit(c,imt,(900,210,510,735));fit(c,imv,(1400,210,510,735));labels=[(60,185,'LIQUID / VAPOR'),(945,185,'TEMPERATURE'),(1450,185,'VELOCITY')]
  for x,y,label in labels:d.text((x,y),label,font=font(22,True),fill=h.MUTED)
  d.text((60,H-55),'ACTUAL SAVED STATES  /  SLOW-MOTION PLAYBACK',font=font(20),fill=h.MUTED)
  path=assets/f'frame_{j:04d}.png';c.save(path);saved.append(c)
  if previous is not None:maxchange=max(maxchange,float(np.max(np.abs(vel-previous))))
  previous=vel.copy();pins.append({'step':int(row['step']),'time_s':row['time_s'],'temperature':tp,'velocity':vp,'VOF':fp,'PLIC':pp,'maxNativeSpeed_m_s':float(np.linalg.norm(vel,axis=0).max())});print(json.dumps({'frame':j,'time_s':row['time_s']}),flush=True)
 target=out/f'CONTINUOUS_XYZ_FLOW_{a.form.upper()}.mp4';fps=24;count=round(a.seconds*fps);native=np.array([x['time_s']for _,x in entries]);indices=np.abs(native[:,None]-np.linspace(native[0],native[-1],count)[None,:]).argmin(axis=0)
 with(out/'encoder.log').open('wb')as log:
  proc=subprocess.Popen([a.ffmpeg,'-hide_banner','-loglevel','error','-n','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r',str(fps),'-i','pipe:0','-an','-c:v','libx264','-threads','2','-preset','veryfast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(target)],stdin=subprocess.PIPE,stdout=subprocess.DEVNULL,stderr=log)
  try:
   for i in indices:proc.stdin.write(saved[i].tobytes())
   proc.stdin.close();assert proc.wait(timeout=60)==0
  finally:
   if proc.poll()is None:proc.kill();proc.wait()
 saved[-1].save(out/'POSTER.png');subprocess.run([a.ffmpeg,'-hide_banner','-v','error','-i',str(target),'-f','null','-'],check=True,timeout=60)
 receipt={'actualSavedStates':len(entries),'firstTime_s':float(native[0]),'lastTime_s':float(native[-1]),'frameFieldsInterpolated':False,'CFDstepsExecuted':0,'presentationFramesRepeatNearestActualSavedState':True,'sourceHistory':hpin,'sourceCase':cpin,'sourcePLICCase':ppin,'independentCSV':csvpins,'independentExportByteComparisons':verified,'gasVolumeError_m3':gas-crow['gas_volume_m3'],'maxAdjacentSavedVelocityChangeLinf_m_s':maxchange,'fixedTemperatureRange_C':[45,63],'fixedSpeedRange_m_s':[0,.6],'completeCoolingCycleQualified':False,'sourceSnapshots':pins,'video':{'name':target.name,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'bytes':target.stat().st_size,'decodedCompletely':True,'frames':count,'fps':fps}};(out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'output':str(target),'states':len(entries),'lastTime_s':float(native[-1])}))
if __name__=='__main__':main()
