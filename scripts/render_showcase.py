#!/usr/bin/env python3
"""Render a frozen actual XYZ PLIC/field selection into two clean 24 s layouts."""
from pathlib import Path
import argparse,io,csv,hashlib,json,os,subprocess,time,wave
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
from PIL import Image,ImageDraw,ImageFont
BG='#08131f';FG='#f1f6fa';MUTED='#92a9b9';TEAL='#4dddc6';COPPER='#ef9a4b';BLUE='#4d9bd1'
CUTS=[0,3,7,10,15,20,24]
def stable(p):
 a=p.stat();raw=p.read_bytes();b=p.stat();assert(a.st_size,a.st_mtime_ns)==(b.st_size,b.st_mtime_ns),'Input is changing: '+p.name
 return raw,{'name':p.name,'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
def plic(p):
 raw,pin=stable(p);f=io.BytesIO(raw)
 def txt():return f.read(80).decode('ascii').strip()
 def ints(n):return np.frombuffer(f.read(4*n),'<i4').copy()
 def reals(n):return np.frombuffer(f.read(4*n),'<f4').astype(float)
 assert txt()=='C Binary';[txt()for _ in range(4)];assert txt()=='extents';reals(6);out=[]
 while f.tell()<len(raw):
  assert txt()=='part';ints(1);txt();assert txt()=='coordinates';nv=int(ints(1)[0]);v=np.array([reals(nv),reals(nv),reals(nv)]).T
  assert txt()=='nsided';nf=int(ints(1)[0]);ss=ints(nf);conn=ints(int(ss.sum()))-1;offset=0
  for n in ss:
   if n>=3:out.append(v[conn[offset:offset+n]])
   offset+=n
 assert out;pin['polygons']=len(out);return out,pin
def box(ax,b,color):
 x0,x1,z0,z1,y0,y1=b;q=np.array([[x,z,y]for x in[x0,x1]for z in[z0,z1]for y in[y0,y1]])
 faces=[q[i]for i in[[0,1,3,2],[4,5,7,6],[0,1,5,4],[2,3,7,6],[0,2,6,4],[1,3,7,5]]]
 rgb=np.array(to_rgb(color));colors=[tuple(np.clip(rgb*x,0,1))for x in[.65,.86,.76,.94,.7,1.0]]
 ax.add_collection3d(Poly3DCollection(faces,facecolors=colors,edgecolors=(.6,.8,.9,.15),linewidths=.35))
def render3d(polys=None,grid=False,azim=-52):
 fig=plt.figure(figsize=(9,8.4),dpi=150,facecolor=BG);ax=fig.add_axes([-.04,-.05,1.08,1.1],projection='3d',facecolor=BG)
 if polys:
  q=[p[:,[0,2,1]]*1000 for p in polys];colors=[]
  light=np.array([.25,-.5,.82]);light/=np.linalg.norm(light)
  for p in q:
   norm=np.cross(p[1]-p[0],p[2]-p[0]);norm/=max(np.linalg.norm(norm),1e-30)
   shade=.58+.42*abs(np.dot(norm,light));rgb=np.array(to_rgb(TEAL))*shade
   colors.append((*rgb,.45 if p[:,2].mean()<10 else .8))
  ax.add_collection3d(Poly3DCollection(q,facecolors=colors,edgecolors='none',zsort='average'))
 else:
  q=np.array([[-10,-10,8],[10,-10,8],[10,10,8],[-10,10,8]])
  ax.add_collection3d(Poly3DCollection([q],facecolors=[(.12,.72,.65,.23)],edgecolors=TEAL,linewidths=.65))
 box(ax,(-4,4,-4,4,0,1.5),COPPER);box(ax,(-6,6,-6,6,19,20),BLUE)
 # Wireframe domain, exact 20 mm cube. Mesh display samples every eighth grid line.
 lines=[]
 for y in[0,20]:
  lines += [np.array([[-10,-10,y],[10,-10,y],[10,10,y],[-10,10,y],[-10,-10,y]])]
 for x in[-10,10]:
  for z in[-10,10]:lines.append(np.array([[x,z,0],[x,z,20]]))
 ax.add_collection3d(Line3DCollection(lines,colors=(.45,.65,.8,.4),linewidths=.65))
 if grid:
  lines=[]
  for v in np.linspace(-10,10,11):
   for y in[0,8,20]:lines.extend([np.array([[-10,v,y],[10,v,y]]),np.array([[v,-10,y],[v,10,y]])])
   lines.extend([np.array([[v,10,0],[v,10,20]]),np.array([[-10,v,0],[-10,v,20]])])
  ax.add_collection3d(Line3DCollection(lines,colors=(.33,.65,.85,.2),linewidths=.45))
 ax.set(xlim=(-11,11),ylim=(-11,11),zlim=(-1,21));ax.set_box_aspect((1,1,1));ax.view_init(elev=22,azim=azim);ax.set_axis_off()
 b=io.BytesIO();fig.savefig(b,format='png',facecolor=BG);plt.close(fig);b.seek(0);return Image.open(b).convert('RGB')
def slice_image(sl,poly,flow=False):
 fig=plt.figure(figsize=(7.2,7.5),dpi=150,facecolor=BG);ax=fig.add_axes([.09,.08,.74,.87],facecolor=BG)
 X=sl[:,:,0]*1000;Y=sl[:,:,1]*1000;S=sl[:,:,8]>.5;speed=np.linalg.norm(sl[:,:,5:8],axis=2)
 data=np.ma.masked_where(S,speed)if flow else sl[:,:,3]-273.15
 im=ax.imshow(data,extent=(-10,10,0,20),origin='lower',cmap='viridis'if flow else'inferno',vmin=0 if flow else 45,vmax=.6 if flow else 63,interpolation='bilinear',aspect='equal')
 z=sl[0,0,2];segs=[]
 for p in poly:
  pts=[]
  for u,v in zip(p,np.roll(p,-1,axis=0)):
   d1=u[2]-z;d2=v[2]-z
   if abs(d1)<1e-10:pts.append(u[:2])
   if d1*d2<0:pts.append((u+d1/(d1-d2)*(v-u))[:2])
  if len(pts)>=2:
   pp=np.array(pts);d=((pp[:,None]-pp[None,:])**2).sum(2);i,j=np.unravel_index(d.argmax(),d.shape)
   if d[i,j]>1e-20:segs.append(pp[[i,j]]*1000)
 ax.add_collection(LineCollection(segs,colors=TEAL,linewidths=.65))
 if flow:
  ix=(slice(2,None,5),slice(2,None,5));mask=S[ix]|(speed[ix]<.003)
  ax.quiver(X[ix],Y[ix],np.ma.masked_where(mask,sl[:,:,5][ix]),np.ma.masked_where(mask,sl[:,:,6][ix]),color='white',angles='xy',scale_units='xy',scale=.16,width=.0035)
 ax.set_xticks([-10,0,10]);ax.set_yticks([0,10,20]);ax.tick_params(colors=MUTED,labelsize=12)
 for s in ax.spines.values():s.set_edgecolor('#315061')
 ax.set_xlabel('x [mm]',color=MUTED,fontsize=13);ax.set_ylabel('y [mm]',color=MUTED,fontsize=13)
 cb=fig.colorbar(im,cax=fig.add_axes([.87,.16,.03,.7]));cb.ax.tick_params(colors=FG,labelsize=12);cb.set_label('Speed [m/s]'if flow else'Temperature [°C]',color=FG,fontsize=13)
 b=io.BytesIO();fig.savefig(b,format='png',facecolor=BG);plt.close(fig);b.seek(0);return Image.open(b).convert('RGB')
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--execution',type=Path,required=True);ap.add_argument('--limit',type=float,default=.02);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--ffmpeg',required=True);ap.add_argument('--audio-dir',type=Path);ap.add_argument('--max-frames',type=int,default=12);ap.add_argument('--fps',type=int,default=24);ap.add_argument('--forms',nargs='+',default=['landscape','portrait']);a=ap.parse_args()
 out=a.output;out.mkdir(parents=True,exist_ok=False);assets=out/'assets';assets.mkdir();work=a.run
 eraw,epin=stable(a.execution);ex=json.loads(eraw);assert ex['grid']==[80,80,80]
 hraw,hpin=stable(work/'history.csv');hist=[{k:float(v)for k,v in r.items()}for r in csv.DictReader(hraw.decode().splitlines())];assert hist
 case=work/'ensight/VaporizingDrop/plic.case';raw,cpin=stable(case);times=[float(t)for t in raw.decode().split('time values:')[1].split()]
 eligible=[]
 for i,t in enumerate(times):
  p=case.parent/'plic'/f'plic.{i+1:06d}'
  if t>a.limit+1e-8 or not p.exists()or time.time()-p.stat().st_mtime<10:continue
  h=min(hist,key=lambda r:abs(r['time_s']-t))
  if abs(h['time_s']-t)<1e-7:eligible.append((t,p,h))
 assert eligible
 if len(eligible)>a.max_frames:eligible=[eligible[i]for i in sorted(set(np.linspace(0,len(eligible)-1,a.max_frames).round().astype(int)))]
 native=[];pins=[];finalpoly=None
 for i,(t,p,h)in enumerate(eligible):
  poly,pin=plic(p);im=render3d(poly);im.save(assets/f'PLIC_{i:02d}.png');native.append(im);pins.append({'physicalTime_s':t,'step':int(h['step']),'nativePLIC':pin});finalpoly=poly
  print(json.dumps({'renderedPLIC':i,'time_s':t}),flush=True)
 finaltime=eligible[-1][0];step=int(eligible[-1][2]['step']);fields=sorted(work.glob(f'field_step{step:08d}_rank????.csv'));assert len(fields)==8,'Need all 8 field packets at final native PLIC time'
 arrays=[];fpins=[]
 for p in fields:
  raw,pin=stable(p);arrays.append(np.loadtxt(io.BytesIO(raw),delimiter=',',skiprows=1));fpins.append(pin)
 vals=np.concatenate(arrays);assert vals.shape==(512000,12)and np.isfinite(vals).all();vals=vals[np.lexsort((vals[:,0],vals[:,1],vals[:,2]))].reshape(80,80,80,12)
 assert np.all((vals[:,:,:,4]>=0)&(vals[:,:,:,4]<=1));gv=float(np.sum((1-vals[:,:,:,4])*(vals[:,:,:,8]<.5)))*.00025**3
 assert abs(gv-eligible[-1][2]['gas_volume_m3'])<1e-13
 geometry=render3d();mesh=render3d(grid=True);temp=slice_image(vals[40],finalpoly);velocity=slice_image(vals[40],finalpoly,flow=True)
 for name,im in [('SETUP',geometry),('GRID',mesh),('TEMPERATURE',temp),('VELOCITY',velocity)]:im.save(assets/(name+'.png'))
 fontdir=Path(matplotlib.get_data_path())/'fonts/ttf'
 def font(n,bold=False):return ImageFont.truetype(str(fontdir/('DejaVuSans-Bold.ttf'if bold else'DejaVuSans.ttf')),n)
 outputs=[]
 for form in a.forms:
  portrait=form=='portrait';W,H=(1080,1920)if portrait else(1920,1080);scale=W/1080 if portrait else H/1080
  def f(n,bold=False):return font(round(n*scale),bold)
  bg=Image.new('RGB',(W,H),BG);d=ImageDraw.Draw(bg)
  # Restrained geometric panel lines: same styling in both compositions.
  d.line((round(W*.06),round(H*.08),round(W*.94),round(H*.08)),fill='#1c3c50',width=2)
  titles=['IMMERSION\nCOOLING'if portrait else'IMMERSION COOLING','HEAT IN.\nHEAT OUT.'if portrait else'HEAT IN. HEAT OUT.','A CARTESIAN\nFOUNDATION'if portrait else'A CARTESIAN FOUNDATION','FOLLOW THE\nINTERFACE'if portrait else'FOLLOW THE INTERFACE','HEAT MEETS\nFLOW'if portrait else'HEAT MEETS FLOW','EARLY 3D\nEVOLUTION'if portrait else'EARLY 3D EVOLUTION']
  def paste_fit(canvas,im,rect):
   x,y,w,h=map(round,rect);q=im.copy();q.thumbnail((w,h),Image.Resampling.LANCZOS);canvas.paste(q,(x+(w-q.width)//2,y+(h-q.height)//2))
  def frame(t):
   scene=next(i for i in range(6)if t<CUTS[i+1]);u=(t-CUTS[scene])/(CUTS[scene+1]-CUTS[scene]);canvas=bg.copy();draw=ImageDraw.Draw(canvas)
   draw.text((W*.06,H*.035),'PHASE CHANGE  /  HEAT TRANSFER',font=f(20,True),fill=TEAL)
   titleY=H*.11 if portrait else H*.12;draw.multiline_text((W*.06,titleY),titles[scene],font=f(76 if portrait else 64,True),spacing=8,fill=FG)
   rect=(W*.02,H*.28,W*.96,H*.49)if portrait else(W*.36,H*.24,W*.63,H*.65)
   im=geometry;actual=None
   if scene==2:im=mesh
   if scene==3:
    index=min(len(native)-1,int(u*len(native)));im=native[index];actual=eligible[index][0]
   if scene==4:im=temp if u<.5 else velocity;actual=finaltime
   if scene==5:im=native[-1];actual=finaltime
   if scene==0:
    factor=1.0+.035*u;size=(round(im.width*factor),round(im.height*factor));q=im.resize(size,Image.Resampling.LANCZOS);x=(q.width-im.width)//2;y=(q.height-im.height)//2;im=q.crop((x,y,x+im.width,y+im.height))
   if scene==2:im=Image.blend(geometry,mesh,min(1.,u*1.5))
   paste_fit(canvas,im,rect)
   textY=H*.80 if portrait else H*.43;textX=W*.065
   rows=[]
   if scene==0:rows=[('3D','Liquid + vapor'),('CHT','Solid + fluid')]
   if scene==1:rows=[('3.79 W','Heated chip'),('45 °C','Cold plate'),('20 mm','Domain width')]
   if scene==2:rows=[('512,000','Cartesian cells'),('0.25 mm','Cell spacing')]
   if scene==3:rows=[('PLIC','Native XYZ interface')]
   if scene==4:rows=[('VELOCITY'if u>=.5 else'TEMPERATURE','Calculated fields')]
   if scene==5:rows=[(f'{1000*finaltime:.0f} ms','Accepted evolution'),('8 MPI','CPU calculation')]
   if portrait:
    gap=W*.88/max(len(rows),1)
    for j,(big,small)in enumerate(rows):
     x=textX+j*gap;draw.text((x,textY),big,font=f(35,True),fill=COPPER if scene==1 else TEAL);draw.text((x,textY+52),small,font=f(20),fill=MUTED)
   else:
    for j,(big,small)in enumerate(rows):
     y=textY+j*H*.15;draw.text((textX,y),big,font=f(48,True),fill=COPPER if scene==1 else TEAL);draw.text((textX,y+67),small,font=f(23),fill=MUTED)
   if actual is not None:
    txt=f't = {1000*actual:.3f} ms';tw=draw.textlength(txt,font=f(24,True));draw.rounded_rectangle((W*.94-tw-24,H*.091,W*.94+6,H*.091+50),radius=12,fill='#102c3b');draw.text((W*.94-tw-12,H*.091+8),txt,font=f(24,True),fill=FG)
   # Result qualification is short; full scientific scope accompanies the video.
   foot='EARLY-TIME NUMERICAL STUDY'if scene!=3 else'ACTUAL PLIC SNAPSHOTS  •  SLOW-MOTION PLAYBACK'
   draw.text((W*.06,H*.925),foot,font=f(19),fill=MUTED)
   for i in range(6):
    x0=W*(.06+i*.148);x1=x0+W*.128;y=H*.972;draw.line((x0,y,x1,y),fill=TEAL if i<=scene else'#244154',width=5)
   return canvas
  poster=frame(21);poster.save(out/f'POSTER_{form.upper()}.png')
  sheet=Image.new('RGB',(3*480,2*round(480*H/W)),BG)
  for j,t in enumerate([1,4.5,8.5,12.5,17,22]):
   im=frame(t);im.thumbnail((480,round(480*H/W)),Image.Resampling.LANCZOS);sheet.paste(im,((j%3)*480,(j//3)*round(480*H/W)))
  sheet.save(out/f'CONTACT_{form.upper()}.png')
  target=out/f'IMMERSION_COOLING_{form.upper()}_DRAFT01.mp4';silent=out/f'_{form}_silent.mp4'if a.audio_dir else target
  command=[a.ffmpeg,'-hide_banner','-loglevel','error','-n','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r',str(a.fps),'-i','pipe:0','-an','-c:v','libx264','-threads','2','-preset','veryfast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(silent)]
  with(out/f'encoder_{form}.log').open('wb')as log:
   p=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=subprocess.DEVNULL,stderr=log)
   try:
    for i in range(24*a.fps):p.stdin.write(frame(i/a.fps).tobytes())
    p.stdin.close();assert p.wait(timeout=60)==0,'Encoder failed'
   finally:
    if p.poll()is None:p.kill();p.wait()
  if a.audio_dir:
   args=[];filters=[];voice=[]
   for i in range(6):
    path=a.audio_dir/f'VOICE_{i:02d}.wav'
    with wave.open(str(path),'rb')as z:duration=z.getnframes()/z.getframerate()
    assert duration<CUTS[i+1]-CUTS[i]-.15,(i,duration,'Re-synthesize shorter or faster natively; do not time-stretch')
    args+=['-i',str(path)];filters.append(f'[{i+1}:a]aresample=48000,adelay={round((CUTS[i]+.08)*1000)}[v{i}]');voice.append({'scene':i,'duration_s':duration,'postTimeStretch':False})
   filters.append(''.join(f'[v{i}]'for i in range(6))+'amix=inputs=6:normalize=0,apad,atrim=0:24,loudnorm=I=-16:TP=-1.5:LRA=7[a]')
   subprocess.run([a.ffmpeg,'-hide_banner','-loglevel','error','-n','-i',str(silent),*args,'-filter_complex',';'.join(filters),'-map','0:v','-map','[a]','-c:v','copy','-c:a','aac','-b:a','160k','-t','24','-movflags','+faststart',str(target)],check=True,timeout=60)
   silent.unlink() # only this script's just-created temporary encode
  # Full decode catches corrupt/truncated MP4s, independently of encode exit.
  subprocess.run([a.ffmpeg,'-hide_banner','-loglevel','error','-i',str(target),'-f','null','-'],check=True,timeout=90)
  outputs.append({'name':target.name,'width':W,'height':H,'fps':a.fps,'frames':24*a.fps,'duration_s':24,'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'fullDecodePassed':True,'audio':bool(a.audio_dir)})
  print(json.dumps({'encoded':target.name,'bytes':target.stat().st_size}),flush=True)
 provenance={'sourceRunLeaf':work.parent.name+'/'+work.name,'executionSHA256':epin['sha256'],'binarySHA256':ex['binarySHA256'],'history':hpin,'nativeCase':cpin,'frames':pins,'fieldsAtFinalTime':fpins,'fieldShape':[80,80,80,12],'centralSliceZ_m':float(vals[40,0,0,2]),'temperatureRange_C':[45,63],'speedRange_m_s':[0,.6],'displayInterpolation':'Bilinear image pixels only; native PLIC polygons are not smoothed into new geometry. No intermediate CFD states generated.','initialColdFilmIncluded':True,'reservoirVaporPortsIncluded':True,'completeCoolingCycleQualified':False,'detachedBubbleRiseClaimed':False,'voiceTimeStretch':False,'outputs':outputs}
 (out/'PROVENANCE.json').write_text(json.dumps(provenance,indent=2)+'\n');print(json.dumps({'output':str(out),'nativeFrames':len(native),'lastPhysicalTime_s':finaltime}))
if __name__=='__main__':main()
