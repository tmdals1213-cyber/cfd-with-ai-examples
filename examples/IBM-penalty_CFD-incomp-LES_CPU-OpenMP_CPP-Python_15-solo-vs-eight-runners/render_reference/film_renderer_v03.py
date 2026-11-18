"""English intro with growing 3D grid, actual paired flow, pressure and signed loads."""
from film_renderer_v02 import *
from functools import lru_cache
from matplotlib import colormaps
V3=R/'v03';SECTION=json.loads((V3/'assets/sections/SECTION_ANIMATION_INDEX.json').read_text())
SERIES=np.load(WS/'work/runner_pressure_load_surface95_20261113/PRESSURE_SERIES_SIGNED.npz',allow_pickle=False)
SUMMARY=json.loads((WS/'work/runner_pressure_load_surface95_20261113/PRESSURE_COMPARISON_SIGNED_SUMMARY.json').read_text())
STAGES=[(0,4.5),(4.5,11),(11,18),(18,26),(26,34),(34,40),(40,46),(46,51)]
NEW_TITLES=['','','','','Same view.\nSame scales.','Pressure changes.\nForce has a direction.','From surface colour\nto pressure load.','Lower pressure load.\nA model estimate.']
LUT=(colormaps['RdBu_r'](np.linspace(0,1,256))[:,:3]*255).astype('u1')
@lru_cache(maxsize=5)
def npz_cache(path):
 with np.load(path,allow_pickle=False)as z:return {k:z[k].copy()for k in z.files}
def section(case,time_s):
 row=min(SECTION['cases'][case],key=lambda r:abs(r['time_s']-time_s));return row,npz_cache(row['archive']['path'])
@lru_cache(maxsize=2)
def flow_index(case):
 root=V3/'assets/smooth_flow'/case;rows=json.loads((V3/'assets/smooth_flow/SMOOTH_FLOW_INDEX.json').read_text())['cases'][case];return root,rows
def actual_flow(size,case,time_s):
 root,rows=flow_index(case);row=min(rows,key=lambda r:abs(r['time_s']-time_s));f=npz_cache(str(root/row['file']));w,h=size;center=(.65,.83,0);scale=min(w/7.4,(h-80)/3.3)
 im=Image.new('RGBA',size);d=ImageDraw.Draw(im);ground(d,size,center,scale,wide=True)
 mesh,depth=bodies(size,0,OFFSETS[:1]if case=='solo'else OFFSETS,center,scale,body=f['body'],normals=f['normals']);im.alpha_composite(mesh)
 palette=[(0.,.3,1.),(0.,1.,.9),(.8,1.,.1),(1.,.12,0.)]
 for j in range(4):
  vv=f['surface%d'%j]
  if not len(vv):continue
  projected_v=projected(vv,center,scale,size);fac=f['surfaceFaces%d'%j];shade=.65+.35*np.abs(f['surfaceNormals%d'%j]@np.array([.45,.8,.38]));color=np.asarray(palette[j])[None,:]*shade[:,None];layer,z=draw_mesh(size,projected_v,fac,color)
  rgba=np.array(layer);rgba[...,3]=np.where((rgba[...,3]>0)&(z>depth),[25,35,46,64][j],0).astype('u1');im.alpha_composite(Image.fromarray(rgba))
 d=ImageDraw.Draw(im)
 for x,y,z in projected(f['particles'],center,scale,size):
  ix,iy=int(x),int(y)
  if 0<=ix<w and 0<=iy<h and z>depth[iy,ix]:d.ellipse((x-1,y-1,x+1,y+1),fill=(160,220,235,190))
 label(d,(24,h-74),'t = %.3f s'%row['time_s'],25);large_legend(d,size,40)
 return im
def lattice(size,t,local):
 w,h=size;center=(0,.87,0);scale=min(w/2.9,h/2.7);az=40+3*math.sin(local*.3);im=Image.new('RGBA',size)
 layer=Image.new('RGBA',size);d=ImageDraw.Draw(layer);growth=smooth(local/4.2);bounds=[(-1.05,1.05),(.02,2.05),(-.77,.77)]
 for axis in range(3):
  others=[b for b in range(3)if b!=axis];coords=[np.arange(bounds[b][0],bounds[b][1]+.001,.07)for b in others]
  for i,q in enumerate(coords[0]):
   for j,r in enumerate(coords[1]):
    coarse=i%3==0 and j%3==0
    if not coarse and ((i*29+j*17+axis*11)%101)/100>growth:continue
    p=np.zeros((2,3));p[:,axis]=bounds[axis];p[:,others[0]]=q;p[:,others[1]]=r
    p3(d,p,center,scale,size,az,10,(55,171,194,66 if coarse else 32),1)
 im.alpha_composite(layer)
 body,_=hero(size,t*.48,OFFSETS[:1],center,scale,az,10);im.alpha_composite(body)
 layer.putalpha(layer.getchannel('A').point(lambda a:a//3));im.alpha_composite(layer)
 d=ImageDraw.Draw(im);label(d,(24,h-50),'SCHEMATIC',27,MUTED)
 if local>4.5:
  panel=Image.new('RGBA',size);pd=ImageDraw.Draw(panel);pd.rounded_rectangle((12,h-320 if h>w else h-235,w-12,h-75),radius=12,fill=(5,10,15,238))
  equation(panel,'CONTINUITY',(25,h-320 if h>w else h-235,w-25,h-220 if h>w else h-165));equation(panel,'MOMENTUM',(25,h-220 if h>w else h-165,w-25,h-92));im.alpha_composite(panel)
 return im
@lru_cache(maxsize=5)
def scalar_display(size,case,index,top):
 row=SECTION['cases'][case][index];a=npz_cache(row['archive']['path']);j=str(row['localIndex']);p=a[('topP'if top else'sideP')+j];fluid=a[('topFluid'if top else'sideFluid')+j]
 w,h=size;xx=np.linspace(-1.65,4.4,w);yy=np.linspace(2.55 if top else 2.25,-2.55 if top else 0,h);cent=[a['x'],a['z'if top else'y']]
 ix=np.clip(np.searchsorted(cent[0],xx)-1,0,len(cent[0])-2);iy=np.clip(np.searchsorted(cent[1],yy)-1,0,len(cent[1])-2);sx=np.clip((xx-cent[0][ix])/(cent[0][ix+1]-cent[0][ix]),0,1);sy=np.clip((yy-cent[1][iy])/(cent[1][iy+1]-cent[1][iy]),0,1)
 value=np.zeros((h,w));weight=np.zeros((h,w))
 for a0 in(0,1):
  for a1 in(0,1):
   wt=np.where(a1,sy,1-sy)[:,None]*np.where(a0,sx,1-sx)[None,:];valid=fluid[ix[None,:]+a0,iy[:,None]+a1];weight+=wt*valid;value+=wt*valid*p[ix[None,:]+a0,iy[:,None]+a1]
 value/=np.maximum(weight,1e-12);rgb=LUT[np.clip(((value+50)/100*255).astype(int),0,255)].astype('f8');alpha=np.clip(weight/.5,0,1);rgb=rgb*alpha[...,None]+np.array([77,92,104])*(1-alpha[...,None]);return Image.fromarray(rgb.astype('u1')).convert('RGBA')
def pressure_panel(size,case,time_s,top=True):
 row,a=section(case,time_s);j=str(row['localIndex']);w,h=size;im=Image.new('RGBA',size);d=ImageDraw.Draw(im)
 p=a[('topP'if top else'sideP')+j];u=a[('topU'if top else'sideU')+j];fluid=a[('topFluid'if top else'sideFluid')+j];vaxis=2 if top else 1
 xf=a['faceX'];yf=a['faceZ'if top else'faceY'];xc=a['x'];yc=a['z'if top else'y'];xmin,xmax=-1.65,4.4;ymin,ymax=(-2.55,2.55)if top else(0,2.25)
 topmargin=64;bottom=102;scale=min((w-25)/(xmax-xmin),(h-topmargin-bottom)/(ymax-ymin));left=(w-scale*(xmax-xmin))/2;upper=topmargin+((h-topmargin-bottom)-scale*(ymax-ymin))/2
 def screen(x,y):return(left+(x-xmin)*scale,upper+(ymax-y)*scale)
 index=SECTION['cases'][case].index(row);raster=scalar_display((max(1,round(scale*(xmax-xmin))),max(1,round(scale*(ymax-ymin)))),case,index,top);im.alpha_composite(raster,(round(left),round(upper)));d=ImageDraw.Draw(im)
 for ix in range(0,len(xc),3):
  for iy in range(0,len(yc),3):
   x,y=xc[ix],yc[iy]
   if not(xmin+.1<x<xmax-.1 and ymin+.1<y<ymax-.1 and fluid[ix,iy]):continue
   start=screen(x,y);end=screen(x+u[ix,iy,0]/20,y+u[ix,iy,vaxis]/20);arrow(d,start,end,(13,131,135),2,6)
 label(d,(18,12),'SOLO'if case=='solo'else'PACED',31,GOLD if case=='solo'else CYAN)
 txt(d,(w-18,15),'t = %.3f s'%row['time_s'],23,MUTED,False,'ra')
 # Axial force direction is separate from velocity: Fx = -Dp.
 fy=h-74;force=row['forceOnRunnerX_N'];sign=1 if force>=0 else-1;cx=w/2;length=abs(force)*8
 if length>=2:arrow(d,(cx-sign*length/2,fy),(cx+sign*length/2,fy),GOLD,5,min(15,length*.4))
 label(d,(w/2,h-42),'PRESSURE FORCE Fx  %+.2f N'%force,25,WHITE,'mm')
 txt(d,(w/2,h-12),'Fx = −Dp   /   running +X →',20,MUTED,False,'mm')
 return im
def pressure_surface(size,case,time_s):
 row,a=section(case,time_s);j=str(row['localIndex']);w,h=size;body=a['body'+j];normal=a['normal'+j];value=a['pnx'+j];center=(0,.86,0)
 color=LUT[np.clip(((value+50)/100*255).astype(int),0,255)]/255;shade=.48+.52*np.maximum(0,normal@np.array([.45,.8,.38]));color=np.clip(color*shade[:,None],0,1)
 v=projected(body+OFFSETS[0],center,min(w/1.35,(h-100)/2.05),size,40,10);im,_=draw_mesh(size,v,TRI,color);d=ImageDraw.Draw(im)
 label(d,(w/2,8),'SOLO'if case=='solo'else'PACED',32,GOLD if case=='solo'else CYAN,'ma');label(d,(w/2,h-24),'%+.2f N'%row['load_N'],38,WHITE,'mm')
 return im
def history_plot(size,time_s):
 w,h=size;im=Image.new('RGBA',size);d=ImageDraw.Draw(im);margin=[80,65,30,65];x0,y0,x1,y1=margin[0],margin[1],w-margin[2],h-margin[3]
 tx=SERIES['commonTime'];solo=SERIES['soloCommon'];paced=SERIES['pacedCommon'];lo=math.floor(min(solo.min(),paced.min())/5)*5;hi=math.ceil(max(solo.max(),paced.max())/5)*5
 def xy(t,f):return(x0+(t-.625)/1.875*(x1-x0),y1-(f-lo)/(hi-lo)*(y1-y0))
 for yy in np.linspace(lo,hi,5):
  y=xy(.625,yy)[1];d.line((x0,y,x1,y),fill=EDGE,width=1);txt(d,(x0-12,y),'%g'%yy,23,MUTED,False,'rm')
 d.line((*xy(.625,0),*xy(2.5,0)),fill=(80,110,125),width=2)
 for x in(.625,1.25,1.875,2.5):
  xx=xy(x,lo)[0];txt(d,(xx,y1+15),'%.2f'%x,23,MUTED,False,'ma')
 for values,color in((solo,GOLD),(paced,CYAN)):
  sel=tx<=time_s;pts=[xy(t,v)for t,v in zip(tx[sel],values[sel])]
  if len(pts)>1:d.line(pts,fill=color,width=4)
  v=float(np.interp(time_s,tx,values));x,y=xy(time_s,v);d.ellipse((x-6,y-6,x+6,y+6),fill=color)
 x=xy(time_s,lo)[0];d.line((x,y0,x,y1),fill=(150,174,184),width=1)
 label(d,(x0,8),'PRESSURE LOAD  Dp [N]',27);txt(d,(x1,10),'time [s]',24,MUTED,False,'ra')
 label(d,(w/2,h-12),'SOLO',23,GOLD,'mm');label(d,(w/2+110,h-12),'PACED',23,CYAN,'mm')
 return im
def pressure_scale(d,size,surface=False):
 w,h=size;bw=260;x=w-bw-28;y=h-40
 for j in range(bw):d.line((x+j,y,x+j,y+12),fill=tuple(LUT[int(j/(bw-1)*255)]))
 label(d,(x,y-39),'p nₓ  [Pa]'if surface else'p − p_ref  [Pa]',25,MUTED)
 txt(d,(x,y+15),'−50',21,MUTED);txt(d,(x+bw,y+15),'+50',21,MUTED,False,'ra')
def frame_v03(t,form='portrait',scale_output=1):
 size=(1080,1920)if form=='portrait'else(1920,1080);w,h=size;portrait=h>w;m=64 if portrait else 84
 idx=next((i for i,(a,b)in enumerate(STAGES)if a<=t<b),7);start,end=STAGES[idx];local=t-start;a=local/(end-start)
 box=(44,415,w-44,1720)if portrait else(65,300,w-65,980);aw,ah=box[2]-box[0],box[3]-box[1]
 if idx<4:
  im=frame_v02(t,form).convert('RGBA');background=base(size);rect=(m+44,35 if not portrait else 55,m+400,95);im.paste(background.crop(rect),rect)
  if idx==2:
   im.paste(background.crop(box),box);im.alpha_composite(lattice((aw,ah),t,local),(box[0],box[1]))
 else:
  im=base(size);d=ImageDraw.Draw(im);d.line((m,74 if portrait else 54,m+35,74 if portrait else 54),fill=CYAN,width=3)
  textlines(d,(m,180 if portrait else 110),NEW_TITLES[idx],78 if portrait else 76,w-2*m,WHITE,True,1.08)
  art=Image.new('RGBA',(aw,ah));ad=ImageDraw.Draw(art);time_s=.625+a*1.875
  if idx in(4,5):
   gap=28;boxes=[(0,0,aw,(ah-gap)//2),(0,(ah+gap)//2,aw,ah)]if portrait else[(0,0,(aw-gap)//2,ah),((aw+gap)//2,0,aw,ah)]
   for j,b in enumerate(boxes):
    pw,ph=b[2]-b[0],b[3]-b[1];case='solo'if j==0 else'group8'
    panel=actual_flow((pw,ph),case,time_s)if idx==4 else pressure_panel((pw,ph-60),case,time_s,top=local<3)
    art.alpha_composite(panel,(b[0],b[1]));label(ad,(b[0]+22,b[1]+20),'SOLO'if j==0 else'PACED',31,GOLD if j==0 else CYAN)if idx==4 else None
   if idx==5:
    label(ad,(26,ah-36),'AIR VELOCITY  /  PRESSURE FORCE',23,CYAN);pressure_scale(ad,(aw,ah))
  elif idx==6:
   if portrait:sw,sh=(aw-24)//2,ah*54//100;chart=(0,sh+90,aw,ah-25)
   else:sw,sh=(aw//2-30)//2,ah-85;chart=(aw//2,55,aw,ah-25)
   for j,case in enumerate(('solo','group8')):art.alpha_composite(pressure_surface((sw,sh),case,time_s),(j*(sw+24),0))
   art.alpha_composite(history_plot((chart[2]-chart[0],chart[3]-chart[1]),time_s),(chart[0],chart[1]));pressure_scale(ad,(2*sw+24,sh+67),True)
  else:
   if portrait:
    label(ad,(aw/2,70),'MEAN PRESSURE LOAD',35,MUTED,'mm');label(ad,(aw*.27,225),'SOLO',34,GOLD,'mm');label(ad,(aw*.73,225),'PACED',34,CYAN,'mm');label(ad,(aw*.27,335),'8.37 N',76,WHITE,'mm');label(ad,(aw*.73,335),'−0.25 N',76,WHITE,'mm')
    label(ad,(aw/2,550),'DIFFERENCE',35,MUTED,'mm');label(ad,(aw/2,690),'8.6 N',132,GOLD,'mm');label(ad,(aw/2,930),'ONE RUNNER.',39,WHITE,'mm');label(ad,(aw/2,1005),'DIFFERENT FLOW.',39,CYAN,'mm');txt(ad,(aw/2,ah-100),'Mean: t = 1.25–2.50 s  /  two gait cycles',27,MUTED,False,'mm')
   else:
    label(ad,(aw*.22,100),'SOLO',31,GOLD,'mm');label(ad,(aw*.22,185),'8.37 N',76,WHITE,'mm');label(ad,(aw*.22,315),'PACED',31,CYAN,'mm');label(ad,(aw*.22,400),'−0.25 N',76,WHITE,'mm');label(ad,(aw*.67,110),'MEAN DIFFERENCE',34,MUTED,'mm');label(ad,(aw*.67,245),'8.6 N',132,GOLD,'mm');label(ad,(aw*.67,410),'ONE RUNNER. DIFFERENT FLOW.',32,CYAN,'mm');txt(ad,(aw/2,ah-65),'Mean: t = 1.25–2.50 s  /  two gait cycles',27,MUTED,False,'mm')
   label(ad,(aw/2,ah-20),'PRESSURE-ONLY ESTIMATE  ≠  TOTAL DRAG',28,MUTED,'mm')
  im.alpha_composite(art,(box[0],box[1]))
 bg=base(size);rect=(w-m-180,35,w-m+5,95);im.paste(bg.crop(rect),rect);d=ImageDraw.Draw(im);txt(d,(w-m,61 if portrait else 41),'%02d / 08'%(idx+1),22,CYAN,False,'ra');py=1830 if portrait else 1038
 # Replace inherited six-chapter progress with the eight-chapter revision.
 bg=base(size);rect=(m,py-50,w-m,py+10);im.paste(bg.crop(rect),rect);d=ImageDraw.Draw(im)
 d.line((m,py,w-m,py),fill=EDGE,width=2);d.line((m,py,m+(w-2*m)*(idx+a)/8,py),fill=CYAN,width=3)
 if idx==4:txt(d,(m,py-44),'CFD  /  SMOOTHED DISPLAY  /  PLAYBACK 0.23x',25,MUTED)
 elif idx in(5,6):txt(d,(m,py-44),'SMOOTHED DISPLAY  /  PRESSURE-ONLY MODEL',25,MUTED)
 fade=min(smooth(local/.22),smooth((end-t)/.18))
 if idx>=4 and fade<1:im=Image.blend(base(size),im,fade)
 if scale_output!=1:im=im.resize((int(w*scale_output),int(h*scale_output)),Image.Resampling.LANCZOS)
 return im.convert('RGB')
