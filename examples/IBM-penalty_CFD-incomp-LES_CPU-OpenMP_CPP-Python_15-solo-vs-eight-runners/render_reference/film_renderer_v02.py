"""Large titles, physical quantities, brief equations and paired surface views."""
from film_renderer import *
V=R/'v02';ASSETS=V/'assets'
TIMES=[(0,4.5),(4.5,11),(11,18),(18,26),(26,37),(37,44)]
TITLES=["Can a team change\none runner's drag?",'Same speed.\nDifferent shelter.','The body moves.\nThe grid stays.','A moving surface.\nA force on the fluid.','Same view.\nSame scales.','Compare the\nsurface loads.']
CHAPTERS=['THE QUESTION','THE PHYSICS','THE CARTESIAN GRID','THE IMMERSED BOUNDARY','THE FLOW FIELD','THE COMPARISON']
PC=np.load(V/'assets_pressure50/PRESSURE_COMPARISON.npz',allow_pickle=False)
EQUATIONS={name:Image.open(ASSETS/(name+'.png')).convert('RGBA')for name in('CONTINUITY','MOMENTUM','IBM')}

def label(d,xy,value,size=32,color=WHITE,anchor=None):txt(d,xy,value,size,color,True,anchor)
def equation(im,name,box,alpha=1):
 image=EQUATIONS[name].copy();x,y,X,Y=box;factor=min((X-x)/image.width,(Y-y)/image.height);image=image.resize((int(image.width*factor),int(image.height*factor)),Image.Resampling.LANCZOS)
 if alpha<1:image.putalpha(image.getchannel('A').point(lambda v:int(v*alpha)))
 im.alpha_composite(image,(int((x+X-image.width)/2),int((y+Y-image.height)/2)))
def large_legend(d,size,limit,kind='vorticity'):
 w,h=size;bw=min(w-60,260);x=w-bw-28;y=h-40
 for j in range(bw):
  u=j/(bw-1)
  if kind=='vorticity':c=np.clip(1.5-np.abs(4*u-np.array([3,2,1])),0,1)
  else:
   from matplotlib import colormaps
   c=np.array(colormaps['RdBu_r'](u)[:3])
  d.line((x+j,y,x+j,y+12),fill=tuple((255*c).astype(int)))
 label(d,(x,y-39),'|ω|  [s⁻¹]'if kind=='vorticity'else'p nₓ  [Pa]',25,MUTED)
 txt(d,(x,y+15),'0'if kind=='vorticity'else'-%g'%limit,21,MUTED);txt(d,(x+bw,y+15),'%g'%limit,21,MUTED,False,'ra')

def clean_flow(size,a):
 cache=json.loads((FLOWROOT/'FLOW_CACHE_PASS.json').read_text());frames=cache['frames'];meta=frames[int(round(max(0,min(1,a))*(len(frames)-1)))];f=flow_frame(meta['sourceIndex'])
 w,h=size;center=(.65,.83,0);scale=min(w/7.4,(h-70)/3.3);im=Image.new('RGBA',size);d=ImageDraw.Draw(im);ground(d,size,center,scale,wide=True)
 mesh,depth=bodies(size,0,OFFSETS,center,scale,body=f['body'],normals=f['normals']);im.alpha_composite(mesh)
 palette=[(0.,.3,1.),(0.,1.,.9),(.8,1.,.1),(1.,.12,0.)]
 for j in range(4):
  v=f['surface%d'%j]
  if not len(v):continue
  p=projected(v,center,scale,size);fac=np.arange(len(v),dtype='i4').reshape(-1,3);color=np.tile(np.array(palette[j],dtype='f4'),(len(v),1));layer,z=draw_mesh(size,p,fac,color);rgba=np.array(layer);rgba[...,3]=np.where((rgba[...,3]>0)&(z>depth),[25,35,46,64][j],0).astype('u1');im.alpha_composite(Image.fromarray(rgba))
 d=ImageDraw.Draw(im)
 for x,y,z in projected(f['particles'],center,scale,size):
  ix,iy=int(x),int(y)
  if 0<=ix<w and 0<=iy<h and z>depth[iy,ix]:d.ellipse((x-1,y-1,x+1,y+1),fill=(160,220,235,190))
 label(d,(24,h-72),'t = %.2f s'%meta['time_s'],25,WHITE);large_legend(d,size,40)
 return im

def frame_v02(t,form='portrait',scale_output=1):
 size=(1080,1920)if form=='portrait'else(1920,1080);w,h=size;portrait=h>w;im=base(size);d=ImageDraw.Draw(im);m=64 if portrait else 84
 idx=next((i for i,(a,b)in enumerate(TIMES)if a<=t<b),5);start,end=TIMES[idx];local=t-start;a=local/(end-start)
 d.line((m,74 if portrait else 54,m+35,74 if portrait else 54),fill=CYAN,width=3);txt(d,(m+53,61 if portrait else 41),'LICA',22,MUTED)
 txt(d,(w-m,61 if portrait else 41),'%02d / %02d'%(idx+1,6),22,CYAN,False,'ra')
 textlines(d,(m,180 if portrait else 110),TITLES[idx],78 if portrait else 76,w-2*m,WHITE,True,1.08)
 box=(44,415,w-44,1720)if portrait else(65,300,w-65,980);aw,ah=box[2]-box[0],box[3]-box[1];art=Image.new('RGBA',(aw,ah));ad=ImageDraw.Draw(art);gait=t*.48
 if idx==0:
  boxes=[(0,0,aw,ah//2-10),(0,ah//2+10,aw,ah)]if portrait else[(0,0,aw//2-20,ah),(aw//2+20,0,aw,ah)]
  for j,b in enumerate(boxes):
   pw,ph=b[2]-b[0],b[3]-b[1];view,_=hero((pw,ph),gait,OFFSETS[:1]if j==0 else OFFSETS,center=(0,.8,0)if j==0 else(1.4,.8,0),scale=min(pw/(3.1 if j==0 else 7.0),ph/2.7),az=40,el=12);art.alpha_composite(view,(b[0],b[1]));label(ad,(b[0]+25,b[1]+22),'SOLO'if j==0 else'PACED',35,GOLD if j==0 else CYAN)
   for yy in(.34,.50,.66):arrow(ad,(b[2]-40,b[1]+ph*yy),(b[2]-120,b[1]+ph*yy),CYAN,2)
 elif idx==1:
  viewsize=(aw,ah-330)if portrait else(aw-500,ah);center=(-1,1.35,0);scale=min(viewsize[0]/19,viewsize[1]/10.5);az=40;el=22;vd=Image.new('RGBA',viewsize);dd=ImageDraw.Draw(vd);ground(dd,viewsize,center,scale,az,el,True)
  for x in(-8,6):
   for y in(0,4):p3(dd,[[x,y,-6],[x,y,6]],center,scale,viewsize,az,el,(40,106,119),2)
   for z in(-6,6):p3(dd,[[x,0,z],[x,4,z]],center,scale,viewsize,az,el,(40,106,119),2)
  for y in(0,4):
   for z in(-6,6):p3(dd,[[-8,y,z],[6,y,z]],center,scale,viewsize,az,el,(40,106,119),2)
  body,_=bodies(viewsize,gait,OFFSETS,center,scale,az,el);vd.alpha_composite(body);art.alpha_composite(vd,(0,0))
  specs=[('21.10','km/h'),('192','spm'),('14 x 12 x 4','m')]
  for j,(value,unit)in enumerate(specs):
   if portrait:x=25 if j!=1 else aw*.54;y=ah-320 if j<2 else ah-155
   else:x=aw-455;y=80+j*170
   label(ad,(x,y),value,59 if j<2 else 46,GOLD if j==0 else WHITE);txt(ad,(x,y+72),unit,31,MUTED)
 elif idx in(2,3):
  center=(0,.87,0);scale=min(aw/2.9,ah/2.7);az=40+3*math.sin(local*.3);view,depth=hero((aw,ah),gait,OFFSETS[:1],center,scale,az,10);gd=ImageDraw.Draw(art)
  spacing=.075 if idx==2 else.14
  for x in np.arange(-1.05,1.06,spacing):p3(gd,[[x,.02,0],[x,1.95,0]],center,scale,(aw,ah),az,10,(39,105,123,110))
  for y in np.arange(.02,1.96,spacing):p3(gd,[[-1.05,y,0],[1.05,y,0]],center,scale,(aw,ah),az,10,(39,105,123,110))
  for z in np.arange(-.8,.81,spacing):p3(gd,[[-1.05,.02,z],[1.05,.02,z]],center,scale,(aw,ah),az,10,(39,105,123,90))
  art.alpha_composite(view)
  if idx==3:
   mask=depth>-1e20;phase=local*.85%1
   for j,yy in enumerate(range(int(ah*.23),int(ah*.79),max(22,ah//24))):
    xs=np.flatnonzero(mask[yy])
    if not len(xs):continue
    xx=int(xs[-1]);tip=(xx+56,yy-15);ad.line(((xx,yy),tip),fill=(37,100,114),width=2);ad.ellipse((xx-5,yy-5,xx+5,yy+5),fill=GOLD);ad.ellipse((tip[0]-4,tip[1]-4,tip[0]+4,tip[1]+4),fill=CYAN)
    q=np.array([xx,yy])*(1-phase)+np.array(tip)*phase;ad.ellipse((q[0]-2.5,q[1]-2.5,q[0]+2.5,q[1]+2.5),fill=WHITE)
  label(ad,(24,ah-50),'SCHEMATIC',27,MUTED)
  eq_begin=4.5 if idx==2 else 5.0
  if local>eq_begin:
   alpha=smooth((local-eq_begin)/.2);overlay=Image.new('RGBA',(aw,ah));od=ImageDraw.Draw(overlay);od.rounded_rectangle((12,ah-320 if portrait else ah-235,aw-12,ah-75),radius=12,fill=(5,10,15,238))
   if idx==2:
    equation(overlay,'CONTINUITY',(25,ah-320 if portrait else ah-235,aw-25,ah-220 if portrait else ah-165));equation(overlay,'MOMENTUM',(25,ah-220 if portrait else ah-165,aw-25,ah-92))
   else:equation(overlay,'IBM',(25,ah-265 if portrait else ah-230,aw-25,ah-100))
   if alpha<1:overlay.putalpha(overlay.getchannel('A').point(lambda v:int(v*alpha)))
   art.alpha_composite(overlay)
 elif idx==4:
  gap=28;boxes=[(0,0,aw,(ah-gap)//2),(0,(ah+gap)//2,aw,ah)]if portrait else[(0,0,(aw-gap)//2,ah),((aw+gap)//2,0,aw,ah)]
  for j,b in enumerate(boxes):
   pw,ph=b[2]-b[0],b[3]-b[1];outline(ad,(b[0],b[1],b[2]-1,b[3]-1),EDGE)
   if j==0:
    panel=Image.new('RGBA',(pw-4,ph-4));pd=ImageDraw.Draw(panel);corners(pd,(16,16,pw-20,ph-20),(44,69,82));label(pd,((pw-4)/2,(ph-4)/2-22),'SOLO',54,WHITE,'mm');label(pd,((pw-4)/2,(ph-4)/2+50),'PENDING',28,CYAN,'mm')
   else:panel=clean_flow((pw-4,ph-4),a)
   art.alpha_composite(panel,(b[0]+2,b[1]+2));label(ad,(b[0]+24,b[1]+20),'SOLO'if j==0 else'PACED',32,CYAN if j else MUTED)
 elif idx==5:
  # Direct comparison of actual pressure contributions; N values remain reserved.
  pw=(aw-25)//2;ph=ah-220;pos=(PC['body']+OFFSETS[0]);normal=PC['normals'];center=(0,.86,0)
  for j in range(2):
   origin=(j*(pw+25),0);v=projected(pos,center,min(pw/1.3,ph/2.1),(pw,ph),40,10);color=PC['soloColor'if j==0 else'pacedColor']
   shade=.48+.52*np.maximum(0,normal@np.array([.45,.8,.38]));color=np.clip(color*shade[:,None],0,1);body,_=draw_mesh((pw,ph),v,TRI,color);art.alpha_composite(body,origin)
   label(ad,(origin[0]+pw/2,25),'SOLO'if j==0 else'PACED',37,WHITE,'ma');label(ad,(origin[0]+pw/2,ph-20),'-- N',50,GOLD,'ma')
  # One shared pressure-contribution scale; no per-case normalization.
  limit=float(PC['range_Pa']);legendwidth=min(aw-120,620);lx=(aw-legendwidth)//2;ly=ah-145
  from matplotlib import colormaps
  for j in range(legendwidth):ad.line((lx+j,ly,lx+j,ly+15),fill=tuple((255*np.array(colormaps['RdBu_r'](j/(legendwidth-1))[:3])).astype(int)))
  label(ad,(aw/2,ly-60),'p nₓ  [Pa]',29,WHITE,'mm');txt(ad,(lx,ly+25),'-%g'%limit,27,MUTED);txt(ad,(lx+legendwidth,ly+25),'+%g'%limit,27,MUTED,False,'ra');label(ad,(aw/2,ah-25),'DRAG VALUES PENDING',27,CYAN,'mm')
 im.alpha_composite(art,(box[0],box[1]));d=ImageDraw.Draw(im);py=1830 if portrait else 1038
 d.line((m,py,w-m,py),fill=EDGE,width=2);d.line((m,py,m+(w-2*m)*(idx+a)/6,py),fill=CYAN,width=3)
 if idx==4:txt(d,(m,py-48),'ACTUAL CFD',27,MUTED)
 fade=min(smooth(local/.22),smooth((end-t)/.18))
 if fade<1:im=Image.blend(base(size),im,fade)
 if scale_output!=1:im=im.resize((int(w*scale_output),int(h*scale_output)),Image.Resampling.LANCZOS)
 return im.convert('RGB')
