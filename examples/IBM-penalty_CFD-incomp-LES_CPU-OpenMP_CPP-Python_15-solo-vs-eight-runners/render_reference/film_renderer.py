"""English editorial animation; illustrative explanation precedes actual CFD."""
from media_common import *
import math

WHITE=(234,242,247);MUTED=(135,158,172);CYAN=(68,211,222);GOLD=(245,181,83);EDGE=(26,51,65);BG=(5,10,15)
SCENES=[(0,7,'THE QUESTION',"Can a team change\none runner's drag?",'Same athlete. Same speed. A different aerodynamic environment.'),
 (7,18,'THE PHYSICS','Define a fair\ncomparison.','Match speed, gait, ground motion and the surrounding air.'),
 (18,29,'THE CARTESIAN GRID','The body moves.\nThe grid stays.','Resolve the fluid on a fixed Cartesian grid.'),
 (29,41,'THE IMMERSED BOUNDARY','A moving surface.\nA force on the fluid.','Couple the prescribed body motion to the surrounding fluid.'),
 (41,53,'THE FLOW FIELD','Same time. Same view.\nSame scales.','Compare the incoming flow around the protected runner.'),
 (53,60,'THE COMPARISON','Measure the change\non one runner.','Compare the protected runner, not the sum of the whole team.')]
LABELS=['QUESTION','PHYSICS','GRID','IBM','FLOW','COMPARE']
FONTS={};BASE={};FLOW={}
GEOMETRY=np.load(R/'assets/INTRO_GEOMETRY.npz',allow_pickle=False)
POSES=GEOMETRY['poses'];NORMALS=GEOMETRY['normals'];TRI=np.ascontiguousarray(GEOMETRY['triangles'],dtype='i4');OFFSETS=GEOMETRY['offsets'];CANON=GEOMETRY['canonical']
FLOWROOT=R/'assets/group8_flow_v2'

def font(size,bold=False):
 key=(size,bold)
 if key not in FONTS:FONTS[key]=ImageFont.truetype('C:/Windows/Fonts/arialbd.ttf'if bold else'C:/Windows/Fonts/arial.ttf',size)
 return FONTS[key]
def txt(draw,xy,text,size=24,fill=WHITE,bold=False,anchor=None):draw.text(xy,text,font=font(size,bold),fill=fill,anchor=anchor)
def wrap(draw,text,size,width,bold=False):
 words=text.split();lines=[];line=''
 for word in words:
  candidate=(line+' '+word).strip()
  if line and draw.textlength(candidate,font=font(size,bold))>width:lines.append(line);line=word
  else:line=candidate
 if line:lines.append(line)
 return lines
def smooth(x):x=max(0.,min(1.,x));return x*x*(3-2*x)
def textlines(draw,xy,text,size,width,fill=WHITE,bold=False,leading=1.2):
 lines=[]
 for part in text.split('\n'):lines.extend(wrap(draw,part,size,width,bold))
 for i,line in enumerate(lines):txt(draw,(xy[0],xy[1]+i*size*leading),line,size,fill,bold)
 return len(lines)*size*leading
def base(size):
 if size not in BASE:
  w,h=size;yy,xx=np.mgrid[0:h,0:w];glow=np.exp(-((xx/w-.66)**2/.35+(yy/h-.55)**2/.25));arr=np.empty((h,w,3),dtype='u1')
  for a,c in enumerate(BG):arr[...,a]=(c+glow*[2,5,8][a]).astype('u1')
  BASE[size]=Image.fromarray(arr).convert('RGBA')
 return BASE[size].copy()
def outline(draw,rect,color=EDGE,width=2):draw.rounded_rectangle(rect,radius=16,outline=color,width=width)
def corners(draw,box,color=CYAN):
 x,y,X,Y=box;s=28
 for a,b,c,d in[(x,y,x+s,y),(x,y,x,y+s),(X,y,X-s,y),(X,y,X,y+s),(x,Y,x+s,Y),(x,Y,x,Y-s),(X,Y,X-s,Y),(X,Y,X,Y-s)]:draw.line((a,b,c,d),fill=color,width=2)
def p3(draw,points,center,scale,size,az=40,el=10,fill=EDGE,width=1):
 q=projected(np.array(points,dtype='f8'),center,scale,size,az,el);draw.line([tuple(x[:2])for x in q],fill=fill,width=width)
def arrow(draw,a,b,color=CYAN,width=2,head=10):
 draw.line((tuple(a),tuple(b)),fill=color,width=width);v=np.array(a,dtype=float)-b;v/=max(1e-10,np.linalg.norm(v));u=np.array([-v[1],v[0]]);p=np.asarray(b);draw.polygon([tuple(p),tuple(p+head*v+head*.45*u),tuple(p+head*v-head*.45*u)],fill=color)
def ground(draw,size,center,scale,az=40,el=10,wide=False):
 lim=7 if wide else 2
 for x in np.arange(-lim,lim+.001,1 if wide else .25):p3(draw,[[x,0,-lim],[x,0,lim]],center,scale,size,az,el,(15,33,43))
 for z in np.arange(-lim,lim+.001,1 if wide else .25):p3(draw,[[-lim,0,z],[lim,0,z]],center,scale,size,az,el,(15,33,43))
def bodies(size,time_s,offsets,center,scale,az=40,el=10,body=None,normals=None):
 if body is None:
  j=int(time_s/.625*30)%30;body=POSES[j];normals=NORMALS[j]
 vertices=[];colors=[];tris=[]
 for j,offset in enumerate(offsets):
  vertices.append(projected(body+offset,center,scale,size,az,el));colors.append(lighting(normals,np.array(GOLD if j==0 else(62,132,146))/255,az,el));tris.append(TRI+j*len(body))
 return draw_mesh(size,np.vstack(vertices),np.vstack(tris),np.vstack(colors))
def hero(size,t,offsets,center=(0,.9,0),scale=None,az=40,el=10):
 w,h=size;scale=scale or min(w/3,h/2.5);im=Image.new('RGBA',size);d=ImageDraw.Draw(im);ground(d,size,center,scale,az,el,wide=len(offsets)>1)
 mesh,dep=bodies(size,t,offsets,center,scale,az,el);im.alpha_composite(mesh);return im,dep
def panel_title(draw,box,title,minor,color=CYAN,size=24):
 txt(draw,(box[0]+24,box[1]+20),title,size,color,True);txt(draw,(box[0]+24,box[1]+54),minor,max(15,size-6),MUTED)
def placeholder(size):
 w,h=size;im=Image.new('RGBA',size);d=ImageDraw.Draw(im);corners(d,(18,18,w-18,h-18),(44,69,82))
 for y in range(120,h-100,55):d.line((65,y,w-65,y),fill=(13,28,37),width=1)
 txt(d,(w/2,h/2-50),'SOLO',min(58,int(w*.1)),WHITE,True,'mm');txt(d,(w/2,h/2+22),'RESULT SLOT',min(23,int(w*.045)),CYAN,True,'mm')
 txt(d,(w/2,h/2+65),'Awaiting matched flow data',min(20,int(w*.04)),MUTED,False,'mm');return im
def flow_frame(k):
 if k not in FLOW:
  with np.load(FLOWROOT/('FLOW_%04d.npz'%k),allow_pickle=False)as z:FLOW[k]={name:z[name].copy()for name in z.files}
  if len(FLOW)>4:del FLOW[next(iter(FLOW))]
 return FLOW[k]
def result(size,a):
 cache=json.loads((FLOWROOT/'FLOW_CACHE_PASS.json').read_text());frames=cache['frames'];frame=frames[int(round(min(1,max(0,a))*(len(frames)-1)))];f=flow_frame(frame['sourceIndex'])
 w,h=size;center=(.65,.83,0);scale=min(w/7.4,h/3.4);im=Image.new('RGBA',size);d=ImageDraw.Draw(im);ground(d,size,center,scale,wide=True)
 mesh,bodydepth=bodies(size,0,OFFSETS,center,scale,body=f['body'],normals=f['normals']);im.alpha_composite(mesh)
 palette=[(0.,.30,1.),(0.,1.,.9),(.8,1.,.1),(1.,.12,0.)]
 for j in range(4):
  v=f['surface%d'%j]
  if not len(v):continue
  screen=projected(v,center,scale,size);faces=np.arange(len(v),dtype='i4').reshape(-1,3);color=np.tile(np.array(palette[j],dtype='f4'),(len(v),1));layer,z=draw_mesh(size,screen,faces,color)
  ar=np.array(layer);ar[...,3]=np.where((ar[...,3]>0)&(z>bodydepth),[25,35,46,64][j],0).astype('u1');im.alpha_composite(Image.fromarray(ar))
 d=ImageDraw.Draw(im);points=projected(f['particles'],center,scale,size)
 for x,y,z in points:
  ix,iy=int(x),int(y)
  if 0<=ix<w and 0<=iy<h and z>bodydepth[iy,ix]:d.ellipse((x-1.1,y-1.1,x+1.1,y+1.1),fill=(160,220,235,190))
 txt(d,(20,h-59),'ACTUAL CFD  /  t = %.3f s'%frame['time_s'],17,WHITE,True);txt(d,(20,h-32),'~0.10 m grid  /  slow motion x6.4',15,MUTED)
 # Shared fixed whole-set vorticity scale; values do not change per frame.
 bw=min(200,w//3);bx=w-bw-20;by=h-41
 for j in range(bw):
  u=j/(bw-1);c=np.clip(1.5-np.abs(4*u-np.array([3,2,1])),0,1);d.line((bx+j,by,bx+j,by+7),fill=tuple((255*c).astype(int)),width=1)
 txt(d,(bx,by-25),'|omega|  [1/s]',15,MUTED);txt(d,(bx,by+11),'0',13,MUTED);txt(d,(bx+bw,by+11),'40',13,MUTED,False,'ra')
 return im

def frame(t,format='portrait',scale_output=1.):
 size=(1080,1920)if format=='portrait'else(1920,1080);w,h=size;portrait=h>w;im=base(size);d=ImageDraw.Draw(im)
 margin=64 if portrait else 84;titley=180 if portrait else 110;titlefont=76 if portrait else 74
 idx=next((i for i,(a,b,*_)in enumerate(SCENES)if a<=t<b),5);start,end,chapter,title,subtitle=SCENES[idx];local=t-start;duration=end-start;a=local/duration
 d.line((margin,75 if portrait else 55,margin+35,75 if portrait else 55),fill=CYAN,width=3)
 txt(d,(margin+52,62 if portrait else 42),'LICA  /  AERODYNAMICS',20,MUTED)
 txt(d,(w-margin,62 if portrait else 42),'%02d / %s'%(idx+1,chapter),18,CYAN,False,'ra')
 textlines(d,(margin,titley),title,titlefont,w-2*margin,WHITE,True,1.08)
 art=(44,410,w-44,1530)if portrait else(65,295,w-65,855);aw,ah=art[2]-art[0],art[3]-art[1]
 canvas=Image.new('RGBA',(aw,ah));cd=ImageDraw.Draw(canvas)
 label='';gait=t*.36
 if idx==0:
  if portrait:
   bw=aw;bh=ah//2-15;boxes=[(0,0,bw,bh),(0,bh+25,bw,ah)]
  else:boxes=[(0,0,aw//2-20,ah),(aw//2+20,0,aw,ah)]
  for j,box in enumerate(boxes):
   pw,ph=box[2]-box[0],box[3]-box[1];view,_=hero((pw,ph),gait,OFFSETS[:1]if j==0 else OFFSETS,center=(0,.8,0)if j==0 else(1.4,.8,0),scale=min(pw/(3.1 if j==0 else 7.0),ph/2.8),az=40,el=12);canvas.alpha_composite(view,(box[0],box[1]));panel_title(cd,box,'SOLO'if j==0 else'WITH PACEMAKERS','Protected runner P0'if j==0 else'7 pacemakers + protected runner',GOLD if j==0 else CYAN,23)
   for yy in [.32,.46,.60]:arrow(cd,(box[2]-45,box[1]+ph*yy),(box[2]-120,box[1]+ph*yy),CYAN,2,9)
  label='FORMATION & AIRFLOW ILLUSTRATION'
 elif idx==1:
  # First establish actual synchronized bodies, then pull back to the real domain.
  zoom=smooth((local-3)/2.5);viewsize=(aw,ah-190)if portrait else(aw-490,ah);center=np.array([1.35,.9,0])*(1-zoom)+np.array([-1,1.4,0])*zoom
  sc=(min(viewsize[0]/7,viewsize[1]/2.8)*(1-zoom)+min(viewsize[0]/20,viewsize[1]/12)*zoom)
  az=40;el=12+10*zoom;vd=Image.new('RGBA',viewsize);vdraw=ImageDraw.Draw(vd);ground(vdraw,viewsize,center,sc,az,el,True)
  if zoom>.01:
   for x in(-8,6):
    for y in(0,4):p3(vdraw,[[x,y,-6],[x,y,6]],center,sc,viewsize,az,el,(35,99,113),2)
    for z in(-6,6):p3(vdraw,[[x,0,z],[x,4,z]],center,sc,viewsize,az,el,(35,99,113),2)
   for y in(0,4):
    for z in(-6,6):p3(vdraw,[[-8,y,z],[6,y,z]],center,sc,viewsize,az,el,(35,99,113),2)
   for z in(-3,0,3):
    p=projected(np.array([[6,1,z],[4.6,1,z]],float),center,sc,viewsize,az,el);arrow(vdraw,p[0,:2],p[1,:2],CYAN,2)
  m,_=bodies(viewsize,gait,OFFSETS,center,sc,az,el);vd.alpha_composite(m);canvas.alpha_composite(vd,(0,0))
  cards=[('21.0975','km/h  /  prescribed speed'),('192','spm  /  synchronized gait'),('14 x 12 x 4','m  /  common domain'),('MOVING','ground  /  runner frame')]
  if portrait:
   for j,(value,unit)in enumerate(cards):
    x=(j%2)*(aw//2);y=ah-210+(j//2)*104;txt(cd,(x+18,y),value,37,GOLD if j==0 else WHITE,True);txt(cd,(x+18,y+48),unit,18,MUTED)
  else:
   for j,(value,unit)in enumerate(cards):x=aw-440;y=22+j*128;txt(cd,(x,y),value,49,GOLD if j==0 else WHITE,True);txt(cd,(x,y+64),unit,21,MUTED)
  label='ACTUAL DOMAIN & PRESCRIBED CONDITIONS'
 elif idx in(2,3):
  viewsize=(aw,ah);center=(0,.87,0);sc=min(aw/3.1,ah/2.7);az=40+3*math.sin(local*.22);el=10
  view,dep=hero(viewsize,gait,OFFSETS[:1],center,sc,az,el);grid=Image.new('RGBA',viewsize);gd=ImageDraw.Draw(grid)
  spacing=.09 if idx==2 else.14
  for x in np.arange(-1.05,1.06,spacing):p3(gd,[[x,.02,0],[x,1.95,0]],center,sc,viewsize,az,el,(39,105,123,105),1)
  for y in np.arange(.02,1.96,spacing):p3(gd,[[-1.05,y,0],[1.05,y,0]],center,sc,viewsize,az,el,(39,105,123,105),1)
  for z in np.arange(-.8,.81,spacing):p3(gd,[[-1.05,.02,z],[1.05,.02,z]],center,sc,viewsize,az,el,(39,105,123,100),1)
  canvas.alpha_composite(grid);canvas.alpha_composite(view)
  if idx==2:
   # A second grid plane is a deliberately enhanced explanatory display.
   canvas.alpha_composite(grid);txt(cd,(25,ah-95),'FIXED FLUID GRID',24,CYAN,True);txt(cd,(25,ah-52),'Fine grid shown schematically',21,MUTED);label='ILLUSTRATIVE GRID / DISPLAY DENSITY ENHANCED'
  else:
   mask=dep>-1.e20;valid=[]
   for yy in range(int(ah*.23),int(ah*.82),max(14,ah//28)):
    xs=np.flatnonzero(mask[yy]);
    if len(xs):valid.append((int(xs[-1]),yy))
   phase=(local*.65)%1
   for j,(xx,yy)in enumerate(valid):
    tip=(xx+45+12*math.sin(j),yy-18);cd.line(((xx,yy),tip),fill=(37,100,114),width=1);cd.ellipse((xx-4,yy-4,xx+4,yy+4),fill=GOLD);cd.ellipse((tip[0]-3,tip[1]-3,tip[0]+3,tip[1]+3),fill=CYAN)
    if j%3==0:
     q=np.array([xx,yy])*(1-phase)+np.array(tip)*phase;cd.ellipse((q[0]-2,q[1]-2,q[0]+2,q[1]+2),fill=WHITE)
   txt(cd,(25,ah-128),'BODY VELOCITY  ->  FLUID FORCE',22,CYAN,True);txt(cd,(25,ah-82),'Surface-penalty coupling',26,WHITE);txt(cd,(25,ah-42),'The moving surface acts on a fixed grid.',20,MUTED);label='ILLUSTRATIVE IBM COUPLING'
 elif idx==4:
  if portrait:gap=26;ph=(ah-gap)//2;boxes=[(0,0,aw,ph),(0,ph+gap,aw,ah)]
  else:gap=36;pw=(aw-gap)//2;boxes=[(0,0,pw,ah),(pw+gap,0,aw,ah)]
  for j,box in enumerate(boxes):
   pw,ph=box[2]-box[0],box[3]-box[1];outline(cd,(box[0],box[1],box[2]-1,box[3]-1),EDGE)
   panel=placeholder((pw-4,ph-4))if j==0 else result((pw-4,ph-4),a);canvas.alpha_composite(panel,(box[0]+2,box[1]+2))
   if j:txt(cd,(box[0]+22,box[1]+18),'WITH PACEMAKERS',22,CYAN,True)
  label='SOLO RESERVED / PACED PANEL USES ACTUAL SAVED CFD'
 else:
  boxes=[]
  if portrait:
   boxw=aw;boxh=190
   for j in range(3):boxes.append((16,45+j*238,aw-16,235+j*238))
  else:
   for j in range(3):x=18+j*(aw//3);boxes.append((x,70,x+aw//3-36,390))
  for j,box in enumerate(boxes):
   outline(cd,box,(44,66,76));txt(cd,(box[0]+25,box[1]+25),['SOLO DRAG','PACED DRAG','REDUCTION'][j],23,MUTED,True)
   txt(cd,(box[0]+28,box[1]+70),'--'if j<2 else'-- %',70 if portrait else 90,GOLD if j==2 else WHITE,True)
   txt(cd,(box[0]+25,box[3]-39),'P0 only  /  N'if j<2 else'Matched time window',20,MUTED)
  fy=ah-150 if portrait else ah-105;txt(cd,(aw/2,fy),'(D_solo - D_paced) / D_solo x 100%',27 if portrait else 34,WHITE,False,'mm');txt(cd,(aw/2,fy+62),'Matched surface-force evaluation pending',23,CYAN,False,'mm');label='COMPARISON SLOT / NO FORCE VALUES INSERTED'
 im.alpha_composite(canvas,(art[0],art[1]));d=ImageDraw.Draw(im)
 statusy=1563 if portrait else 874;txt(d,(margin,statusy),label,16 if portrait else 17,CYAN if idx not in(2,3)else MUTED)
 sy=1620 if portrait else 920;textlines(d,(margin,sy),subtitle,34 if portrait else 32,w-2*margin,WHITE,False,1.3)
 if idx==1:txt(d,(margin,1734 if portrait else 974),'rho = 1.225 kg/m3   /   nu = 1.5e-5 m2/s',20 if portrait else 22,MUTED)
 if idx==5:txt(d,(margin,1734 if portrait else 974),'Exploratory CFD  /  Surface-penalty IBM + Vreman SGS',19 if portrait else 21,MUTED)
 py=1835 if portrait else 1028;d.line((margin,py-20,w-margin,py-20),fill=EDGE,width=1)
 step=(w-2*margin)/6
 for j,label in enumerate(LABELS):
  txt(d,(margin+j*step,py),label,16 if portrait else 18,CYAN if j==idx else(70,91,104),j==idx)
  if j==idx:d.line((margin+j*step,py-20,margin+(j+1)*step-15,py-20),fill=CYAN,width=3)
 # Quiet editorial dip at cuts; no fabricated blending of CFD time states.
 fade=min(smooth(local/.32),smooth((end-t)/.25))
 if fade<1:im=Image.blend(base(size),im,fade)
 if scale_output!=1:im=im.resize((int(w*scale_output),int(h*scale_output)),Image.Resampling.LANCZOS)
 return im.convert('RGB')
