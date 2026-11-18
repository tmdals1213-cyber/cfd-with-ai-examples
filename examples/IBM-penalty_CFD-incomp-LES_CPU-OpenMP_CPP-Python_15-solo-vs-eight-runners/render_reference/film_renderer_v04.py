"""Public-facing edit: high resolution pressure explanation, no production labels."""
from film_renderer_v03_synced import *
V4=R/'v04';HR=V4/'assets/pressure';PRESSURE_ASSETS={}
for case in('solo','group8'):
 PRESSURE_ASSETS[case]=(Image.open(HR/(case+'_PRESSURE.png')).convert('RGBA'),Image.open(HR/(case+'_BODY.png')).convert('RGBA'),json.loads((HR/(case+'_SOURCE.json')).read_text()))
def pressure_explainer(size,case,local):
 w,h=size;background,body,meta=PRESSURE_ASSETS[case];im=background.copy();d=ImageDraw.Draw(im)
 for k,path in enumerate(meta['streamlinePaths']):
  points=np.array(path);d.line([tuple(p)for p in points],fill=(34,126,133,150),width=2)
  length=np.r_[0,np.cumsum(np.linalg.norm(np.diff(points,axis=0),axis=1))]
  if length[-1]<30:continue
  for frac in(.18,.63):
   u=((local*.16+k*.071+frac)%1)*length[-1];q=np.array([np.interp(u,length,points[:,j])for j in range(2)]);r=np.array([np.interp(min(u+18,length[-1]),length,points[:,j])for j in range(2)])
   if np.linalg.norm(r-q)>3:arrow(d,q,r,(20,117,126,210),3,8)
 im.alpha_composite(body)
 if case=='group8':
  # Relative pressure ordering around P0 is supported by the native plane probes.
  # These are explanatory callouts, not new pressure or force operands.
  for x0,word,color in((-.55,'HIGHER',(181,87,64)),(.58,'LOWER',(25,100,148))):
   x=(x0+.95)/3.7*im.width;y=.16/2.25*im.height
   txt(d,(x,y),word,46,color,True,'mm');arrow(d,(x,y+35),(x,((2.25-1.47)/2.25)*im.height),color,3,12)
 im=im.resize(size,Image.Resampling.LANCZOS)
 return im
def frame_v04(t,form='portrait',scale_output=1):
 size=(1080,1920)if form=='portrait'else(1920,1080);w,h=size;portrait=h>w;m=64 if portrait else 84;idx=next((i for i,(a,b)in enumerate(STAGES)if a<=t<b),7);start,end=STAGES[idx];local=t-start;bg=base(size)
 if idx==5:
  im=bg.copy();d=ImageDraw.Draw(im);d.line((m,74 if portrait else 54,m+35,74 if portrait else 54),fill=CYAN,width=3);txt(d,(w-m,61 if portrait else 41),'06 / 08',22,CYAN,False,'ra')
  textlines(d,(m,180 if portrait else 110),'Pressure explains\nthe difference.',78 if portrait else 76,w-2*m,WHITE,True,1.08)
  box=(44,415,w-44,1720)if portrait else(65,300,w-65,980);aw,ah=box[2]-box[0],box[3]-box[1];art=Image.new('RGBA',(aw,ah));ad=ImageDraw.Draw(art)
  boxes=[(0,0,aw,(ah-34)//2),(0,(ah+34)//2,aw,ah)]if portrait else[(0,0,(aw-28)//2,ah),((aw+28)//2,0,aw,ah)]
  for case,b in zip(('solo','group8'),boxes):
   pw,ph=b[2]-b[0],b[3]-b[1];meta=PRESSURE_ASSETS[case][2];factor=min(pw/1800,(ph-130)/1094);rw,rh=int(1800*factor),int(1094*factor);x=b[0]+(pw-rw)//2;y=b[1]+48
   art.alpha_composite(pressure_explainer((rw,rh),case,local),(x,y));label(ad,(b[0]+20,b[1]+7),'SOLO'if case=='solo'else'PACED',32,GOLD if case=='solo'else CYAN);txt(ad,(b[2]-20,b[1]+10),'t = %.3f s'%meta['time_s'],23,MUTED,False,'ra')
   fx=meta['forceX_N'];cy=y+rh+25;cx=b[0]+pw/2;length=abs(fx)*13;sign=1 if fx>=0 else-1;arrow(ad,(cx-sign*length/2,cy),(cx+sign*length/2,cy),GOLD,5,13)
   label(ad,(cx,cy+35),'PRESSURE FORCE  %+.2f N'%fx,28,WHITE,'mm')
  im.alpha_composite(art,(box[0],box[1]));d=ImageDraw.Draw(im);py=1830 if portrait else 1038
  lx=w-315;ly=1788 if portrait else 1008;bw=230
  for k in range(bw):d.line((lx+k,ly,lx+k,ly+9),fill=tuple(LUT[int(k/(bw-1)*255)]))
  txt(d,(lx,ly-29),'p − p_ref  [Pa]',22,MUTED,True);txt(d,(lx,ly+12),'−50',18,MUTED);txt(d,(lx+bw,ly+12),'+50',18,MUTED,False,'ra')
  d.line((m,py,w-m,py),fill=EDGE,width=2);d.line((m,py,m+(w-2*m)*(idx+local/(end-start))/8,py),fill=CYAN,width=3)
 else:
  im=frame_v03_synced(t,form).convert('RGBA');py=1830 if portrait else 1038
  # Remove production descriptors, leaving geometry, equations and quantities.
  rect=(m,py-51,w-m,py-5);im.paste(bg.crop(rect),rect)
  if idx in(2,3):
   rect=(60 if portrait else 78,1660 if portrait else 920,350,1720 if portrait else 980);im.paste(bg.crop(rect),rect)
  if idx==7:
   rect=(m,150 if portrait else 100,w-m,380 if portrait else 288);im.paste(bg.crop(rect),rect);d=ImageDraw.Draw(im)
   textlines(d,(m,180 if portrait else 110),'One runner.\nDifferent flow.',78 if portrait else 76,w-2*m,WHITE,True,1.08)
   rect=(44 if portrait else 65,1660 if portrait else 930,w-44 if portrait else w-65,1720 if portrait else 980);im.paste(bg.crop(rect),rect);d=ImageDraw.Draw(im);label(d,(w/2,1695 if portrait else 958),'PRESSURE CONTRIBUTION',28,MUTED,'mm')
 if scale_output!=1:im=im.resize((int(w*scale_output),int(h*scale_output)),Image.Resampling.LANCZOS)
 return im.convert('RGB')
