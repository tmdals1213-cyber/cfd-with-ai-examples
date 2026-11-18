"""Presentation fix: graph markers use the exact saved surface sample timestamps.

Original v03 source and segments remain immutable. Both cases' original N values
are unchanged. The nearest saved poses may differ in physical time by milliseconds.
"""
import film_renderer_v03 as M
from film_renderer_v03 import *
original_surface=M.pressure_surface
def synced_surface(size,case,time_s):
 im=original_surface(size,case,time_s);row,_=section(case,time_s);d=ImageDraw.Draw(im);txt(d,(size[0]/2,50),'t = %.3f s'%row['time_s'],22,MUTED,False,'ma');return im
def synced_history(size,time_s):
 w,h=size;im=Image.new('RGBA',size);d=ImageDraw.Draw(im);x0,y0,x1,y1=80,65,w-30,h-65
 lo=math.floor(min(SERIES['soloLoad'].min(),SERIES['pacedLoad'].min())/5)*5;hi=math.ceil(max(SERIES['soloLoad'].max(),SERIES['pacedLoad'].max())/5)*5
 def xy(t,f):return(x0+(t-.625)/1.875*(x1-x0),y1-(f-lo)/(hi-lo)*(y1-y0))
 for yy in np.linspace(lo,hi,5):
  y=xy(.625,yy)[1];d.line((x0,y,x1,y),fill=EDGE,width=1);txt(d,(x0-12,y),'%g'%yy,23,MUTED,False,'rm')
 d.line((*xy(.625,0),*xy(2.5,0)),fill=(80,110,125),width=2)
 for x in(.625,1.25,1.875,2.5):txt(d,(xy(x,lo)[0],y1+15),'%.2f'%x,23,MUTED,False,'ma')
 for case,prefix,color in(('solo','solo',GOLD),('group8','paced',CYAN)):
  row,_=section(case,time_s);ts=SERIES[prefix+'Time'];values=SERIES[prefix+'Load'];sel=(ts>.625)&(ts<row['time_s']);xx=np.r_[.625,ts[sel],row['time_s']];vv=np.r_[np.interp(.625,ts,values),values[sel],row['load_N']]
  assert abs(float(np.interp(row['time_s'],ts,values))-row['load_N'])<1e-8
  pts=[xy(t,v)for t,v in zip(xx,vv)]
  if len(pts)>1:d.line(pts,fill=color,width=4)
  x,y=pts[-1];d.ellipse((x-6,y-6,x+6,y+6),fill=color)
 x=xy(time_s,lo)[0];d.line((x,y0,x,y1),fill=(150,174,184),width=1)
 label(d,(x0,8),'PRESSURE LOAD  Dp [N]',27);txt(d,(x1,10),'time [s]',24,MUTED,False,'ra');label(d,(w/2-60,h-12),'SOLO',23,GOLD,'mm');label(d,(w/2+70,h-12),'PACED',23,CYAN,'mm');return im
M.pressure_surface=synced_surface;M.history_plot=synced_history
def frame_v03_synced(t,form='portrait',scale_output=1):return M.frame_v03(t,form,scale_output)
