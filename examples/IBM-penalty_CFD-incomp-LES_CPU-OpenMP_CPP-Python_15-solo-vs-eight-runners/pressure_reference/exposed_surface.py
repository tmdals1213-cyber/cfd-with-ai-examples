"""Fresh, finer exposed-union quadrature of a saved original body pose."""
from pressure_kernel import *
import ctypes as C,struct
from media_common import DP,IP,ptr
class Event(C.Structure):_fields_=[('info',C.c_int*4),('data',C.c_double*13)]
def exposed(cfg,pose_pin,q=95):
 raw=frozen(pose_pin);magic,nv,t,epoch=struct.unpack_from('<8sidQ',raw);assert magic==b'RDXPOS01'
 xyz=np.frombuffer(raw,'<f8',nv*3,28).reshape(-1,3).copy();vel=np.frombuffer(raw,'<f8',nv*3,28+nv*24).reshape(-1,3).copy();top=frozen(cfg['topology']);context=top[20:84]
 assert pin(Path(cfg['geometryDLL']['path']))['SHA']==cfg['geometryDLL']['SHA'];lib=C.CDLL(cfg['geometryDLL']['path'])
 lib.rnu1_rd_load.argtypes=[C.c_char_p,C.c_char_p,C.POINTER(C.c_void_p)];lib.rnu1_rd_load.restype=C.c_int
 lib.rnu1_rd_update.argtypes=[C.c_void_p,C.c_int,DP,DP,C.c_double,C.c_uint64,C.c_char_p];lib.rnu1_rd_update.restype=C.c_int
 lib.rnu1_rd_segment.argtypes=[C.c_void_p,DP,DP,C.c_double,C.c_uint64,C.c_int,C.c_int,IP,IP,IP,IP,C.POINTER(Event)];lib.rnu1_rd_segment.restype=C.c_int;lib.rnu1_rd_release.argtypes=[C.c_void_p]
 handle=C.c_void_p();assert lib.rnu1_rd_load(cfg['topology']['path'].encode('utf8'),context,C.byref(handle))==0
 assert lib.rnu1_rd_update(handle,nv,ptr(xyz,DP),ptr(vel,DP),t,1,context)==0
 lo=xyz.min(0);hi=xyz.max(0);points=[];normals=[];areas=[];vol=[];queries=0
 try:
  for axis in range(3):
   other=[a for a in range(3)if a!=axis];pa=(hi[other[0]]-lo[other[0]])*(hi[other[1]]-lo[other[1]])/(q*q);volume=0.
   for i in range(q):
    for j in range(q):
     a=np.zeros(3);b=np.zeros(3)
     for k,ax in enumerate(other):a[ax]=b[ax]=lo[ax]+((i if k==0 else j)+.5)*(hi[ax]-lo[ax])/q
     a[axis]=lo[axis]-.125;b[axis]=hi[axis]+.125;nr=C.c_int();ne=C.c_int();si=C.c_int();ei=C.c_int();events=(Event*256)()
     code=lib.rnu1_rd_segment(handle,ptr(a,DP),ptr(b,DP),t,1,1,256,C.byref(nr),C.byref(ne),C.byref(si),C.byref(ei),events)
     assert code==0 and si.value==ei.value==0 and ne.value%2==0,('New exposed95 geometry refusal',axis,i,j,code)
     queries+=1
     for k in range(ne.value):
      d=events[k].data;normal=np.array(d[7:10]);assert abs(normal@normal-1)<1e-10 and events[k].info[3]==(-1 if k%2==0 else 1)
      if k%2==1:volume+=(d[1+axis]-events[k-1].data[1+axis])*pa
      if int(np.argmax(abs(normal)))==axis:points.append(list(d[1:4]));normals.append(normal);areas.append(pa/abs(normal[axis]))
   vol.append(volume)
 finally:lib.rnu1_rd_release(handle)
 points=np.array(points);normals=np.array(normals);areas=np.array(areas);assert np.isfinite(points).all()and np.all(areas>0)
 return points,normals,areas,float(np.mean(vol)),dict(queries=queries,pointCount=len(points),axisVolume_m3=vol),xyz
