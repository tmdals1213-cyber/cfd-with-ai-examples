from pathlib import Path
import ctypes as C,hashlib,io,json,os,sys
sys.dont_write_bytecode=True
R=Path(__file__).resolve().parent.parent;WS=R.parents[1]
sys.path.insert(0,str(WS/'work/plot_deps'))
import numpy as np
from PIL import Image,ImageDraw,ImageFont,ImageFilter
DLLDIR=os.add_dll_directory('__LOCAL_PATH_REQUIRES_CONFIGURATION__')
LIB=C.CDLL(str(R/'assets/mesh_media.dll'))
FP=C.POINTER(C.c_float);DP=C.POINTER(C.c_double);IP=C.POINTER(C.c_int);BP=C.POINTER(C.c_ubyte)
LIB.raster.argtypes=[C.c_int,C.c_int,C.c_int,FP,FP,C.c_int,IP,BP,FP,C.c_int];LIB.raster.restype=C.c_int
LIB.normals.argtypes=[C.c_int,DP,C.c_int,IP,IP,FP]
LIB.isosurface.argtypes=[C.c_int]*3+[DP]*3+[DP,BP,C.c_double,C.c_int,FP];LIB.isosurface.restype=C.c_int
def ptr(a,t=FP):return a.ctypes.data_as(t)
def pin(q):
 q=Path(q).resolve();b=q.read_bytes();return dict(path=str(q),bytes=len(b),SHA=hashlib.sha256(b).hexdigest())
def frozen(p):
 b=Path(p['path']).read_bytes();assert len(b)==p['bytes']and hashlib.sha256(b).hexdigest()==p['SHA'];return b
def save(q,x):
 with Path(q).open('x',encoding='utf8')as f:json.dump(x,f,indent=2)
def basis(az=40,el=10):
 a,e=np.deg2rad([az,el]);eye=np.array([np.cos(a)*np.cos(e),np.sin(e),np.sin(a)*np.cos(e)]);right=np.array([np.sin(a),0,-np.cos(a)]);up=np.cross(eye,right);return np.stack([right,up,eye],axis=1)
def projected(x,center,scale,size,az=40,el=10):
 p=(x-np.asarray(center))@basis(az,el);p[:,0]=size[0]/2+scale*p[:,0];p[:,1]=size[1]/2-scale*p[:,1];return np.ascontiguousarray(p,dtype='f4')
def draw_mesh(size,verts,faces,colors,threads=8):
 w,h=size;p=np.ascontiguousarray(verts,dtype='f4');c=np.ascontiguousarray(colors,dtype='f4');tri=np.ascontiguousarray(faces,dtype='i4');rgba=np.zeros((h,w,4),dtype='u1');depth=np.empty((h,w),dtype='f4')
 assert LIB.raster(w,h,len(p),ptr(p),ptr(c),len(tri),ptr(tri,IP),ptr(rgba,BP),ptr(depth),threads)==0
 return Image.fromarray(rgba),depth
def vertex_normals(p,faces,canon):
 p=np.ascontiguousarray(p,dtype='f8');tri=np.ascontiguousarray(faces,dtype='i4');canon=np.ascontiguousarray(canon,dtype='i4');out=np.empty(p.shape,dtype='f4');LIB.normals(len(p),ptr(p,DP),len(tri),ptr(tri,IP),ptr(canon,IP),ptr(out));return out
def lighting(n,color,az=40,el=10):
 light=np.array([.45,.8,.38]);light/=np.linalg.norm(light);shade=.26+.67*np.maximum(0,n@light);eye=basis(az,el)[:,2];rim=np.maximum(0,1-np.abs(n@eye))**3;return np.minimum(1,np.asarray(color)*shade[:,None]+rim[:,None]*.14).astype('f4')
def grid_info(p):
 with np.load(io.BytesIO(frozen(p)),allow_pickle=False)as z:faces=[np.ascontiguousarray(z['face%d'%a])for a in range(3)]
 n=tuple(len(f)-1 for f in faces);h=[np.diff(f)for f in faces];cent=[.5*(f[1:]+f[:-1])for f in faces];shapes=[tuple(n[b]+int(a==b)for b in range(3))for a in range(3)];sizes=[int(np.prod(s))for s in shapes];off=np.r_[0,np.cumsum(sizes)];m=[]
 for a in range(3):
  v=np.ones(shapes[a]);q=[1]*3
  for b in range(3):
   if b!=a:q=[1]*3;q[b]=n[b];v*=h[b].reshape(q)
  span=np.r_[h[a][0]/2,(h[a][:-1]+h[a][1:])/2,h[a][-1]/2];q=[1]*3;q[a]=n[a]+1;m.append(1.225*(v*span.reshape(q)))
 md=np.concatenate([q.ravel()for q in m]);assert hashlib.sha256(md.tobytes()).hexdigest()==p['sourcePin']['MdSHA'];return faces,cent,n,shapes,off,md
def velocity(momentum,shapes,off,md):
 f=momentum/md;v=[]
 for a in range(3):
  u=f[off[a]:off[a+1]].reshape(shapes[a]);l=[slice(None)]*3;r=l.copy();l[a]=slice(None,-1);r[a]=slice(1,None);v.append(.5*(u[tuple(l)]+u[tuple(r)]))
 return np.stack(v,axis=-1)
def trilinear(v,points,cent):
 if not len(points):return np.empty((0,3))
 ii=[np.clip(np.searchsorted(cent[a],points[:,a],side='right')-1,0,len(cent[a])-2)for a in range(3)];s=[(points[:,a]-cent[a][ii[a]])/(cent[a][ii[a]+1]-cent[a][ii[a]])for a in range(3)];out=np.zeros((len(points),3))
 for a in(0,1):
  for b in(0,1):
   for c in(0,1):out+=v[ii[0]+a,ii[1]+b,ii[2]+c]*(np.where(a,s[0],1-s[0])*np.where(b,s[1],1-s[1])*np.where(c,s[2],1-s[2]))[:,None]
 return out
