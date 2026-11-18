import sys,io,json,hashlib,time,ctypes as C
from pathlib import Path
sys.dont_write_bytecode=True
H=Path(__file__).resolve().parent
READ=0
def need(x,m):
 if not x:raise RuntimeError(m)
def sha(b):return hashlib.sha256(b).hexdigest()
def frozen(pin):
 global READ
 q=Path(pin['path']);a=q.stat();need(q.is_file()and not q.is_symlink()and a.st_size==pin['bytes']and a.st_size<=134217728,'Exact bounded input file');b=q.read_bytes();z=q.stat();READ+=len(b);need((a.st_size,a.st_mtime_ns)==(z.st_size,z.st_mtime_ns)and sha(b)==pin['SHA'],'Stable exact input SHA');return b
def save(name,x):
 with(H/name).open('x',encoding='utf8')as f:json.dump(x,f,indent=2)
def pinfile(name):
 q=H/name;b=q.read_bytes();need(len(b)<=134217728,'Bounded output');return dict(path=str(q),bytes=len(b),SHA=sha(b))
def np_save(name,x,np):
 with(H/name).open('xb')as f:np.save(f,x,allow_pickle=False)
 return pinfile(name)
def start():
 cfg=json.loads((H/'BASE_LOCAL_INPUT.json').read_text('utf8'))
 for pin in cfg['codePins']:frozen(pin)
 for pin in json.loads(frozen(cfg['numpyManifest']))['files']:frozen(pin)
 sys.path.insert(0,str(H.parents[1]/'plot_deps'))
 import numpy as np
 need(np.dtype('f8').itemsize==8 and C.sizeof(C.c_double)==8,'Portable f64 ABI')
 return np,cfg,time.monotonic()
def grid(np,pin):
 with np.load(io.BytesIO(frozen(pin)),allow_pickle=False)as z:faces=[np.ascontiguousarray(z['face%d'%a])for a in range(3)]
 n=tuple(len(f)-1 for f in faces);h=[np.diff(f)for f in faces];shapes=[tuple(n[b]+int(a==b)for b in range(3))for a in range(3)];sizes=[int(np.prod(q))for q in shapes];off=np.r_[0,np.cumsum(sizes)];D=int(off[-1]);N=int(np.prod(n));mass=[];areas=[]
 for a in range(3):
  ar=np.ones(shapes[a])
  for b in range(3):
   if b!=a:q=[1,1,1];q[b]=n[b];ar*=h[b].reshape(q)
  span=np.r_[h[a][0]/2,(h[a][:-1]+h[a][1:])/2,h[a][-1]/2];q=[1,1,1];q[a]=n[a]+1;mass.append(1.225*(ar*span.reshape(q)));areas.append(ar)
 md=np.concatenate([q.ravel()for q in mass]);need(sha(md.tobytes())==pin['sourcePin']['MdSHA'],'Exact accepted face dual masses')
 return faces,n,N,h,shapes,sizes,off,D,mass,areas,md
def dll(pin):
 frozen(pin);return C.CDLL(pin['path'])
DP=C.POINTER(C.c_double);IP=C.POINTER(C.c_int);LP=C.POINTER(C.c_longlong)
def ptr(x):return x.ctypes.data_as(DP)

import mmap,gc
OPEN_MAPS=[]
def shared_arrays(np,D,N,K,desc,role):
 need(role in('prepare','pressure','audit','control'),'Declared shared-buffer role');L=N+K
 guard=mmap.mmap(-1,64,tagname=desc['guard']['name'],access=mmap.ACCESS_READ);OPEN_MAPS.append(guard);need(guard[:64]==bytes.fromhex(desc['guardHex']),'Exact live buffer owner')
 shapes=dict(prepared=(10*D+N,),rhs=(L,),state=(3*L+32,),trace=(385,9),output=(3*D+2*L,),stats=(12,),predictorStats=(6,));result={}
 for name,shape in shapes.items():
  size=8
  for dim in shape:size*=dim
  need(desc[name]['bytes']==size,'Exact typed shared buffer size')
  ro=role=='audit'or(role=='pressure'and name in('prepared','rhs','predictorStats'))
  mm=mmap.mmap(-1,size,tagname=desc[name]['name'],access=mmap.ACCESS_READ if ro else mmap.ACCESS_WRITE);OPEN_MAPS.append(mm);result[name]=np.ndarray(shape,dtype='f8',buffer=mm)
 return result
def close_shared():
 gc.collect()
 for mm in reversed(OPEN_MAPS):mm.close()
 OPEN_MAPS.clear()
