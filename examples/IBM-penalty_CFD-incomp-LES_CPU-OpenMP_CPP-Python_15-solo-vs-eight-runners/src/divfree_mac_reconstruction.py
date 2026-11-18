"""Boundary-aware H(div) polynomial: interior anchor selection, actual wall point unchanged; no velocity clipping."""
"""Affine-exact conservative H(div) MAC reconstruction, derived from shared linear face traces."""
import numpy as np

def layout(faces):
 n=[len(x)-1 for x in faces];shapes=[tuple(n[b]+int(a==b)for b in range(3))for a in range(3)];sizes=[int(np.prod(s))for s in shapes];off=np.r_[0,np.cumsum(sizes)];cent=[.5*(x[:-1]+x[1:])for x in faces];return n,shapes,off,cent

def rows(faces,points,component,normal_derivative=False):
 n,shapes,off,cent=layout(faces);points=np.asarray(points,dtype='f8');a=int(component);K=len(points);cell=np.column_stack([np.searchsorted(faces[b],points[:,b],side='right')-1 for b in range(3)]);physicalCell=cell.copy();cell=np.column_stack([np.minimum(np.maximum(cell[:,b],1),n[b]-2)for b in range(3)]);assert 0<=a<3 and points.shape==(K,3)and np.isfinite(points).all();assert all(np.all((physicalCell[:,b]>=0)&(physicalCell[:,b]<n[b]))for b in range(3));width=np.column_stack([faces[b][cell[:,b]+1]-faces[b][cell[:,b]]for b in range(3)]);lower=np.column_stack([faces[b][cell[:,b]]for b in range(3)]);f=(points-lower)/width;xi=f[:,a]-.5;terms={}
 def add(comp,delta,value):
  key=(comp,tuple(delta));terms[key]=terms.get(key,0.)+value
 def derivative(b):
  j=cell[:,b];hm=cent[b][j]-cent[b][j-1];hp=cent[b][j+1]-cent[b][j];return(-hp/(hm*(hm+hp)),(hp-hm)/(hm*hp),hm/(hp*(hm+hp)))
 for plane in(0,1):
  delta=[0,0,0];delta[a]=plane;v=np.full(K,1. if plane else-1.)/width[:,a]if normal_derivative else f[:,a]if plane else 1-f[:,a];add(a,delta,v)
 for b in range(3):
  if b==a:continue
  db=derivative(b)
  for plane in(0,1):
   own=(np.full(K,1. if plane else-1.)/width[:,a])if normal_derivative else(f[:,a]if plane else 1-f[:,a])
   for shift,d in zip((-1,0,1),db):
    delta=[0,0,0];delta[a]=plane;delta[b]=shift;add(a,delta,own*(f[:,b]-.5)*width[:,b]*d)
  da=derivative(a);factor=(-xi*width[:,a]/width[:,b])if normal_derivative else(-.5*width[:,a]**2/width[:,b]*(xi*xi-.25))
  for plane in(0,1):
   for shift,d in zip((-1,0,1),da):
    delta=[0,0,0];delta[b]=plane;delta[a]=shift;add(b,delta,factor*(1 if plane else-1)*d)
 assert len(terms)==22;indices=np.empty((K,24),dtype='i8');weights=np.zeros((K,24))
 for column,((comp,delta),value)in enumerate(terms.items()):
  q=cell+np.asarray(delta);shape=shapes[comp];assert all(np.all((q[:,b]>=0)&(q[:,b]<shape[b]))for b in range(3));indices[:,column]=off[comp]+(q[:,0]*shape[1]+q[:,1])*shape[2]+q[:,2];weights[:,column]=value
 indices[:,22:]=indices[:,:1];assert np.isfinite(weights).all();return indices,weights

def coordinates(faces,indices):
 n,shapes,off,cent=layout(faces);flat=indices.ravel();components=np.searchsorted(off[1:],flat,side='right');position=np.empty((len(flat),3))
 for a in range(3):
  mask=components==a;q=flat[mask]-off[a];shape=shapes[a];ijk=[q//(shape[1]*shape[2]),q//shape[2]%shape[1],q%shape[2]]
  for b in range(3):position[mask,b]=(faces[b]if a==b else cent[b])[ijk[b]]
 return components.reshape(indices.shape),position.reshape(indices.shape+(3,))

def sample(faces,velocity,points):
 return np.column_stack([np.sum(w*velocity[ix],axis=1)for ix,w in[rows(faces,points,a)for a in range(3)]])
