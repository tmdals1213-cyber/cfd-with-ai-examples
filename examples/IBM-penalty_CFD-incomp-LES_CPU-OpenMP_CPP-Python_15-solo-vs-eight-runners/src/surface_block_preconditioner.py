"""SPD pressure/body block factorization; preconditioner only, not a physical operator."""
def make_blocks(np,faces,points,person,indices,weights,md,dt,ks,area):
 K=len(indices);assert K==3*len(points)and K<=49152;cells=np.column_stack([np.searchsorted(faces[a],points[:,a],side='right')-1 for a in range(3)]);groups={}
 for point in range(len(points)):
  for axis in range(3):groups.setdefault((int(person[point]),)+tuple(int(v)for v in cells[point])+(axis,),[]).append(3*point+axis)
 batches={};epsilon=1/(dt*ks*area)
 for key,allrows in sorted(groups.items()):
  for start in range(0,len(allrows),48):
   rows=np.array(allrows[start:start+48],dtype='i4');unique=np.unique(indices[rows].ravel());A=np.zeros((len(rows),len(unique)));loc=np.searchsorted(unique,indices[rows]);np.add.at(A,(np.repeat(np.arange(len(rows)),indices.shape[1]),loc.ravel()),weights[rows].ravel());Z=A/np.sqrt(md[unique])[None,:];gram=Z.dot(Z.T);gram[np.diag_indices(len(rows))]+=epsilon[rows];factor=np.linalg.cholesky(gram);assert np.isfinite(factor).all()and np.all(np.diag(factor)>0);batches.setdefault(len(rows),[]).append((rows,factor))
 out={}
 for size,items in batches.items():out[size]=(np.stack([v[0]for v in items]),np.stack([v[1]for v in items]))
 allrows=np.concatenate([v[0].ravel()for v in out.values()]);assert np.array_equal(np.sort(allrows),np.arange(K));return dict(K=K,batches=out)
def solve_blocks(np,blocks,rhs):
 result=np.zeros(blocks['K'])
 for size,(rows,L)in blocks['batches'].items():
  z=np.linalg.solve(L,np.asarray(rhs,dtype='f8')[rows][...,None]);z=np.linalg.solve(np.swapaxes(L,1,2),z)[...,0];result[rows]=z
 return result
def save_blocks(np,blocks,path):
 data={'K':np.array([blocks['K']]),'sizes':np.array(sorted(blocks['batches']),dtype='i4')}
 for size,(rows,L)in blocks['batches'].items():data['rows%d'%size]=rows;data['factor%d'%size]=L
 with path.open('xb')as stream:np.savez(stream,**data)
def load_blocks(np,archive):
 K=int(archive['K'][0]);batches={}
 for size in archive['sizes']:
  size=int(size);rows=archive['rows%d'%size].copy();L=archive['factor%d'%size].copy();assert 1<=size<=48 and rows.ndim==2 and rows.shape[1]==size and L.shape==rows.shape+(size,)and np.isfinite(L).all()and np.all(np.diagonal(L,axis1=1,axis2=2)>0)and np.all(np.triu(L,1)==0);batches[size]=(rows,L)
 allrows=np.concatenate([v[0].ravel()for v in batches.values()]);assert np.array_equal(np.sort(allrows),np.arange(K));return dict(K=K,batches=batches)
def make_action(np,pressure_pc,blocks,n,shapes,off,areas,md,indices,weights):
 from tensor_pressure_pc import apply as pressure_apply
 D=len(md);N=int(np.prod(n));K=len(indices)
 def gradient(pressure):
  out=np.zeros(D);pr=np.asarray(pressure,dtype='f8').reshape(n)
  for a in range(3):
   mid=[slice(None)]*3;mid[a]=slice(1,-1);left=[slice(None)]*3;right=left.copy();left[a]=slice(None,-1);right[a]=slice(1,None);out[off[a]:off[a+1]].reshape(shapes[a])[tuple(mid)]=areas[a][tuple(mid)]*(pr[tuple(right)]-pr[tuple(left)])
  out[:off[1]].reshape(shapes[0])[0]=areas[0][0]*pr[0];return out
 def divergence(value):
  out=np.zeros(n)
  for a in range(3):out-=np.diff(areas[a]*value[off[a]:off[a+1]].reshape(shapes[a]),axis=a)
  return out.ravel()
 def transpose_C(pressure):
  velocity=gradient(pressure)/md;out=np.zeros(K)
  for k in range(indices.shape[1]):out+=weights[:,k]*velocity[indices[:,k]]
  return out
 def C_action(reaction):
  spread=np.zeros(D);np.add.at(spread,indices.ravel(),(weights*np.asarray(reaction,dtype='f8')[:,None]).ravel());return divergence(spread/md)
 def apply(rhs):
  rhs=np.asarray(rhs,dtype='f8');assert rhs.shape==(N+K,);first=pressure_apply(np,pressure_pc,rhs[:N]);reduced=rhs[N:]-transpose_C(first);last=solve_blocks(np,blocks,reduced);out=np.r_[first-pressure_apply(np,pressure_pc,C_action(last)),last];assert np.isfinite(out).all();return out
 return apply,C_action,transpose_C
