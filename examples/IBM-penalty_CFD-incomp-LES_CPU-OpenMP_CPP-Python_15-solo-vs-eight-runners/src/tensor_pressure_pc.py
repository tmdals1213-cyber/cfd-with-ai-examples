"""SPD tensor inverse preconditioner for constant-mobility MAC pressure blocks.

This is a preconditioner. Actual solver acceptance still requires the original
true operator, strict raw residual, BOTHdiv, momentum/work and max384 gates.
"""
def assemble(np,widths,rho,bases,eigenvalues):
 shape=tuple(len(h)for h in widths);den=eigenvalues[0][:,None,None]+eigenvalues[1][None,:,None]+eigenvalues[2][None,None,:];assert np.isfinite(den).all()and float(den.min())>0
 scale=1/np.sqrt(widths[0][:,None,None]*widths[1][None,:,None]*widths[2][None,None,:]);return dict(shape=shape,widths=widths,rho=rho,bases=bases,eigenvalues=eigenvalues,den=den,scale=scale)

def make(np,faces,rho):
 assert rho>0;hs=[];bases=[];eigenvalues=[]
 for a in range(3):
  h=np.diff(np.asarray(faces[a],dtype='f8'));assert len(h)>=2 and np.all(h>0);n=len(h);K=np.zeros((n,n))
  for j in range(1,n):
   w=1/(.5*(h[j-1]+h[j]));K[j-1,j-1]+=w;K[j,j]+=w;K[j-1,j]-=w;K[j,j-1]-=w
  if a==0:K[0,0]+=2/h[0]
  inv=1/np.sqrt(h);symmetric=inv[:,None]*K*inv[None,:];lam,V=np.linalg.eigh(symmetric);assert float(abs(V.T.dot(V)-np.eye(n)).max())<=2e-12
  hs.append(h);bases.append(V);eigenvalues.append(lam)
 return assemble(np,hs,rho,bases,eigenvalues)

def apply(np,pc,residual):
 residual=np.asarray(residual,dtype='f8');assert residual.size==int(np.prod(pc['shape']))and np.isfinite(residual).all();value=residual.reshape(pc['shape'])*pc['scale']
 for a in range(3):value=np.moveaxis(np.tensordot(pc['bases'][a].T,value,axes=(1,a)),0,a)
 value=value/pc['den']
 for a in(2,1,0):value=np.moveaxis(np.tensordot(pc['bases'][a],value,axes=(1,a)),0,a)
 out=pc['rho']*pc['scale']*value;assert np.isfinite(out).all();return np.ascontiguousarray(out).ravel()

def save(np,pc,path):
 data={'rho':np.array([pc['rho']])}
 for a in range(3):data['width%d'%a]=pc['widths'][a];data['basis%d'%a]=pc['bases'][a];data['eigen%d'%a]=pc['eigenvalues'][a]
 with path.open('xb')as stream:np.savez(stream,**data)

def load(np,archive):
 hs=[archive['width%d'%a].copy()for a in range(3)];basis=[archive['basis%d'%a].copy()for a in range(3)];eig=[archive['eigen%d'%a].copy()for a in range(3)];rho=float(archive['rho'][0]);assert rho>0
 for a in range(3):
  n=len(hs[a]);assert np.all(hs[a]>0)and basis[a].shape==(n,n)and eig[a].shape==(n,)and np.isfinite(basis[a]).all()and np.isfinite(eig[a]).all()and float(abs(basis[a].T.dot(basis[a])-np.eye(n)).max())<=2e-12
 return assemble(np,hs,rho,basis,eig)
