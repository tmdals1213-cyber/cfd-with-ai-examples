"""Dimensional runner-frame kinematics. No fluid solver or surface repair.

The finite mode uses all43 original samples. The separate periodic C1 mode
passes through original samples0..41 and constructs the last interval to0;
source endpoint42 is retained for provenance but excluded from that loop.
"""
import math
import numpy as np

def running_speed(cadence_spm,step_length_m):
 if not all(math.isfinite(float(x)) and x>0 for x in (cadence_spm,step_length_m)):raise ValueError('positive cadence steps/min and one-step length m required')
 return cadence_spm*step_length_m/60.

class Motion:
 def __init__(self,cache,cadence_spm,steps_per_cycle=2,mode='periodic_c1'):
  if mode not in ('periodic_c1','finite_linear'):raise ValueError('explicit motion mode')
  if cadence_spm<=0 or steps_per_cycle<=0:raise ValueError('positive cadence and steps/cycle')
  with np.load(cache,allow_pickle=False) as z:
   self.p=z['positions_m'].astype('f8');self.faces=z['triangles'].copy();self.times=z['source_times_s'].astype('f8')
  if len(self.p)!=43 or self.p.shape[2]!=3 or len(self.times)!=43:raise ValueError('qualified43sample layout')
  self.dt=1/60.;self.N=42;self.T_source=self.N*self.dt
  if not np.allclose(self.times,np.arange(43)/60,atol=1e-11,rtol=0):raise ValueError('source timestamps')
  self.cadence=float(cadence_spm);self.steps_per_cycle=steps_per_cycle;self.T=60.*steps_per_cycle/cadence_spm;self.rate=self.T_source/self.T;self.mode=mode
  q=self.p[:self.N];self.tangent=(np.roll(q,-1,axis=0)-np.roll(q,1,axis=0))/(2*self.dt)
 def sample(self,time_s,phase_cycles=0.,offset_m=(0.,0.,0.)):
  if not math.isfinite(float(time_s)) or not math.isfinite(float(phase_cycles)):raise ValueError('finite time/phase')
  source_t=(time_s/self.T+phase_cycles)*self.T_source
  if self.mode=='finite_linear':
   if not 0<=source_t<=self.T_source:raise ValueError('finite raw clip ends; no implicit wrap')
   if source_t==0.:return self.p[0].copy()+np.array(offset_m),(self.p[1]-self.p[0])*(self.rate/self.dt)
   if source_t==self.T_source:return self.p[-1].copy()+np.array(offset_m),(self.p[-1]-self.p[-2])*(self.rate/self.dt)
   j=min(int(source_t/self.dt),self.N-1);s=(source_t-j*self.dt)/self.dt;p=(1-s)*self.p[j]+s*self.p[j+1];v=(self.p[j+1]-self.p[j])*(self.rate/self.dt)
  else:
   a=(source_t%self.T_source)/self.dt;j=int(a)%self.N;s=a-math.floor(a);k=(j+1)%self.N;p0,p1=self.p[j],self.p[k];m0,m1=self.tangent[j],self.tangent[k]
   p=(2*s**3-3*s*s+1)*p0+(s**3-2*s*s+s)*self.dt*m0+(-2*s**3+3*s*s)*p1+(s**3-s*s)*self.dt*m1
   v=((6*s*s-6*s)*p0/self.dt+(3*s*s-4*s+1)*m0+(-6*s*s+6*s)*p1/self.dt+(3*s*s-2*s)*m1)*self.rate
  return p+np.array(offset_m),v
 def exact_position_envelope(self):
  # Coordinate extrema of EACH cubic interval, including derivative roots.
  lo=np.full(3,np.inf);hi=-lo;velocity_component_bound=0.
  for j in range(self.N):
   k=(j+1)%self.N;p0,p1=self.p[j],self.p[k];m0=self.dt*self.tangent[j];m1=self.dt*self.tangent[k]
   A=2*p0-2*p1+m0+m1;B=-3*p0+3*p1-2*m0-m1;C=m0;D=p0
   lo=np.minimum(lo,np.minimum(p0,p1).min(0));hi=np.maximum(hi,np.maximum(p0,p1).max(0))
   # Handle linear derivative and quadratic derivative without state clipping.
   a=3*A;b=2*B;c=C;disc=b*b-4*a*c;roots=[]
   linear=(abs(a)<1e-14)&(abs(b)>1e-14);r=np.full(a.shape,np.nan);np.divide(-c,b,out=r,where=linear);roots.append(r)
   quad=(abs(a)>=1e-14)&(disc>=0)
   for sign in (-1,1):
    r=np.full(a.shape,np.nan);np.divide(-b+sign*np.sqrt(np.maximum(0.,disc)),2*a,out=r,where=quad);roots.append(r)
   for s in roots:
    valid=(s>0)&(s<1);value=((A*s+B)*s+C)*s+D
    for axis in range(3):
     z=value[:,axis][valid[:,axis]]
     if len(z):lo[axis]=min(lo[axis],float(z.min()));hi[axis]=max(hi[axis],float(z.max()))
   velocity_component_bound=max(velocity_component_bound,float(np.max((abs(3*A)+abs(2*B)+abs(C))*(self.rate/self.dt))))
  return lo,hi,velocity_component_bound

def relative_air_and_ground(cadence_spm,step_length_m,air_lab_m_s=(0.,0.,0.)):
 V=running_speed(cadence_spm,step_length_m);runner=np.array([V,0.,0.]);return np.array(air_lab_m_s)-runner,-runner
