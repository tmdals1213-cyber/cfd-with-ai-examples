#if !defined(_WIN32)
#ifndef __declspec
#define __declspec(x) __attribute__((visibility("default")))
#endif
#endif
// Exploratory policy c0d43aea15ba2212c5c9e95dcc99bac7ea8ed3ef444df3ac4411c243a839943f; stability targets, no total-drag qualification.
#include <algorithm>
#include <omp.h>
#include <atomic>
#include <cstdint>
namespace parallel_evidence {
std::atomic<unsigned long long> masks[3];
inline void note(int role){
 thread_local unsigned seen=0;
 unsigned bit=1u<<role;
 if(!(seen&bit)){int id=omp_get_thread_num();if(id<64)masks[role].fetch_or(1ULL<<id,std::memory_order_relaxed);seen|=bit;}
}

}

#include <cmath>
#include <cstddef>
#include <vector>
#include <array>
#include <limits>
#include <cstdio>
#ifndef NX
#define NX 5
#define NY 7
#define NZ 9
#endif
namespace parallel_evidence {
inline double fixed_dot(const std::vector<double>&a,const std::vector<double>&b){
 const size_t block=4096,count=(a.size()+block-1)/block;
 std::vector<double> partial(count,0.);
 #pragma omp parallel for schedule(static)
 for(size_t k=0;k<count;k++){
  note(2);double v=0;size_t stop=std::min(a.size(),(k+1)*block);
  for(size_t i=k*block;i<stop;i++)v+=a[i]*b[i];partial[k]=v;
 }
 double result=0;for(double value:partial)result+=value;return result;
}
}
namespace pacemaker_air_v1 {
const int n[3]={NX,NY,NZ};
const size_t N=size_t(NX)*NY*NZ;
static size_t ix(const int* s,const int* t){return (size_t(t[0])*s[1]+t[1])*s[2]+t[2];}
struct Node {int a,t[3];double V,M;};
struct Edge {size_t z,l,r;double area,dist,w;};
struct Grid {
 std::vector<double> h[3];int ds[3][3];size_t off[3],D;std::vector<Node> nodes;
 double area(int a,const int* t)const{double v=1.;for(int b=0;b<3;b++)if(a!=b)v*=h[b][t[b]];return v;}
 double span(int a,int i)const{return i==0?.5*h[a][0]:i==n[a]?.5*h[a].back():.5*(h[a][i-1]+h[a][i]);}
 bool boundary(const Node& q)const{return q.t[q.a]==0||q.t[q.a]==n[q.a];}
 bool closed(const Node& q)const{return q.a!=0&&boundary(q);}
 bool inlet(const Node& q)const{return q.a==0&&q.t[0]==n[0];}
 bool outlet(const Node& q)const{return q.a==0&&q.t[0]==0;}
 double dualArea(const Node& q,int b)const{double v=1.;for(int c=0;c<3;c++)if(c!=b)v*=c==q.a?span(c,q.t[c]):h[c][q.t[c]];return v;}
 double Q(int a,const int* t,const std::vector<double>& u)const{return area(a,t)*u[off[a]+ix(ds[a],t)];}
 double lifted(int a,int b,const int* t,const std::vector<double>& u)const{
  int l[3]={t[0],t[1],t[2]},r[3]={t[0],t[1],t[2]};
  if(a==b){int j=t[a];if(j==0)return Q(b,t,u);if(j==n[a]+1){l[a]=n[a];return Q(b,l,u);}l[a]=j-1;return .5*(Q(b,l,u)+Q(b,r,u));}
  double v=0.;if(t[a]>0){l[a]--;v+=.5*Q(b,l,u);}if(t[a]<n[a])v+=.5*Q(b,r,u);return v;
 }
};
static double dot(const std::vector<double>& a,const std::vector<double>& b){double z=0.;for(size_t i=0;i<a.size();i++)z+=(double)a[i]*b[i];return z;}
static void apply(const std::vector<Edge>& e,const std::vector<double>& p,std::vector<double>& v){std::fill(v.begin(),v.end(),0.);for(const Edge& f:e){if(f.l==N){v[f.r]+=f.w*p[f.r];}else{double q=f.w*(p[f.l]-p[f.r]);v[f.l]+=q;v[f.r]-=q;}}}
static bool finite(const double* v,size_t n){for(size_t j=0;j<n;j++)if(!std::isfinite(v[j]))return false;return true;}
}

namespace pacemaker_air_v1 {
static bool build(Grid&g,const double*xf,const double*yf,const double*zf,double rho){
 const double*f[3]={xf,yf,zf};g.D=0;
 for(int a=0;a<3;a++){for(int j=0;j<=n[a];j++)if(!std::isfinite(f[a][j])||(j&&f[a][j]<=f[a][j-1]))return false;for(int j=0;j<n[a];j++)g.h[a].push_back(f[a][j+1]-f[a][j]);}
 for(int a=0;a<3;a++){for(int b=0;b<3;b++)g.ds[a][b]=n[b];g.ds[a][a]++;g.off[a]=g.D;g.D+=size_t(g.ds[a][0])*g.ds[a][1]*g.ds[a][2];
  for(int i=0;i<g.ds[a][0];i++)for(int j=0;j<g.ds[a][1];j++)for(int k=0;k<g.ds[a][2];k++){Node q;q.a=a;q.t[0]=i;q.t[1]=j;q.t[2]=k;q.V=g.area(a,q.t)*g.span(a,q.t[a]);q.M=rho*q.V;if(!(q.M>0)||!std::isfinite(q.M))return false;g.nodes.push_back(q);}}
 return true;
}
static bool divs(const Grid&g,const double*P,double&dc,double&dv){
 dc=dv=0;for(int i=0;i<n[0];i++)for(int j=0;j<n[1];j++)for(int k=0;k<n[2];k++){int t[3]={i,j,k};double c=0,f=0,V=g.h[0][i]*g.h[1][j]*g.h[2][k];for(int a=0;a<3;a++){int r[3]={i,j,k};r[a]++;size_t l=g.off[a]+ix(g.ds[a],t),q=g.off[a]+ix(g.ds[a],r);double lo=P[l]/g.nodes[l].M,hi=P[q]/g.nodes[q].M,ar=g.area(a,t);c+=(hi-lo)/g.h[a][t[a]];f+=ar*hi-ar*lo;}if(!std::isfinite(c)||!std::isfinite(f))return false;dc=std::max(dc,std::abs(c));dv=std::max(dv,std::abs(f/V));}return true;
}
using LD=double;using VV=std::vector<LD>;
static bool qfinite(LD x){return __builtin_isfinite(x);}
static LD qabs(LD x){return x<0?-x:x;}
static LD qsqrt(LD x){if(x==0)return 0;if(!(x>0)||!qfinite(x))return (LD)std::numeric_limits<double>::quiet_NaN();LD y=(LD)std::sqrt((double)x);for(int k=0;k<2;k++)y=(y+x/y)/(LD)2;return y;}
static LD dotld(const VV&a,const VV&b){LD z=0;for(size_t i=0;i<N;i++)z+=a[i]*b[i];return z;}
static void applyld(const std::vector<Edge>&e,const VV&p,VV&v){std::fill(v.begin(),v.end(),0);for(const Edge&f:e){if(f.l==N)v[f.r]+=(LD)f.w*p[f.r];else{LD q=(LD)f.w*(p[f.l]-p[f.r]);v[f.l]+=q;v[f.r]-=q;}}}
}


#include "vreman_sgs.hpp"
extern "C" __declspec(dllexport) int pm_dimensions(int*out){if(!out)return 1;out[0]=NX;out[1]=NY;out[2]=NZ;return 0;}
// Stage 1 material transport and penalty predictor ONLY. Actual Md/Pstar/b return
// before the distinct pressure API is called, allowing SHA-bound RHS admission.
extern "C" __declspec(dllexport) int pm_prepare(const double*xf,const double*yf,const double*zf,double rho,double dt,double kappa,double nu,double wind,
 const double*oldP,const double*wall,const double*chi,double*out,double*stats){
 using namespace pacemaker_air_v1;
 if(!xf||!yf||!zf||!oldP||!wall||!chi||!out||!stats)return 1;
 if(rho!=1.225||!(dt>0)||kappa<0||nu<0||!(wind<0)||!std::isfinite(dt+kappa+nu+wind))return 2;
 Grid g;if(!build(g,xf,yf,zf,rho)||!finite(oldP,g.D)||!finite(wall,g.D)||!finite(chi,g.D))return 3;
 double dc,dv,hmin=1e99;for(int a=0;a<3;a++)for(double h:g.h[a])hmin=std::min(hmin,h);if(!divs(g,oldP,dc,dv)||std::max(dc,dv)*hmin/std::abs(wind)>1e-4||dt*std::max(dc,dv)>1e-4)return 4;
 std::vector<double>u(g.D),adv(g.D),vis(g.D),pr(g.D),pt(g.D),sr(g.D),ss(g.D),jt(g.D),js(g.D),md(g.D),mob(g.D),b(N,0.);
 for(size_t z=0;z<g.D;z++){const Node&q=g.nodes[z];if(chi[z]<0||chi[z]>1||(g.boundary(q)&&chi[z]!=0)||(g.closed(q)&&oldP[z]!=0)||(g.inlet(q)&&oldP[z]!=q.M*wind)||(g.outlet(q)&&oldP[z]>0))return 4;u[z]=oldP[z]/q.M;md[z]=q.M;mob[z]=1/(rho+dt*kappa*chi[z]);}
 double cflMax=0,diffMax=0;
 int transportError=0;
 #pragma omp parallel for schedule(static) reduction(max:cflMax,diffMax,transportError)
 for(size_t z=0;z<g.D;z++){parallel_evidence::note(0);const Node&q=g.nodes[z];int a=q.a;double outgoing=0,diag=0;
  for(int ax=0;ax<3;ax++){int lo[3]={q.t[0],q.t[1],q.t[2]},hi[3]={q.t[0],q.t[1],q.t[2]},upper[3]={q.t[0],q.t[1],q.t[2]};lo[ax]--;hi[ax]++;upper[ax]++;
   double gl=g.lifted(a,ax,q.t,u),gh=g.lifted(a,ax,upper,u);bool hl=q.t[ax]>0,hr=q.t[ax]+1<g.ds[a][ax];double ul=hl?u[g.off[a]+ix(g.ds[a],lo)]:0,ur=hr?u[g.off[a]+ix(g.ds[a],hi)]:0;
   if(ax==0&&!hl)ul=u[z];if(ax==0&&!hr)ur=a==0?wind:0;
   adv[z]-=(gh*(gh>=0?u[z]:ur)-gl*(gl>=0?ul:u[z]))/q.V;outgoing+=std::max(gh,0.)+std::max(-gl,0.);double ar=g.dualArea(q,ax);
   if(hl){double dis=ax==a?g.h[ax][q.t[ax]-1]:.5*(g.h[ax][q.t[ax]-1]+g.h[ax][q.t[ax]]);double c=ar/dis;vis[z]+=nu*c*(ul-u[z])/q.V;diag+=c;}
   if(hr){double dis=ax==a?g.h[ax][q.t[ax]]:.5*(g.h[ax][q.t[ax]]+g.h[ax][q.t[ax]+1]);double c=ar/dis;vis[z]+=nu*c*(ur-u[z])/q.V;diag+=c;}}
  double cf=dt*outgoing/q.V,df=dt*nu*diag/q.V;cflMax=std::max(cflMax,cf);diffMax=std::max(diffMax,df);if(!std::isfinite(adv[z]+vis[z])||cf+df>.5)transportError=5;
  pr[z]=oldP[z]+dt*q.M*(adv[z]+vis[z]);jt[z]=g.closed(q)?-pr[z]:g.inlet(q)?q.M*wind-pr[z]:0.;pt[z]=g.closed(q)?0:g.inlet(q)?q.M*wind:pr[z];
  sr[z]=(pt[z]+dt*q.V*kappa*chi[z]*wall[z])*(rho/(rho+dt*kappa*chi[z]));js[z]=g.closed(q)?-sr[z]:g.inlet(q)?q.M*wind-sr[z]:0.;ss[z]=g.closed(q)?0:g.inlet(q)?q.M*wind:sr[z];if(!std::isfinite(ss[z]))transportError=5;}
 if(transportError)return transportError;
 #pragma omp parallel for collapse(3) schedule(static)
 for(int i=0;i<NX;i++)for(int j=0;j<NY;j++)for(int k=0;k<NZ;k++){int t[3]={i,j,k};double sum=0;for(int a=0;a<3;a++){int r[3]={i,j,k};r[a]++;size_t l=g.off[a]+ix(g.ds[a],t),q=g.off[a]+ix(g.ds[a],r);double ar=g.area(a,t);double qlo=(ss[l]/md[l])*ar,qhi=(ss[q]/md[q])*ar;sum+=qhi-qlo;}b[ix(n,t)]=-sum/dt;}
 const std::vector<double>*v[]={&adv,&vis,&pr,&pt,&sr,&ss,&jt,&js,&md,&mob};for(int i=0;i<10;i++)std::copy(v[i]->begin(),v[i]->end(),out+size_t(i)*g.D);std::copy(b.begin(),b.end(),out+10*g.D);stats[0]=cflMax;stats[1]=diffMax;stats[2]=dc;stats[3]=dv;stats[4]=double(N);stats[5]=double(g.D);return 0;
}
// Variational surface coupling. This is an explicitly new physical operator.
namespace surface_coupled_v1 {
using namespace pacemaker_air_v1;
using PCFunction=int(*)(const double*,double*,void*);
static LD scalar(const VV&a,const VV&b){return parallel_evidence::fixed_dot(a,b);}
struct Operator {
 Grid g;bool valid;int K;size_t L;const double*md;const long long*index;const double*weights;std::vector<Edge>edges;VV epsilon,gp,ib,force;
 std::vector<size_t>rowptr;std::vector<int>surface_row;std::vector<double>surface_weight;
 Operator(const double*xf,const double*yf,const double*zf,double rho,double dt,double penalty,const double*mass,int count,const long long*ixs,const double*ws,const double*areas):valid(false),K(count),md(mass),index(ixs),weights(ws){
  if(!mass||!ixs||!ws||!areas||K<3||K%3||K>49152||!(dt>0)||!(penalty>0)||rho!=1.225||!qfinite(dt+penalty)||!build(g,xf,yf,zf,rho))return;
  L=N+size_t(K);epsilon.resize(K);gp.assign(g.D,0);ib.assign(g.D,0);force.assign(g.D,0);
  for(size_t i=0;i<g.D;i++)if(mass[i]!=g.nodes[i].M)return;
  for(int r=0;r<K;r++){if(!(areas[r]>0)||!qfinite(areas[r]))return;if(areas[r]!=areas[(r/3)*3])return;epsilon[r]=1/((LD)dt*penalty*areas[r]);LD sum=0;LD componentSums[3]={0,0,0};
   for(int k=0;k<24;k++){size_t t=size_t(r)*24+k;long long q=ixs[t];double w=ws[t];if(q<0||size_t(q)>=g.D||!qfinite(w))return;if(w!=0&&(g.closed(g.nodes[q])||g.inlet(g.nodes[q])))return;sum+=w;componentSums[g.nodes[q].a]+=w;}
   if(qabs(sum-1)>4e-15L)return;for(int a=0;a<3;a++)if(qabs(componentSums[a]-(a==r%3?1:0))>4e-15L)return;}
  for(size_t v=0;v<g.D;v++){const Node&q=g.nodes[v];int a=q.a,j=q.t[a];if(g.outlet(q)){int r[3]={0,q.t[1],q.t[2]};Edge e={v,N,ix(n,r),g.area(a,q.t),g.span(a,0),0};edges.push_back(e);continue;}if(g.inlet(q)||j==0||j==n[a])continue;int l[3]={q.t[0],q.t[1],q.t[2]},r[3]={q.t[0],q.t[1],q.t[2]};l[a]--;Edge e={v,ix(n,l),ix(n,r),g.area(a,q.t),g.span(a,j),0};edges.push_back(e);}
  rowptr.assign(g.D+1,0);
  for(int r=0;r<K;r++)for(int k=0;k<24;k++){size_t t=size_t(r)*24+k;if(weights[t]!=0)rowptr[size_t(index[t])+1]++;}
  for(size_t i=0;i<g.D;i++)rowptr[i+1]+=rowptr[i];
  surface_row.resize(rowptr.back());surface_weight.resize(rowptr.back());
  std::vector<size_t>cursor=rowptr;
  for(int r=0;r<K;r++)for(int k=0;k<24;k++){
   size_t t=size_t(r)*24+k;if(weights[t]==0)continue;
   size_t slot=cursor[size_t(index[t])]++;surface_row[slot]=r;surface_weight[slot]=weights[t];
  }
  valid=true;
 }
 void pressure_rows(const VV&v,VV&out,const double*base=nullptr,double sign=1.)const{
  #pragma omp parallel for schedule(static)
  for(size_t cell=0;cell<N;cell++){
   parallel_evidence::note(1);
   int t[3]={int(cell/(NY*NZ)),int((cell/NZ)%NY),int(cell%NZ)};
   double result=base?base[cell]:0.;
   for(int a=0;a<3;a++){
    double ar=g.area(a,t);
    if(a==0||t[a]>0){size_t face=g.off[a]+ix(g.ds[a],t);result+=sign*ar*(v[face]/md[face]);}
    if(t[a]+1<n[a]){int upper[3]={t[0],t[1],t[2]};upper[a]++;size_t face=g.off[a]+ix(g.ds[a],upper);result-=sign*ar*(v[face]/md[face]);}
   }
   out[cell]=result;
  }
 }
 void apply(const VV&x,VV&y){
  std::fill(gp.begin(),gp.end(),0.);
  #pragma omp parallel for schedule(static)
  for(size_t eidx=0;eidx<edges.size();eidx++){
   parallel_evidence::note(1);const Edge&e=edges[eidx];gp[e.z]=(LD)e.area*(x[e.r]-(e.l==N?0:x[e.l]));
  }
  #pragma omp parallel for schedule(static)
  for(size_t face=0;face<g.D;face++){
   parallel_evidence::note(1);double value=0.;
   for(size_t t=rowptr[face];t<rowptr[face+1];t++)value+=surface_weight[t]*x[N+surface_row[t]];
   ib[face]=value;force[face]=gp[face]+value;
  }
  pressure_rows(force,y);
  #pragma omp parallel for schedule(static)
  for(int r=0;r<K;r++){
   parallel_evidence::note(1);LD v=0.;
   for(int k=0;k<24;k++){size_t t=size_t(r)*24+k;v+=(LD)weights[t]*(force[index[t]]/(LD)md[index[t]]);}
   y[N+r]=v+epsilon[r]*x[N+r];
  }
 }
 void effective(const double*b,VV&rhs){pressure_rows(ib,rhs,b,-1.);}

};
}
extern "C" __declspec(dllexport) int pm_surface_probe(const double*xf,const double*yf,const double*zf,double rho,double dt,double penalty,const double*md,int K,const long long*indices,const double*weights,const double*areas,const double*x,double*y){
 using namespace surface_coupled_v1;Operator op(xf,yf,zf,rho,dt,penalty,md,K,indices,weights,areas);if(!op.valid||!x||!y)return 1;VV a(op.L),b(op.L);for(size_t i=0;i<op.L;i++){if(!qfinite(x[i]))return 2;a[i]=x[i];}op.apply(a,b);std::copy(b.begin(),b.end(),y);return 0;
}
extern "C" __declspec(dllexport) int pm_surface_quad_project_slice(const double*xf,const double*yf,const double*zf,double rho,double dt,double penalty,const double*prepared,int K,const long long*indices,const double*weights,const double*areas,const double*wall,const double*rhs,double*out,double*trace,double*stats,double*state,int maxAdvance,surface_coupled_v1::PCFunction precondition,void*context){
 using namespace surface_coupled_v1;
 if(!xf||!yf||!zf||!prepared||!wall||!rhs||!out||!trace||!stats||!state||!precondition||K<3||K%3||K>49152||maxAdvance<1||maxAdvance>64)return 1;
 Grid check;if(!build(check,xf,yf,zf,rho))return 1;size_t D=check.D;const double*ss=prepared+5*D;const double*md=prepared+8*D;const double*mob=prepared+9*D;
 if(!finite(prepared,10*D+N)||!finite(wall,size_t(K))||!finite(rhs,N+size_t(K)))return 2;
 for(size_t i=0;i<D;i++)if(mob[i]!=1/rho)return 2;
 Operator op(xf,yf,zf,rho,dt,penalty,md,K,indices,weights,areas);if(!op.valid)return 2;size_t L=op.L;
 for(size_t i=0;i<N;i++)if(rhs[i]!=prepared[10*D+i])return 2;
 for(int r=0;r<K;r++){double value=0;for(int k=0;k<24;k++){size_t t=size_t(r)*24+k;value+=weights[t]*(ss[indices[t]]/md[indices[t]]);}if(rhs[N+r]!=(value-wall[r])/dt)return 2;}
 VV x(L,0),res(L),direction(L),z(L),Ax(L),raw(L),effective(N),original(L);for(size_t i=0;i<L;i++)original[i]=rhs[i];
 LD bn=qsqrt(scalar(original,original)),tau=std::max((LD)1e-10,(LD)1e-5*bn),rz=0;int iteration=0;double hmin=1e99;for(int a=0;a<3;a++)for(double h:check.h[a])hmin=std::min(hmin,h);const double divLimit=1e-4*(2813./480.)/hmin;
 auto pc=[&](const VV&a,VV&b){std::vector<double>bridgeIn(L),bridgeOut(L);for(size_t i=0;i<L;i++)bridgeIn[i]=(double)a[i];if(precondition(bridgeIn.data(),bridgeOut.data(),context)!=0)return false;for(size_t i=0;i<L;i++)b[i]=(LD)bridgeOut[i];for(size_t i=0;i<L;i++)if(!qfinite(b[i]))return false;state[11]+=1;return true;};
 if(state[0]==0){for(size_t i=0;i<3*L+32;i++)if(state[i]!=0)return 2;res=original;if(!pc(res,z))return 3;direction=z;rz=scalar(res,z);}
 else{
  if(state[0]!=1794426759L||state[1]<0||state[1]>384||std::floor(state[1])!=state[1]||state[5]!=0||state[6]!=NX||state[7]!=NY||state[8]!=NZ||state[9]!=K||state[10]!=7||((LD*)(state+20))[1]!=bn||((LD*)(state+20))[2]!=tau)return 2;
  iteration=int(state[1]);rz=((LD*)(state+20))[0];const LD*stored=(const LD*)(state+32);for(size_t i=0;i<L;i++){x[i]=stored[i];res[i]=stored[L+i];direction[i]=stored[2*L+i];if(!qfinite(x[i])||!qfinite(res[i])||!qfinite(direction[i]))return 2;}}
 const int first=iteration;
 auto preserve=[&](int status){state[0]=1794426759L;state[1]=iteration;state[2]=rz;state[3]=bn;state[4]=tau;state[5]=status;state[6]=NX;state[7]=NY;state[8]=NZ;state[9]=K;state[10]=7;LD*scalars=(LD*)(state+20);scalars[0]=rz;scalars[1]=bn;scalars[2]=tau;LD*stored=(LD*)(state+32);for(size_t i=0;i<L;i++){stored[i]=x[i];stored[L+i]=res[i];stored[2*L+i]=direction[i];}};
 LD rn=bn,rp=bn,rs=bn,tp=0,en=0;double dc=0,dv=0;std::vector<double>pn(D),hi(L),lo(L);
 for(;;){
  LD recursive=qsqrt(scalar(res,res));if(!qfinite(recursive))return 6;bool checkTrue=iteration==first||iteration%16==0||iteration>=384||iteration-first>=maxAdvance||recursive<=tau;
  rn=rp=rs=tp=en=(LD)std::numeric_limits<double>::quiet_NaN();dc=dv=std::numeric_limits<double>::quiet_NaN();
  if(checkTrue){
   op.apply(x,Ax);state[12]+=1;LD sp=0,sk=0;for(size_t i=0;i<L;i++){raw[i]=original[i]-Ax[i];if(i<N)sp+=raw[i]*raw[i];else sk+=raw[i]*raw[i];}rn=qsqrt(sp+sk);rp=qsqrt(sp);rs=qsqrt(sk);if(!qfinite(rn))return 6;
   op.effective(rhs,effective);en=qsqrt(scalar(effective,effective));tp=std::max((LD)1e-10,(LD)1e-5*en);
   if(rn<=tau*.5L&&rp<=tp*.5L){VV candidate=x;for(size_t i=0;i<L;i++){double h=(double)x[i],l=(double)(x[i]-(LD)h);candidate[i]=(LD)h+(LD)l;}op.apply(candidate,Ax);LD cp=0,cs=0;for(size_t i=0;i<L;i++){LD r=original[i]-Ax[i];if(i<N)cp+=r*r;else cs+=r*r;}LD crn=qsqrt(cp+cs),crp=qsqrt(cp);if(crn<=tau*.5L&&crp<=tp*.5L){for(size_t i=0;i<D;i++)pn[i]=(double)((LD)ss[i]-(LD)dt*op.force[i]);if(!divs(op.g,pn.data(),dc,dv))return 6;if(dc<=divLimit&&dv<=divLimit&&dt*std::max(dc,dv)<=1e-4){x.swap(candidate);rn=crn;rp=crp;rs=qsqrt(cs);}else dc=dv=std::numeric_limits<double>::infinity();}else dc=dv=std::numeric_limits<double>::infinity();}
  }
  double*tr=trace+size_t(iteration)*9;tr[0]=iteration;tr[1]=(double)rn;tr[2]=(double)rp;tr[3]=(double)rs;tr[4]=(double)tau;tr[5]=(double)tp;tr[6]=dc;tr[7]=dv;tr[8]=(double)rz;
  if(iteration%64==0)std::fprintf(stderr,"SURFACE_CG iter=%d full=%.9e pressure=%.9e tau=%.9e\n",iteration,(double)rn,(double)rp,(double)tau);
  if(checkTrue&&rn<=tau*.5L&&rp<=tp*.5L&&dc<=divLimit&&dv<=divLimit&&dt*std::max(dc,dv)<=1e-4)break;
  if(iteration>=384||!(rz>0)||!qfinite(rz))return 6;if(iteration-first>=maxAdvance){preserve(0);return 101;}
  op.apply(direction,Ax);LD denominator=scalar(direction,Ax);if(!(denominator>0)||!qfinite(denominator))return 7;LD alpha=rz/denominator;
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<L;i++){parallel_evidence::note(2);x[i]+=alpha*direction[i];res[i]-=alpha*Ax[i];}if(!pc(res,z))return 3;LD rzNew=scalar(res,z);LD beta=rzNew/rz;
  #pragma omp parallel for schedule(static)
  for(size_t i=0;i<L;i++){parallel_evidence::note(2);direction[i]=z[i]+beta*direction[i];}
rz=rzNew;iteration++;
 }
 for(size_t i=0;i<D;i++)if(!qfinite(pn[i])||(op.g.outlet(op.g.nodes[i])&&pn[i]>0))return 8;
 for(size_t i=0;i<L;i++){hi[i]=(double)x[i];lo[i]=(double)(x[i]-(LD)hi[i]);if((LD)hi[i]+(LD)lo[i]!=x[i])return 9;}
 // op's cached forces are from the final true application, not a search vector.
 std::copy(pn.begin(),pn.end(),out);for(size_t i=0;i<D;i++){out[D+i]=(double)op.gp[i];out[2*D+i]=(double)op.ib[i];}std::copy(hi.begin(),hi.end(),out+3*D);std::copy(lo.begin(),lo.end(),out+3*D+L);
 preserve(1);stats[0]=iteration;stats[1]=(double)rn;stats[2]=(double)tau;stats[3]=(double)rp;stats[4]=(double)tp;stats[5]=dc;stats[6]=dv;stats[7]=(double)bn;stats[8]=(double)en;stats[9]=(double)rs;stats[10]=sizeof(LD);stats[11]=(double)state[11];return 0;
}

extern "C" __declspec(dllexport) int pm_surface_pc_probe(int count,const double*x,double*y,surface_coupled_v1::PCFunction pc,void*context){if(count<1||count>1500000||!x||!y||!pc)return 1;return pc(x,y,context);}

extern "C" __declspec(dllexport) int pm_quad_precision(int*out){using namespace pacemaker_air_v1;LD e=1;int bits=1;while((LD)1+e/(LD)2>(LD)1&&bits<200){e/=(LD)2;++bits;}out[0]=sizeof(LD);out[1]=bits;out[2]=1794426759;out[3]=7;return bits==53&&sizeof(LD)==8?0:1;}

// Diagnostic only: same outgoing dual-volume flux and diffusion rates as pm_prepare.
extern "C" __declspec(dllexport) int pm_flow_rates(const double*xf,const double*yf,const double*zf,double rho,double nu,const double*P,double*stats){
 using namespace pacemaker_air_v1;if(!xf||!yf||!zf||!P||!stats||rho!=1.225||nu<0)return 1;
 Grid g;if(!build(g,xf,yf,zf,rho)||!finite(P,g.D))return 2;std::vector<double>u(g.D);double vmax=0,hmin=1e99,ca=0,cd=0,cc=0,dc,dv;
 for(int a=0;a<3;a++)for(double h:g.h[a])hmin=std::min(hmin,h);
 for(size_t z=0;z<g.D;z++){u[z]=P[z]/g.nodes[z].M;vmax=std::max(vmax,std::abs(u[z]));}
 for(size_t z=0;z<g.D;z++){const Node&q=g.nodes[z];int a=q.a;double outgoing=0,diag=0;
  for(int ax=0;ax<3;ax++){int upper[3]={q.t[0],q.t[1],q.t[2]};upper[ax]++;double gl=g.lifted(a,ax,q.t,u),gh=g.lifted(a,ax,upper,u);outgoing+=std::max(gh,0.)+std::max(-gl,0.);double ar=g.dualArea(q,ax);
   if(q.t[ax]>0){double dis=ax==a?g.h[ax][q.t[ax]-1]:.5*(g.h[ax][q.t[ax]-1]+g.h[ax][q.t[ax]]);diag+=ar/dis;}
   if(q.t[ax]+1<g.ds[a][ax]){double dis=ax==a?g.h[ax][q.t[ax]]:.5*(g.h[ax][q.t[ax]]+g.h[ax][q.t[ax]+1]);diag+=ar/dis;}}
  double ra=outgoing/q.V,rd=nu*diag/q.V;ca=std::max(ca,ra);cd=std::max(cd,rd);cc=std::max(cc,ra+rd);}
 if(!divs(g,P,dc,dv))return 3;stats[0]=ca;stats[1]=cd;stats[2]=cc;stats[3]=vmax;stats[4]=hmin;stats[5]=dc;stats[6]=dv;return 0;
}

extern "C" __declspec(dllexport) void pm_parallel_report(unsigned long long*out){
 for(int i=0;i<3;i++)out[i]=parallel_evidence::masks[i].load();out[3]=omp_get_max_threads();
}

extern "C" __declspec(dllexport) int pm_sgs_probe(const double*xf,const double*yf,const double*zf,double rho,double cv,
 const double*momentum,const unsigned char*fluid,double*acceleration,double*nut,double*stats){
 using namespace penalty_les;
 Grid g;if(!build(g,xf,yf,zf,rho)||!momentum||!acceleration||!nut||!stats)return 1;
 std::vector<double>u(g.D),force;
 for(size_t i=0;i<g.D;i++)u[i]=momentum[i]/g.nodes[i].M;
 SGS s(g,u,fluid,cv);if(!s.valid)return 2;s.acceleration(force);
 double work=0.,bound=0.;
 for(size_t i=0;i<g.D;i++){work+=rho*g.nodes[i].V*u[i]*force[i];bound=std::max(bound,s.row_bound(g.nodes[i],s.maximum));}
 std::copy(force.begin(),force.end(),acceleration);std::copy(s.nut.begin(),s.nut.end(),nut);
 stats[0]=s.maximum;stats[1]=work;stats[2]=rho*s.dissipation();stats[3]=bound;return 0;
}
extern "C" __declspec(dllexport) double pm_vreman_gradient(const double*a,const double*h,double cv){
 double A[3][3];for(int i=0;i<3;i++)for(int j=0;j<3;j++)A[i][j]=a[i*3+j];return penalty_les::viscosity(A,h,cv);
}
extern "C" __declspec(dllexport) int pm_sgs_frozen_probe(const double*xf,const double*yf,const double*zf,double rho,
 const double*momentum,const double*nut,double*out,double*bound){
 using namespace penalty_les;Grid g;if(!build(g,xf,yf,zf,rho)||!momentum||!nut||!out||!bound)return 1;
 std::vector<double>u(g.D),force;for(size_t i=0;i<g.D;i++)u[i]=momentum[i]/g.nodes[i].M;
 SGS s(g,u,nullptr,0.);s.maximum=0.;
 for(size_t i=0;i<N;i++){if(!std::isfinite(nut[i])||nut[i]<0.)return 2;s.nut[i]=nut[i];s.maximum=std::max(s.maximum,nut[i]);}
 s.acceleration(force);*bound=0.;for(size_t i=0;i<g.D;i++)*bound=std::max(*bound,s.row_bound(g.nodes[i],s.maximum));
 std::copy(force.begin(),force.end(),out);return 0;
}
extern "C" __declspec(dllexport) int pm_flow_rates_les(const double*xf,const double*yf,const double*zf,double rho,double nu,
 double cv,int ground,const double*momentum,double*stats){
 using namespace penalty_les;
 int code=pm_flow_rates(xf,yf,zf,rho,nu,momentum,stats);if(code)return code;
 Grid g;if(!build(g,xf,yf,zf,rho))return 2;std::vector<double>u(g.D);
 for(size_t i=0;i<g.D;i++)u[i]=momentum[i]/g.nodes[i].M;
 // All-fluid evaluation bounds any viscosity subsequently zeroed in solid cells.
 SGS s(g,u,nullptr,cv);if(!s.valid)return 3;double rb=0.,gb=0.;
 #pragma omp parallel for reduction(max:rb,gb) schedule(static)
 for(size_t i=0;i<g.D;i++){
  const Node&q=g.nodes[i];rb=std::max(rb,s.row_bound(q,s.maximum));
  if(ground&&q.a!=1&&q.t[1]==0)gb=std::max(gb,nu*g.dualArea(q,1)/(.5*g.h[1][0]*q.V));
 }
 stats[1]+=rb+gb;stats[2]=stats[0]+stats[1];stats[7]=s.maximum;stats[8]=rb;stats[9]=gb;stats[10]=cv;return 0;
}
extern "C" __declspec(dllexport) int pm_prepare_les(const double*xf,const double*yf,const double*zf,double rho,double dt,
 double nu,double wind,double cv,int ground,const double*momentum,const unsigned char*fluid,
 double*prepared,double*nut,double*stats,double*sgsstats){
 using namespace penalty_les;
 Grid g;if(!build(g,xf,yf,zf,rho)||!momentum||!fluid||!prepared||!nut||!stats||!sgsstats)return 1;
 for(size_t i=0;i<N;i++)if(fluid[i]>1)return 2;
 std::vector<double>zero(g.D,0.),u(g.D),force;
 int code=pm_prepare(xf,yf,zf,rho,dt,0.,nu,wind,momentum,zero.data(),zero.data(),prepared,stats);if(code)return 10+code;
 for(size_t i=0;i<g.D;i++)u[i]=momentum[i]/g.nodes[i].M;
 SGS s(g,u,fluid,cv);if(!s.valid)return 20;s.acceleration(force);
 double rb=0.,gb=0.,work=0.,mean=0.,volume=0.;
 for(size_t i=0;i<g.D;i++){
  const Node&q=g.nodes[i];rb=std::max(rb,s.row_bound(q,s.maximum));work+=rho*q.V*u[i]*force[i];
  if(ground&&q.a!=1&&q.t[1]==0)gb=std::max(gb,nu*g.dualArea(q,1)/(.5*g.h[1][0]*q.V));
 }
 double diss=rho*s.dissipation();
 if(!std::isfinite(diss+work)||std::abs(work+diss)>1e-9*std::max(1.,diss))return 21;
 stats[1]+=dt*(rb+gb);if(stats[0]+stats[1]>.5)return 22;
 #pragma omp parallel for schedule(static)
 for(size_t i=0;i<g.D;i++){
  const Node&q=g.nodes[i];double extra=force[i];
  if(ground&&q.a!=1&&q.t[1]==0){double ub=q.a==0?wind:0.;extra+=nu*g.dualArea(q,1)*(ub-u[i])/(.5*g.h[1][0]*q.V);}
  prepared[g.D+i]+=extra;
  double raw=momentum[i]+dt*q.M*(prepared[i]+prepared[g.D+i]);
  double constrained=g.closed(q)?0.:g.inlet(q)?q.M*wind:raw;
  prepared[2*g.D+i]=raw;prepared[3*g.D+i]=constrained;prepared[4*g.D+i]=constrained;prepared[5*g.D+i]=constrained;
  prepared[6*g.D+i]=constrained-raw;prepared[7*g.D+i]=0.;
 }
 #pragma omp parallel for collapse(3) schedule(static)
 for(int i=0;i<NX;i++)for(int j=0;j<NY;j++)for(int k=0;k<NZ;k++){
  int t[3]={i,j,k};double divergence=0.;
  for(int a=0;a<3;a++){int hi[3]={i,j,k};hi[a]++;size_t lo=g.off[a]+ix(g.ds[a],t),up=g.off[a]+ix(g.ds[a],hi);
   double area=g.area(a,t);divergence+=area*(prepared[5*g.D+up]/g.nodes[up].M)-area*(prepared[5*g.D+lo]/g.nodes[lo].M);}
  prepared[10*g.D+ix(n,t)]=-divergence/dt;
 }
 for(size_t cell=0;cell<N;cell++){
  int t[3]={int(cell/(NY*NZ)),int((cell/NZ)%NY),int(cell%NZ)};double v=g.h[0][t[0]]*g.h[1][t[1]]*g.h[2][t[2]];mean+=v*s.nut[cell];volume+=v;
 }
 std::copy(s.nut.begin(),s.nut.end(),nut);
 sgsstats[0]=s.maximum;sgsstats[1]=mean/volume;sgsstats[2]=work;sgsstats[3]=diss;sgsstats[4]=rb;sgsstats[5]=gb;sgsstats[6]=cv;
 sgsstats[7]=std::count(fluid,fluid+N,(unsigned char)1);sgsstats[8]=ground;return 0;
}
