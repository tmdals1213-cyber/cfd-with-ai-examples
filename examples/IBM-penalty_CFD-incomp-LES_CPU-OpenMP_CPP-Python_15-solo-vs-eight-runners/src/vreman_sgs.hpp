// Vreman (2004), anisotropic cell widths. Independent implementation of the
// mathematical model; C_v=.035 matches the previously selected LICA setting.
// SGS stress is the negative adjoint of the staggered symmetric gradient:
// this retains shear cross terms and has a nonpositive discrete energy work.
namespace penalty_les {
using namespace pacemaker_air_v1;
inline double viscosity(const double A[3][3],const double width[3],double coefficient){
 double norm=0.,minor=0.;
 for(int a=0;a<3;a++)for(int b=0;b<3;b++)norm+=A[a][b]*A[a][b];
 if(norm==0.)return 0.;
 // Cauchy-Binet sum of squared 2x2 minors of the width-weighted gradient.
 // It equals B_beta without subtracting nearly equal positive numbers.
 for(int a=0;a<3;a++)for(int b=a+1;b<3;b++)
  for(int j=0;j<3;j++)for(int k=j+1;k<3;k++){
   double v=width[j]*width[k]*(A[a][j]*A[b][k]-A[a][k]*A[b][j]);minor+=v*v;
  }
 return coefficient*std::sqrt(minor/norm);
}
struct SGS {
 const Grid&g;const std::vector<double>&u;
 std::vector<std::array<double,3>>center;
 std::vector<double>nut;double maximum=0.;bool valid=true;
 SGS(const Grid&grid,const std::vector<double>&velocity,const unsigned char*fluid,double coefficient):g(grid),u(velocity),center(N),nut(N,0.){
  if(!std::isfinite(coefficient)||coefficient<0.){valid=false;return;}
  #pragma omp parallel for schedule(static)
  for(size_t c=0;c<N;c++){
   parallel_evidence::note(0);
   int t[3]={int(c/(NY*NZ)),int((c/NZ)%NY),int(c%NZ)};
   for(int a=0;a<3;a++){int hi[3]={t[0],t[1],t[2]};hi[a]++;center[c][a]=.5*(u[g.off[a]+ix(g.ds[a],t)]+u[g.off[a]+ix(g.ds[a],hi)]);}
  }
  int bad=0;double maxNut=0.;
  #pragma omp parallel for schedule(static) reduction(max:bad,maxNut)
  for(size_t c=0;c<N;c++){
   if(fluid&&!fluid[c])continue;
   int t[3]={int(c/(NY*NZ)),int((c/NZ)%NY),int(c%NZ)};
   double A[3][3],h[3]={g.h[0][t[0]],g.h[1][t[1]],g.h[2][t[2]]};
   for(int a=0;a<3;a++)for(int b=0;b<3;b++){
    int lo[3]={t[0],t[1],t[2]},hi[3]={t[0],t[1],t[2]};
    if(a==b){hi[b]++;A[a][b]=(u[g.off[a]+ix(g.ds[a],hi)]-u[g.off[a]+ix(g.ds[a],t)])/h[b];}
    else if(t[b]==0){hi[b]++;A[a][b]=(center[ix(n,hi)][a]-center[c][a])/(.5*(h[b]+g.h[b][1]));}
    else if(t[b]==n[b]-1){lo[b]--;A[a][b]=(center[c][a]-center[ix(n,lo)][a])/(.5*(h[b]+g.h[b][t[b]-1]));}
    else{
     lo[b]--;hi[b]++;double dl=.5*(h[b]+g.h[b][t[b]-1]),dr=.5*(h[b]+g.h[b][t[b]+1]);
     A[a][b]=-dr/(dl*(dl+dr))*center[ix(n,lo)][a]+(dr-dl)/(dl*dr)*center[c][a]+dl/(dr*(dl+dr))*center[ix(n,hi)][a];
    }
   }
   double value=viscosity(A,h,coefficient);
   if(!std::isfinite(value)||value<0.)bad=1;else{nut[c]=value;maxNut=std::max(maxNut,value);}
  }
  maximum=maxNut;valid=!bad;
 }
 double normal(int a,const int*t)const{
  int hi[3]={t[0],t[1],t[2]};hi[a]++;
  return 2*nut[ix(n,t)]*(u[g.off[a]+ix(g.ds[a],hi)]-u[g.off[a]+ix(g.ds[a],t)])/g.h[a][t[a]];
 }
 double edge_nut(int a,int b,const int*t)const{
  double value=0.;
  for(int da=0;da<2;da++)for(int db=0;db<2;db++){
   int q[3]={t[0],t[1],t[2]};q[a]-=1-da;q[b]-=1-db;
   double wa=g.h[a][t[a]-da]/(g.h[a][t[a]-1]+g.h[a][t[a]]);
   double wb=g.h[b][t[b]-db]/(g.h[b][t[b]-1]+g.h[b][t[b]]);
   value+=wa*wb*nut[ix(n,q)];
  }
  return value;
 }
 double shear(int a,int b,const int*t)const{
  if(t[a]<=0||t[a]>=n[a]||t[b]<=0||t[b]>=n[b])return 0.;
  int ua[3]={t[0],t[1],t[2]},ub[3]={t[0],t[1],t[2]};ua[b]--;ub[a]--;
  double gradient=(u[g.off[a]+ix(g.ds[a],t)]-u[g.off[a]+ix(g.ds[a],ua)])/g.span(b,t[b])
                 +(u[g.off[b]+ix(g.ds[b],t)]-u[g.off[b]+ix(g.ds[b],ub)])/g.span(a,t[a]);
  return edge_nut(a,b,t)*gradient;
 }
 double row_bound(const Node&q,double maxNu)const{
  int a=q.a,j=q.t[a];double value=0.;double ar=g.area(a,q.t);
  if(j>0)value+=4*maxNu*ar/(g.h[a][j-1]*q.V);
  if(j<n[a])value+=4*maxNu*ar/(g.h[a][j]*q.V);
  if(j>0&&j<n[a])for(int b=0;b<3;b++)if(b!=a){
   for(int side=0;side<2;side++){
    int edge=q.t[b]+side;if(edge<=0||edge>=n[b])continue;
    value+=maxNu/g.h[b][q.t[b]]*(2/g.span(b,edge)+2/g.span(a,j));
   }
  }
  return value;
 }
 void acceleration(std::vector<double>&out)const{
  out.resize(g.D);
  #pragma omp parallel for schedule(static)
  for(size_t z=0;z<g.D;z++){
   parallel_evidence::note(0);
   const Node&q=g.nodes[z];int a=q.a,j=q.t[a];double value=0.;
   if(j<n[a])value+=normal(a,q.t)/g.span(a,j);
   if(j>0){int lo[3]={q.t[0],q.t[1],q.t[2]};lo[a]--;value-=normal(a,lo)/g.span(a,j);}
   for(int b=0;b<3;b++)if(b!=a){int hi[3]={q.t[0],q.t[1],q.t[2]};hi[b]++;value+=(shear(a,b,hi)-shear(a,b,q.t))/g.h[b][q.t[b]];}
   out[z]=value;
  }
 }
 double dissipation()const{
  double sum=0.;
  for(size_t cell=0;cell<N;cell++){
   int t[3]={int(cell/(NY*NZ)),int((cell/NZ)%NY),int(cell%NZ)};
   double vol=g.h[0][t[0]]*g.h[1][t[1]]*g.h[2][t[2]];
   for(int a=0;a<3;a++){int hi[3]={t[0],t[1],t[2]};hi[a]++;double grad=(u[g.off[a]+ix(g.ds[a],hi)]-u[g.off[a]+ix(g.ds[a],t)])/g.h[a][t[a]];sum+=2*nut[cell]*vol*grad*grad;}
  }
  for(int a=0;a<3;a++)for(int b=a+1;b<3;b++){
   int c=3-a-b;
   for(int i=1;i<n[a];i++)for(int j=1;j<n[b];j++)for(int k=0;k<n[c];k++){
    int t[3];t[a]=i;t[b]=j;t[c]=k;double visc=edge_nut(a,b,t);if(visc==0)continue;
    double stress=shear(a,b,t);sum+=stress*stress/visc*g.span(a,i)*g.span(b,j)*g.h[c][k];
   }
  }
  return sum;
 }
};
}
