#include <vector>
#include <cmath>
#include <cstdint>
#include <limits>
#include <algorithm>
struct T3Term {int j; long double a;};
static bool ok(double v){return std::isfinite(v);}
static bool span_overlap(const void* a,std::size_t an,const void* b,std::size_t bn){if(!an||!bn)return false;std::uintptr_t x=(std::uintptr_t)a,y=(std::uintptr_t)b;return x<y+bn&&y<x+an;}
extern "C" int thermal_t3_abi(int* out,int len){if(!out||len!=4)return 1;out[0]=sizeof(int);out[1]=sizeof(double);out[2]=3;out[3]=23;return 0;}
extern "C" int thermal_t3_step(const int* m,int ml,const double* ctl,int cl,const int* kind,const double* cell,const int* edge,const double* geom,const int* adv,const double* mdot,const int* boundary,const double* bc,double* tout,double* qcond,double* qadv,double* qbnd,double* stats){
 if(!m||!ctl||ml!=23||cl!=4)return 1;
 const int n=m[1],ne=m[2],na=m[3],nb=m[4];
 if(m[0]!=3||n<1||n>4096||ne<0||ne>16384||na<0||na>16384||nb<0||nb>8192||m[5]<1||m[5]!=m[6]||m[5]!=m[7]||m[5]!=m[8]||m[9]<1||m[9]>4096)return 2;
 int expected[13]={n,6*n,2*ne,4*ne,2*na,na,nb,4*nb,n,ne,na,nb,8};for(int k=0;k<13;++k)if(m[10+k]!=expected[k])return 1;
 if(!kind||!cell||!tout||!stats||(ne&&(!edge||!geom||!qcond))||(na&&(!adv||!mdot||!qadv))||(nb&&(!boundary||!bc||!qbnd)))return 1;
 if(!ok(ctl[0])||ctl[0]<=0||!ok(ctl[1])||ctl[1]<0||!ok(ctl[2])||ctl[2]<=0||!ok(ctl[3])||ctl[3]<=0)return 3;
 const void* ins[12]={m,ctl,kind,cell,edge,geom,adv,mdot,boundary,bc,0,0};std::size_t inlen[12]={23*sizeof(int),4*sizeof(double),(std::size_t)n*sizeof(int),(std::size_t)6*n*sizeof(double),(std::size_t)2*ne*sizeof(int),(std::size_t)4*ne*sizeof(double),(std::size_t)2*na*sizeof(int),(std::size_t)na*sizeof(double),(std::size_t)nb*sizeof(int),(std::size_t)4*nb*sizeof(double),0,0};
 void* outs[5]={tout,qcond,qadv,qbnd,stats};std::size_t olen[5]={(std::size_t)n*sizeof(double),(std::size_t)ne*sizeof(double),(std::size_t)na*sizeof(double),(std::size_t)nb*sizeof(double),8*sizeof(double)};
 for(int a=0;a<5;++a){for(int b=0;b<10;++b)if(span_overlap(outs[a],olen[a],ins[b],inlen[b]))return 90;for(int b=0;b<a;++b)if(span_overlap(outs[a],olen[a],outs[b],olen[b]))return 90;}
 std::vector<long double> C(n),diag(n),rhs(n),x(n),mass(n,0);std::vector<std::vector<T3Term> > row(n);std::vector<double> G(ne),H(na),HB(nb);double fcp=-1,frho=-1;
 const long double dt=ctl[0];
 for(int i=0;i<n;++i){if(kind[i]!=0&&kind[i]!=1)return 2;for(int k=0;k<6;++k)if(!ok(cell[6*i+k]))return 3;double rho=cell[6*i],cp=cell[6*i+1],v=cell[6*i+2],k=cell[6*i+3],old=cell[6*i+4];if(rho<=0||cp<=0||v<=0||k<=0||old<0)return 2;if(kind[i]){if(fcp<0){fcp=cp;frho=rho;}if(cp!=fcp||rho!=frho)return 2;}C[i]=(long double)rho*cp*v;if(!(C[i]>0)||!std::isfinite(C[i]))return 3;diag[i]=C[i]/dt;rhs[i]=diag[i]*old+cell[6*i+5];x[i]=old;}
 for(int e=0;e<ne;++e){int i=edge[2*e],j=edge[2*e+1];if(i<0||j<0||i>=n||j>=n||i==j)return 2;for(int k=0;k<4;++k)if(!ok(geom[4*e+k]))return 3;double A=geom[4*e],di=geom[4*e+1],dj=geom[4*e+2],r=geom[4*e+3];if(A<=0||di<=0||dj<=0||r<0)return 5;long double den=(long double)di/cell[6*i+3]+r+(long double)dj/cell[6*j+3];G[e]=(double)((long double)A/den);if(!ok(G[e])||G[e]<=0)return 5;diag[i]+=G[e];diag[j]+=G[e];row[i].push_back(T3Term{j,-(long double)G[e]});row[j].push_back(T3Term{i,-(long double)G[e]});}
 for(int a=0;a<na;++a){int i=adv[2*a],j=adv[2*a+1];double f=mdot[a];if(i<0||j<0||i>=n||j>=n||i==j||!kind[i]||!kind[j])return 2;if(!ok(f))return 3;mass[i]+=f;mass[j]-=f;int donor=f>=0?i:j,recv=f>=0?j:i;H[a]=std::fabs(f)*cell[6*donor+1];if(!ok(H[a]))return 3;diag[donor]+=H[a];row[recv].push_back(T3Term{donor,-(long double)H[a]});}
 for(int b=0;b<nb;++b){int i=boundary[b];if(i<0||i>=n)return 2;for(int k=0;k<4;++k)if(!ok(bc[4*b+k]))return 3;double f=bc[4*b],tin=bc[4*b+1],gb=bc[4*b+2],tb=bc[4*b+3];if(tin<0||tb<0||gb<0||(!kind[i]&&f!=0))return 2;mass[i]+=f;HB[b]=std::fabs(f)*cell[6*i+1];diag[i]+=gb;rhs[i]+=(long double)gb*tb;if(f>=0)diag[i]+=HB[b];else rhs[i]+=(long double)HB[b]*tin;}
 for(int i=0;i<n;++i)if(std::fabs(mass[i])>ctl[1])return 4;
 long double bn=0;for(int i=0;i<n;++i){if(!(diag[i]>0)||!std::isfinite(diag[i])||!std::isfinite(rhs[i]))return 3;bn=std::max(bn,std::fabs(rhs[i]));}const long double tol=std::max((long double)ctl[2],(long double)ctl[3]*bn);
 std::vector<double> T(n);long double raw=0;int iterations=0;bool converged=false;
 for(int it=1;it<=m[9];++it){for(int i=0;i<n;++i){long double v=rhs[i];for(std::size_t k=0;k<row[i].size();++k)v-=row[i][k].a*x[row[i][k].j];x[i]=v/diag[i];if(!std::isfinite(x[i]))return 3;}for(int i=0;i<n;++i)T[i]=(double)x[i];raw=0;for(int i=0;i<n;++i){long double r=diag[i]*T[i]-rhs[i];for(std::size_t k=0;k<row[i].size();++k)r+=row[i][k].a*T[row[i][k].j];raw=std::max(raw,std::fabs(r));}iterations=it;if(raw<=tol){converged=true;break;}}
 if(!converged)return 6;for(int i=0;i<n;++i)if(!ok(T[i])||T[i]<0)return 7;
 std::vector<double> qc(ne),qa(na),qb(nb);std::vector<long double> net(n,0);long double de=0,src=0,bo=0,local=0;double tmin=T[0],tmax=T[0];
 for(int e=0;e<ne;++e){int i=edge[2*e],j=edge[2*e+1];qc[e]=G[e]*(T[i]-T[j]);net[i]+=qc[e];net[j]-=qc[e];}
 for(int a=0;a<na;++a){int i=adv[2*a],j=adv[2*a+1];qa[a]=(mdot[a]>=0?1.:-1.)*H[a]*T[mdot[a]>=0?i:j];net[i]+=qa[a];net[j]-=qa[a];}
 for(int b=0;b<nb;++b){int i=boundary[b];qb[b]=(bc[4*b]>=0?HB[b]*T[i]:-HB[b]*bc[4*b+1])+bc[4*b+2]*(T[i]-bc[4*b+3]);net[i]+=qb[b];bo+=qb[b];}
 for(int i=0;i<n;++i){long double d=C[i]*(T[i]-cell[6*i+4]);de+=d;src+=cell[6*i+5];local=std::max(local,std::fabs(d-dt*(cell[6*i+5]-net[i])));tmin=std::min(tmin,T[i]);tmax=std::max(tmax,T[i]);}
 double st[8]={(double)iterations,(double)raw,(double)local,(double)de,(double)src,(double)bo,tmin,tmax};for(int k=0;k<8;++k)if(!ok(st[k]))return 3;
 for(int i=0;i<n;++i)tout[i]=T[i];for(int e=0;e<ne;++e)qcond[e]=qc[e];for(int a=0;a<na;++a)qadv[a]=qa[a];for(int b=0;b<nb;++b)qbnd[b]=qb[b];for(int k=0;k<8;++k)stats[k]=st[k];return 0;
}
