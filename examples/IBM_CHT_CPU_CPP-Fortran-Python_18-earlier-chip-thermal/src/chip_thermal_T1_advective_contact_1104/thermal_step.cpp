#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#define API extern "C" __declspec(dllexport)
static bool overlaps(const void* a,size_t na,const void* b,size_t nb){auto x=(uintptr_t)a,y=(uintptr_t)b;return x<y+nb&&y<x+na;}
API int lica_thermal_step(int n,int m,int nf,const int* edge,const double* G,const double* C,const double* old,const double* source,double H,double inlet,double dt,double* theta,double* Q,double* boundary){
 if(n<1||n>64||m<1||m>128||nf<1||nf>n||!edge||!G||!C||!old||!source||!theta||!Q||!boundary)return 1;
 if(!std::isfinite(H)||std::abs(H)>10000||!std::isfinite(inlet)||std::abs(inlet)>10000||!std::isfinite(dt)||dt<=0||dt>1)return 2;
 const void* ptr[8]={edge,G,C,old,source,theta,Q,boundary};size_t len[8]={size_t(2*m)*sizeof(int),size_t(m)*8,size_t(n)*8,size_t(n)*8,size_t(n)*8,size_t(n)*8,size_t(m)*8,16};
 for(int o=5;o<8;++o)for(int i=0;i<o;++i)if(overlaps(ptr[o],len[o],ptr[i],len[i]))return 3;
 for(int i=0;i<n;++i)if(!std::isfinite(C[i])||C[i]<=0||C[i]>1e6||!std::isfinite(old[i])||std::abs(old[i])>1e4||!std::isfinite(source[i])||std::abs(source[i])>1e6)return 4;
 for(int e=0;e<m;++e)if(edge[2*e]<0||edge[2*e]>=n||edge[2*e+1]<0||edge[2*e+1]>=n||edge[2*e]==edge[2*e+1]||!std::isfinite(G[e])||G[e]<=0||G[e]>1e6)return 5;
 std::vector<double> A(n*n,0),b(n),x(n),q(m);
 for(int i=0;i<n;++i){A[i*n+i]=C[i]/dt;b[i]=(C[i]/dt)*old[i]+source[i];}
 for(int e=0;e<m;++e){int i=edge[2*e],j=edge[2*e+1];A[i*n+i]+=G[e];A[j*n+j]+=G[e];A[i*n+j]-=G[e];A[j*n+i]-=G[e];}
 double h=std::abs(H);
 for(int i=0;i<nf;++i){A[i*n+i]+=h;int upstream=H>=0?i-1:i+1;if(upstream>=0&&upstream<nf)A[i*n+upstream]-=h;else b[i]+=h*inlet;}
 for(int k=0;k<n;++k){int p=k;for(int i=k+1;i<n;++i)if(std::abs(A[i*n+k])>std::abs(A[p*n+k]))p=i;if(!std::isfinite(A[p*n+k])||A[p*n+k]<=0)return 6;if(p!=k){for(int j=k;j<n;++j)std::swap(A[k*n+j],A[p*n+j]);std::swap(b[k],b[p]);}for(int i=k+1;i<n;++i){double f=A[i*n+k]/A[k*n+k];for(int j=k+1;j<n;++j)A[i*n+j]-=f*A[k*n+j];b[i]-=f*b[k];}}
 for(int i=n-1;i>=0;--i){double z=b[i];for(int j=i+1;j<n;++j)z-=A[i*n+j]*x[j];x[i]=z/A[i*n+i];if(!std::isfinite(x[i]))return 7;}
 for(int e=0;e<m;++e)q[e]=G[e]*(x[edge[2*e]]-x[edge[2*e+1]]);
 double out=h*x[H>=0?nf-1:0],in=h*inlet;if(!std::isfinite(out)||!std::isfinite(in))return 7;
 std::copy(x.begin(),x.end(),theta);std::copy(q.begin(),q.end(),Q);boundary[0]=in;boundary[1]=out;return 0;
}
