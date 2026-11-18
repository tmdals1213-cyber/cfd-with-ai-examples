#include "nearest_geometry.hpp"
#include <filesystem>
#include <map>
#include <tuple>
#include <cstddef>
#include <climits>
#include "rnu_exact.hpp"
namespace RunnerNonzeroUnionV1 {
// NEW explicit out-of-band diagnostic context. Original predicates/math remain unchanged.
struct RdDiagRecord{int info[16];double data[20];};
static_assert(sizeof(RdDiagRecord)==224,"diagnostic ABI");
struct RdDiagnosticContext{RdDiagRecord record;const std::vector<int>*orientation;V a,b;double time;uint64_t epoch;int owner;};
static thread_local RdDiagnosticContext*rd_active_diagnostic=nullptr;
static void rd_diag_capture(int site,int component,int count,int orient,int before,V query,V direction,int triangle,int sign,int eventIndex,double eventT)noexcept{
 RdDiagnosticContext*c=rd_active_diagnostic;if(!c||c->record.info[1])return;
 for(int k=0;k<16;k++)c->record.info[k]=0;for(int k=0;k<20;k++)c->record.data[k]=0;
 int*z=c->record.info;z[0]=1;z[1]=1;z[2]=site;z[3]=component;z[4]=count;z[5]=orient;z[6]=before;z[7]=triangle;z[8]=sign;z[9]=eventIndex;z[10]=-1;z[11]=c->owner;z[12]=(int)c->epoch;z[13]=-1;
 bool isA=true,isB=true;for(int k=0;k<3;k++){isA=isA&&query.x[k]==c->a.x[k];isB=isB&&query.x[k]==c->b.x[k];}z[10]=isA?0:isB?1:-1;
 double*x=c->record.data;x[0]=c->time;for(int k=0;k<3;k++){x[1+k]=query.x[k];x[4+k]=c->a.x[k];x[7+k]=c->b.x[k];x[10+k]=direction.x[k];}x[13]=eventT;
}
static void rd_diag_query8(const std::vector<int64_t>&counts,V p,V d)noexcept{
 if(!rd_active_diagnostic)return;for(size_t c=0;c<counts.size();c++)if(std::abs(counts[c])>1){int orientation=0;const auto*v=rd_active_diagnostic->orientation;if(v&&c<v->size())orientation=(*v)[c];rd_diag_capture(1,(int)c,counts[c],orientation,counts[c],p,d,-1,0,-1,-1);return;}
}

// Component-union mask only. Nearest RAW component is not an exposed union wall.
struct UnionMesh{
 Mesh m;int nc;std::vector<int>component,component_actor;
 bool raybox(const Box&b,V p,V d,double eps)const{
  double lo=0,hi=std::numeric_limits<double>::infinity();
  for(int a=0;a<3;a++){double u=(b.lo.x[a]-eps-p.x[a])/d.x[a],v=(b.hi.x[a]+eps-p.x[a])/d.x[a];if(u>v)std::swap(u,v);lo=std::max(lo,u);hi=std::min(hi,v);if(hi<lo)return false;}return true;
 }
 int triangle_ray(int t,V p,V d,double eps,int&crossing)const{
  V a=m.pos(m.faces[t][0]),e1=m.pos(m.faces[t][1])-a,e2=m.pos(m.faces[t][2])-a,n=cross(e1,e2),h=cross(d,e2),s=p-a;
  double det=dot(e1,h),threshold=128*std::numeric_limits<double>::epsilon()*std::sqrt(dot(e1,e1)*dot(e2,e2));
  if(std::abs(det)<=threshold){if(std::abs(dot(s,n))<=eps*std::sqrt(dot(n,n))&&raybox(m.trianglebox(t),p,d,eps))return 2;return 0;}
  double inv=1/det,u=dot(s,h)*inv;V q=cross(s,e1);double v=dot(d,q)*inv,time=dot(e2,q)*inv,tol=128*std::numeric_limits<double>::epsilon();
  if(u< -tol||v< -tol||u+v>1+tol||time< -eps)return 0;
  if(time<=eps||u<=tol||v<=tol||1-u-v<=tol)return 2;
  crossing=dot(n,d)>0?1:-1;return 1;
 }
 bool rayvisit(int node,V p,V d,int owner,double eps,std::vector<int64_t>&counts)const{
  const Node&b=m.nodes[node];if(!raybox(b.box,p,d,eps))return true;
  if(b.left<0){for(int i=b.lo;i<b.hi;i++){int t=m.order[i];if(owner&&m.actor[t]!=owner)continue;int s=0,code=triangle_ray(t,p,d,eps,s);if(code==2)return false;if(code==1)counts[component[t]]+=s;}return true;}
  return rayvisit(b.left,p,d,owner,eps,counts)&&rayvisit(b.right,p,d,owner,eps,counts);
 }
 int query(V p,double time,int owner,int&inside,int&valid,int&components_inside,Hit&hit)const{
  if(m.current<0||time!=m.times[m.current])return 3;if(!finite(p)||owner<0)return 1;Hit h;m.visit(0,p,owner,h);if(h.tri<0)return 4;
  double scale=1;for(int a=0;a<3;a++)scale=std::max(scale,std::max(std::abs(p.x[a]),m.nodes[0].box.hi.x[a]-m.nodes[0].box.lo.x[a]));double eps=128*std::numeric_limits<double>::epsilon()*scale;
  if(h.d2<=eps*eps)return 6; // Any raw component surface: conservative refusal.
  V directions[5]={V(.873,.321,.365),V(-.431,.812,.394),V(.267,-.517,.813),V(-.697,-.379,.609),V(.581,.743,-.333)};
  std::vector<int64_t>reference;int accepted=0;
  for(V d:directions){d=d*(1/std::sqrt(dot(d,d)));std::vector<int64_t>counts(nc,0);if(!rayvisit(0,p,d,owner,eps,counts))continue;
   // NEW explicit nonzero-winding region: preserve every signed integer count.
   if(accepted&&counts!=reference)return 7;reference=counts;accepted++;
  }
  if(accepted<2)return 7;int count=0;for(int64_t n:reference)count+=(n!=0);
  inside=count>0;valid=accepted;components_inside=count;hit=h;return 0;
 }
};

// NEW exposed-union segment crossing provider. No wholegrid mask/IBM/CFD release.
struct EdgeEvent{int info[4];double data[13];};static_assert(sizeof(EdgeEvent)==120,"edge event ABI");

static int rd_diag_event_index(const std::vector<EdgeEvent>&all,const EdgeEvent&e)noexcept{for(size_t k=0;k<all.size();k++)if(all[k].info[0]==e.info[0]&&all[k].info[1]==e.info[1]&&all[k].data[0]==e.data[0])return (int)k;return -1;}
struct SegmentMesh:UnionMesh{
 std::vector<int>orientation;
 void prepare(int snap,double time){m.prepare(snap,time);std::vector<long double>vol(nc,0);std::vector<V>anchor(nc);std::vector<int>set(nc,0);for(int t=0;t<m.nf;t++){int c=component[t];V a=m.pos(m.faces[t][0]),b=m.pos(m.faces[t][1]),d=m.pos(m.faces[t][2]);if(!set[c]){anchor[c]=a;set[c]=1;}a=a-anchor[c];b=b-anchor[c];d=d-anchor[c];vol[c]+=(long double)a.x[0]*((long double)b.x[1]*d.x[2]-(long double)b.x[2]*d.x[1])+(long double)a.x[1]*((long double)b.x[2]*d.x[0]-(long double)b.x[0]*d.x[2])+(long double)a.x[2]*((long double)b.x[0]*d.x[1]-(long double)b.x[1]*d.x[0]);}orientation.resize(nc);for(int c=0;c<nc;c++){if(!std::isfinite(vol[c]))throw std::runtime_error("nonfinite signed component volume");orientation[c]=vol[c]>0?1:vol[c]<0?-1:0;}}
 bool segmentbox(const Box&b,V p,V d,double L,double eps)const{double lo=0,hi=L;for(int a=0;a<3;a++){if(d.x[a]==0){if(p.x[a]<b.lo.x[a]-eps||p.x[a]>b.hi.x[a]+eps)return false;continue;}double v=(b.lo.x[a]-eps-p.x[a])/d.x[a],w=(b.hi.x[a]+eps-p.x[a])/d.x[a];if(v>w)std::swap(v,w);lo=std::max(lo,v);hi=std::min(hi,w);if(hi<lo)return false;}return true;}
 int counts(V p,double time,int owner,double eps,std::vector<int64_t>&out)const{int in,v,n;Hit h;int code=query(p,time,owner,in,v,n,h);if(code)return code;V ds[5]={V(.873,.321,.365),V(-.431,.812,.394),V(.267,-.517,.813),V(-.697,-.379,.609),V(.581,.743,-.333)};for(V d:ds){d=d*(1/std::sqrt(dot(d,d)));std::vector<int64_t>c(nc,0);if(rayvisit(0,p,d,owner,eps,c)){/* NEW nonzero fill: endpoint integer layers preserved, never clamped. */out=c;return 0;}}return 7;}
 int visit_segment(int node,V p,V d,double L,int owner,double eps,std::vector<EdgeEvent>&events)const{
  const Node&z=m.nodes[node];if(!segmentbox(z.box,p,d,L,eps))return 0;if(z.left>=0){int c=visit_segment(z.left,p,d,L,owner,eps,events);return c?c:visit_segment(z.right,p,d,L,owner,eps,events);}
  const double tol=128*std::numeric_limits<double>::epsilon();
  for(int k=z.lo;k<z.hi;k++){int t=m.order[k];if(owner&&m.actor[t]!=owner)continue;V a=m.pos(m.faces[t][0]),e1=m.pos(m.faces[t][1])-a,e2=m.pos(m.faces[t][2])-a,n=cross(e1,e2),h=cross(d,e2),s=p-a;double det=dot(e1,h),limit=tol*std::sqrt(dot(e1,e1)*dot(e2,e2));if(std::abs(det)<=limit){if(std::abs(dot(s,n))<=eps*std::sqrt(dot(n,n))&&segmentbox(m.trianglebox(t),p,d,L,eps))return 6;continue;}double inv=1/det,u=dot(s,h)*inv;V q=cross(s,e1);double v=dot(d,q)*inv,r=dot(e2,q)*inv;if(u< -tol||v< -tol||u+v>1+tol||r< -eps||r>L+eps)continue;if(r<=eps||r>=L-eps||u<=tol||v<=tol||1-u-v<=tol)return 6;if(events.size()>=256)return 10;
   EdgeEvent ev;int c=component[t],sgn=dot(n,d)>0?1:-1;ev.info[0]=t;ev.info[1]=c;ev.info[2]=sgn;ev.info[3]=0;ev.data[0]=r/L;double bary[3]={1-u-v,u,v};V point=p+d*r,wall=m.vel(m.faces[t][0])*bary[0]+m.vel(m.faces[t][1])*bary[1]+m.vel(m.faces[t][2])*bary[2],normal=n*((orientation[c]?orientation[c]:1)/std::sqrt(dot(n,n)));for(int j=0;j<3;j++){ev.data[1+j]=point.x[j];ev.data[4+j]=wall.x[j];ev.data[7+j]=normal.x[j];ev.data[10+j]=bary[j];}events.push_back(ev);
  }return 0;
 }
 int segment(V a,V b,double time,int owner,int&raw,int&start,int&end,std::vector<EdgeEvent>&exposed,std::vector<EdgeEvent>*rawEvents=nullptr,std::vector<int64_t>*startCounts=nullptr,std::vector<int64_t>*endCounts=nullptr,std::vector<int>*eventGroups=nullptr)const{
  if(m.current<0||time!=m.times[m.current])return 3;if(!finite(a)||!finite(b)||owner<0)return 1;V d=b-a;double L=std::sqrt(dot(d,d));if(!(L>0)||!std::isfinite(L))return 1;d=d*(1/L);double scale=std::max(1.,L);for(int j=0;j<3;j++)scale=std::max(scale,std::max(std::abs(a.x[j]),std::abs(b.x[j])));double eps=128*std::numeric_limits<double>::epsilon()*scale;std::vector<int64_t>c0,c1;int code=counts(a,time,owner,eps,c0);if(code)return code;code=counts(b,time,owner,eps,c1);if(code)return code;std::vector<EdgeEvent>all;code=visit_segment(0,a,d,L,owner,eps,all);if(code)return code;std::sort(all.begin(),all.end(),[](const EdgeEvent&a,const EdgeEvent&b){return a.data[0]<b.data[0];});auto same=[&](const EdgeEvent&u,const EdgeEvent&v){auto x=RNUExact::parameter(m.pos(m.faces[u.info[0]][0]),m.pos(m.faces[u.info[0]][1]),m.pos(m.faces[u.info[0]][2]),a,b);auto y=RNUExact::parameter(m.pos(m.faces[v.info[0]][0]),m.pos(m.faces[v.info[0]][1]),m.pos(m.faces[v.info[0]][2]),a,b);return RNUExact::equal(x,y);};for(size_t k=1;k<all.size();k++)if((all[k].data[0]-all[k-1].data[0])*L<=eps&&!same(all[k],all[k-1]))return 9;
  auto occupied=[](const std::vector<int64_t>&c){for(int64_t v:c)if(v!=0)return true;return false;};bool initial=occupied(c0);
  const std::vector<int64_t>initialCounts=c0;std::vector<int>groupID(all.size(),-1);int ng=0;
  for(size_t k=0;k<all.size();){size_t j=k+1;while(j<all.size()&&same(all[k],all[j]))j++;const std::vector<int64_t>before=c0;bool was=occupied(c0);for(size_t z=k;z<j;z++){c0[all[z].info[1]]-=all[z].info[2];groupID[z]=ng;}bool now=occupied(c0);if(was!=now){std::vector<size_t>active;for(size_t z=k;z<j;z++){int c=all[z].info[1];if((before[c]!=0)!=(c0[c]!=0))active.push_back(z);}if(active.empty())return 7;EdgeEvent chosen=all[active[0]];chosen.info[3]=was?1:-1;V chosenN(chosen.data[7],chosen.data[8],chosen.data[9]),chosenV(chosen.data[4],chosen.data[5],chosen.data[6]);if((dot(chosenN,d)>0?1:-1)!=chosen.info[3])chosenN=chosenN*(-1.);for(size_t aidx=1;aidx<active.size();aidx++){size_t z=active[aidx];V n(all[z].data[7],all[z].data[8],all[z].data[9]),v(all[z].data[4],all[z].data[5],all[z].data[6]);if((dot(n,d)>0?1:-1)!=chosen.info[3])n=n*(-1.);double tol=64*std::numeric_limits<double>::epsilon();if(dot(n-chosenN,n-chosenN)>tol*tol||dot(v-chosenV,v-chosenV)>tol*tol*std::max(1.,dot(chosenV,chosenV)))return 11;}for(int z=0;z<3;z++)chosen.data[7+z]=chosenN.x[z];exposed.push_back(chosen);}k=j;ng++;}
  if(c0!=c1)return 7;raw=(int)all.size();start=initial;end=occupied(c1);if(rawEvents)*rawEvents=all;if(startCounts)*startCounts=initialCounts;if(endCounts)*endCounts=c1;if(eventGroups)*eventGroups=groupID;return 0;
 }
};

#include <map>
#include <tuple>
struct DynamicRunner{SegmentMesh body;std::vector<int>canonical;std::string context;uint64_t epoch=0;int builds=0,refits=0;double time=0;};
extern "C" __declspec(dllexport) int rnu1_rd_load(const char*path,const char*context,void**out){try{
 if(!path||!context||!out)return 1;std::ifstream f(std::filesystem::u8path(path),std::ios::binary);char magic[8],ctx[64];int n[3];rd(f,magic,8);rd(f,n,3);rd(f,ctx,64);
 bool fixture=std::string(magic,8)=="RNUW0001"&&n[0]>=4&&n[0]<=128&&n[1]>=4&&n[1]<=256&&n[2]>=1&&n[2]<=16,actual=std::string(magic,8)=="RDYN0001"&&n[0]==28374&&n[1]==49112&&n[2]==95;
 if((!fixture&&!actual)||std::string(ctx,64)!=std::string(context))return 2;
 DynamicRunner d;d.context=std::string(context);SegmentMesh&u=d.body;Mesh&m=u.m;m.nv=n[0];m.nf=n[1];m.ns=1;u.nc=n[2];m.faces.resize(m.nf);u.component.resize(m.nf);m.actor.resize(m.nf);m.part.resize(m.nf);d.canonical.resize(m.nv);u.component_actor.assign(u.nc,-1);m.times.resize(1);m.positions.resize(m.nv);m.velocity.resize(m.nv);
 rd(f,m.faces.data(),m.nf);rd(f,u.component.data(),m.nf);rd(f,m.actor.data(),m.nf);rd(f,m.part.data(),m.nf);rd(f,d.canonical.data(),m.nv);char tail;if(f.read(&tail,1))return 2;
 for(int c:d.canonical)if(c<0||c>=m.nv)return 2;
 std::map<std::tuple<int,int,int>,std::pair<int,int>>edges;
 for(int t=0;t<m.nf;t++){for(int v:m.faces[t])if(v<0||v>=m.nv)return 2;int c=u.component[t];if(c<0||c>=u.nc||m.actor[t]<=0||m.part[t]<0)return 2;if(u.component_actor[c]>=0&&u.component_actor[c]!=m.actor[t])return 2;u.component_actor[c]=m.actor[t];
  for(int j=0;j<3;j++){int a=d.canonical[m.faces[t][j]],b=d.canonical[m.faces[t][(j+1)%3]];if(a==b)return 2;auto&v=edges[std::make_tuple(c,std::min(a,b),std::max(a,b))];v.first++;v.second+=a<b?1:-1;}}
 for(int a:u.component_actor)if(a<0)return 2;for(auto&e:edges)if(e.second.first!=2||e.second.second!=0)return 2;
 *out=new DynamicRunner(std::move(d));return 0;
 }catch(...){return 1;}}
extern "C" __declspec(dllexport) int rnu1_rd_update(void*handle,int nv,const double*xyz,const double*velocity,double time,uint64_t epoch,const char*context){try{
 if(!handle||!xyz||!velocity||!context||!std::isfinite(time)||time<0)return 1;DynamicRunner&d=*(DynamicRunner*)handle;
 if(nv!=d.body.m.nv||std::string(context)!=d.context)return 2;if(epoch!=d.epoch+1||(d.epoch&&time<=d.time))return 3;
 SegmentMesh candidate=d.body;candidate.m.times[0]=time;
 for(int i=0;i<nv;i++){V p(xyz[3*i],xyz[3*i+1],xyz[3*i+2]),v(velocity[3*i],velocity[3*i+1],velocity[3*i+2]);if(!finite(p)||!finite(v))return 2;candidate.m.positions[i]=p;candidate.m.velocity[i]=v;}
 bool initial=candidate.m.nodes.empty();candidate.prepare(0,time);d.body=std::move(candidate);d.epoch=epoch;d.time=time;if(initial)d.builds++;else d.refits++;return 0;
 }catch(...){return 2;}}
extern "C" __declspec(dllexport) int rnu1_rd_state(const void*handle,uint64_t*epoch,double*time,int*builds,int*refits){if(!handle||!epoch||!time||!builds||!refits)return 1;const DynamicRunner&d=*(const DynamicRunner*)handle;if(!d.epoch)return 3;*epoch=d.epoch;*time=d.time;*builds=d.builds;*refits=d.refits;return 0;}
extern "C" __declspec(dllexport) int rnu1_rd_segment(const void*handle,const double*a,const double*b,double time,uint64_t epoch,int owner,int capacity,int*raw,int*number,int*start,int*end,EdgeEvent*out){try{
 if(!handle||!a||!b||!raw||!number||!start||!end||!out||capacity<0)return 1;const DynamicRunner&d=*(const DynamicRunner*)handle;if(!d.epoch||epoch!=d.epoch||time!=d.time)return 3;
 int r=0,s=0,e=0;std::vector<EdgeEvent>x;int code=d.body.segment(V(a[0],a[1],a[2]),V(b[0],b[1],b[2]),time,owner,r,s,e,x);if(code)return code;if((int)x.size()>capacity)return 10;*raw=r;*number=x.size();*start=s;*end=e;std::copy(x.begin(),x.end(),out);return 0;
 }catch(...){return 1;}}
extern "C" __declspec(dllexport) void rnu1_rd_release(void*h){delete (DynamicRunner*)h;}
#include <cstddef>
extern "C" __declspec(dllexport) int rnu1_rd_abi(int*i,int*r,int*e,int*o,int*u){if(!i||!r||!e||!o||!u)return 1;*i=sizeof(int);*r=sizeof(double);*e=sizeof(EdgeEvent);*o=offsetof(EdgeEvent,data);*u=sizeof(uint64_t);return 0;}

// NEW ABI: binds one caller-owned per-thread context, invokes the unchanged atomic rnu1_rd_segment.
extern "C" __declspec(dllexport) int rnu1_rd_segment_diagnostic(const void*handle,const double*a,const double*b,double time,uint64_t epoch,int owner,int capacity,int*raw,int*number,int*start,int*end,EdgeEvent*out,RdDiagRecord*diagnostic){
 if(!diagnostic)return 1;RdDiagnosticContext context{};context.orientation=handle?&((const DynamicRunner*)handle)->body.orientation:nullptr;context.a=a?V(a[0],a[1],a[2]):V();context.b=b?V(b[0],b[1],b[2]):V();context.time=time;context.epoch=epoch;context.owner=owner;
 context.record.info[0]=1;context.record.info[3]=-1;context.record.info[7]=-1;context.record.info[9]=-1;context.record.info[10]=-1;context.record.info[13]=-1;
 struct Bind{RdDiagnosticContext*old;Bind(RdDiagnosticContext*p):old(rd_active_diagnostic){rd_active_diagnostic=p;}~Bind(){rd_active_diagnostic=old;}}bind(&context);
 int code=rnu1_rd_segment(handle,a,b,time,epoch,owner,capacity,raw,number,start,end,out);*diagnostic=context.record;return code;
}


extern "C" __declspec(dllexport) int rnu1_rd_method(int*version,int*countBytes,int*fillRule){if(!version||!countBytes||!fillRule)return 1;*version=1;*countBytes=sizeof(int64_t);*fillRule=1;return 0;}
extern "C" __declspec(dllexport) int rnu1_rd_winding(const void*handle,const double*p,double time,uint64_t epoch,int owner,int countLength,int*inside,int*accepted,int*components_inside,int64_t*counts){try{if(!handle||!p||!inside||!accepted||!components_inside||!counts)return 1;const DynamicRunner&d=*(const DynamicRunner*)handle;if(!d.epoch||epoch!=d.epoch||time!=d.time)return 3;if(countLength!=d.body.nc)return 2;int in=0,valid=0,n=0;Hit h;int code=d.body.query(V(p[0],p[1],p[2]),time,owner,in,valid,n,h);if(code)return code;std::vector<int64_t>w;code=d.body.counts(V(p[0],p[1],p[2]),time,owner,0.,w);if(code)return code;*inside=in;*accepted=valid;*components_inside=n;std::copy(w.begin(),w.end(),counts);return 0;}catch(...){return 1;}}
extern "C" __declspec(dllexport) int rnu1_rd_segment_full(const void*handle,const double*a,const double*b,double time,uint64_t epoch,int owner,int capacity,int countLength,int*raw,int*number,int*start,int*end,EdgeEvent*exposed,EdgeEvent*rawEvents,int64_t*startCounts,int64_t*endCounts,int*groupIndex){try{if(!handle||!a||!b||!raw||!number||!start||!end||!exposed||!rawEvents||!startCounts||!endCounts||!groupIndex||capacity<0)return 1;const DynamicRunner&d=*(const DynamicRunner*)handle;if(!d.epoch||epoch!=d.epoch||time!=d.time)return 3;if(countLength!=d.body.nc)return 2;int r=0,s=0,e=0;std::vector<EdgeEvent>x,all;std::vector<int64_t>ca,cb;std::vector<int>groups;int code=d.body.segment(V(a[0],a[1],a[2]),V(b[0],b[1],b[2]),time,owner,r,s,e,x,&all,&ca,&cb,&groups);if(code)return code;if(r>capacity||int(x.size())>capacity)return 10;*raw=r;*number=x.size();*start=s;*end=e;std::copy(x.begin(),x.end(),exposed);std::copy(all.begin(),all.end(),rawEvents);std::copy(ca.begin(),ca.end(),startCounts);std::copy(cb.begin(),cb.end(),endCounts);std::copy(groups.begin(),groups.end(),groupIndex);return 0;}catch(const std::runtime_error&e){std::string s=e.what();return s.find("PREDICATE_")==0?12:1;}catch(...){return 1;}}
extern "C" __declspec(dllexport) int rnu1_counts_transition(int nc,int ne,const int64_t*before,const int*component,const int*sign,int64_t*after,int*maskBefore,int*maskAfter){try{if(!before||!component||!sign||!after||!maskBefore||!maskAfter)return 1;if(nc<1||nc>16||ne<0||ne>256)return 2;std::vector<int64_t>x(before,before+nc);int mb=0,ma=0;for(int64_t v:x)mb|=v!=0;for(int i=0;i<ne;i++){int c=component[i],s=sign[i];if(c<0||c>=nc||(s!=1&&s!=-1))return 3;if((s==1&&x[c]==INT64_MIN)||(s==-1&&x[c]==INT64_MAX))return 4;x[c]-=s;}for(int64_t v:x)ma|=v!=0;std::copy(x.begin(),x.end(),after);*maskBefore=mb;*maskAfter=ma;return 0;}catch(...){return 1;}}
} // namespace RunnerNonzeroUnionV1
