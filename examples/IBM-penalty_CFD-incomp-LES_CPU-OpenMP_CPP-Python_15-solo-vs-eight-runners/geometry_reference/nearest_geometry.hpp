#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>
// Nearest moving triangle and barycentric wall velocity ONLY. No inside sign,
// union mask, fluid reset, forcing, pressure solve or physical drag computation.
struct V{double x[3];V(double a=0,double b=0,double c=0):x{a,b,c}{}V operator+(V b)const{return V(x[0]+b.x[0],x[1]+b.x[1],x[2]+b.x[2]);}V operator-(V b)const{return V(x[0]-b.x[0],x[1]-b.x[1],x[2]-b.x[2]);}V operator*(double a)const{return V(x[0]*a,x[1]*a,x[2]*a);}};
static_assert(sizeof(V)==24&&sizeof(int)==4&&sizeof(double)==8&&sizeof(std::array<int,3>)==12,"fixed binary ABI");
double dot(V a,V b){return a.x[0]*b.x[0]+a.x[1]*b.x[1]+a.x[2]*b.x[2];}
V cross(V a,V b){return V(a.x[1]*b.x[2]-a.x[2]*b.x[1],a.x[2]*b.x[0]-a.x[0]*b.x[2],a.x[0]*b.x[1]-a.x[1]*b.x[0]);}
bool finite(V a){return std::isfinite(a.x[0])&&std::isfinite(a.x[1])&&std::isfinite(a.x[2]);}
template<class T>void rd(std::ifstream&f,T*p,size_t n){f.read((char*)p,sizeof(T)*n);if(!f)throw std::runtime_error("short sealed input");}
struct Box{V lo,hi;Box():lo(1e300,1e300,1e300),hi(-1e300,-1e300,-1e300){}void add(V p){for(int a=0;a<3;a++){lo.x[a]=std::min(lo.x[a],p.x[a]);hi.x[a]=std::max(hi.x[a],p.x[a]);}}void add(Box b){add(b.lo);add(b.hi);}double lower2(V p)const{double d=0;for(int a=0;a<3;a++){double q=p.x[a]<lo.x[a]?lo.x[a]-p.x[a]:(p.x[a]>hi.x[a]?p.x[a]-hi.x[a]:0);d+=q*q;}return d;}};
struct Node{Box box;int lo,hi,left=-1,right=-1;};
struct Hit{int tri=-1,actor=-1,part=-1;double b[3]={0,0,0};V point,u,normal;double d2=std::numeric_limits<double>::infinity();};
struct Mesh{
 int nv,nf,ns,current=-1;std::string context;std::vector<std::array<int,3>> faces;std::vector<int>actor,part,order;std::vector<double>times;std::vector<V>positions,velocity;std::vector<Node>nodes;
 V pos(int i)const{return positions[(size_t)current*nv+i];}V vel(int i)const{return velocity[(size_t)current*nv+i];}
 Box trianglebox(int t)const{Box b;for(int v:faces[t])b.add(pos(v));return b;}
 V center(int t)const{return (pos(faces[t][0])+pos(faces[t][1])+pos(faces[t][2]))*(1./3);}
 int build(int lo,int hi){
  int id=nodes.size();nodes.push_back(Node());Node q;q.lo=lo;q.hi=hi;for(int i=lo;i<hi;i++)q.box.add(trianglebox(order[i]));
  if(hi-lo>8){int a=0;for(int b=1;b<3;b++)if(q.box.hi.x[b]-q.box.lo.x[b]>q.box.hi.x[a]-q.box.lo.x[a])a=b;int mid=(lo+hi)/2;std::nth_element(order.begin()+lo,order.begin()+mid,order.begin()+hi,[&](int t,int u){double c=center(t).x[a],d=center(u).x[a];return c<d||(c==d&&t<u);});q.left=build(lo,mid);q.right=build(mid,hi);}
  nodes[id]=q;return id;
 }
 Box refit(int n){Node&q=nodes[n];Box b;if(q.left<0){for(int i=q.lo;i<q.hi;i++)b.add(trianglebox(order[i]));}else{b.add(refit(q.left));b.add(refit(q.right));}q.box=b;return b;}
 void prepare(int snap,double time){
  if(snap<0||snap>=ns||!std::isfinite(time)||time!=times[snap])throw std::runtime_error("snapshot-time mismatch");current=snap;
  for(int t=0;t<nf;t++){V n=cross(pos(faces[t][1])-pos(faces[t][0]),pos(faces[t][2])-pos(faces[t][0]));if(!(dot(n,n)>0)||!finite(n))throw std::runtime_error("degenerate triangle refused, not dropped");}
  if(nodes.empty()){order.resize(nf);std::iota(order.begin(),order.end(),0);nodes.reserve(2*nf);build(0,nf);}else refit(0);
 }
 void closest_triangle(int t,V p,Hit&h)const{
  V a=pos(faces[t][0]),b=pos(faces[t][1]),c=pos(faces[t][2]),ab=b-a,ac=c-a,ap=p-a;double d1=dot(ab,ap),d2=dot(ac,ap),w0,w1,w2;
  if(d1<=0&&d2<=0){w0=1;w1=w2=0;}
  else{V bp=p-b;double d3=dot(ab,bp),d4=dot(ac,bp);
   if(d3>=0&&d4<=d3){w1=1;w0=w2=0;}
   else{double vc=d1*d4-d3*d2;
    if(vc<=0&&d1>=0&&d3<=0){w1=d1/(d1-d3);w0=1-w1;w2=0;}
    else{V cp=p-c;double d5=dot(ab,cp),d6=dot(ac,cp);
     if(d6>=0&&d5<=d6){w2=1;w0=w1=0;}
     else{double vb=d5*d2-d1*d6;
      if(vb<=0&&d2>=0&&d6<=0){w2=d2/(d2-d6);w0=1-w2;w1=0;}
      else{double va=d3*d6-d5*d4;
       if(va<=0&&(d4-d3)>=0&&(d5-d6)>=0){w2=(d4-d3)/((d4-d3)+(d5-d6));w1=1-w2;w0=0;}
       else{double den=1/(va+vb+vc);w1=vb*den;w2=vc*den;w0=1-w1-w2;}
      }
     }
    }
   }
  }
  V q=a*w0+b*w1+c*w2;double dd=dot(p-q,p-q);if(dd>h.d2||(dd==h.d2&&h.tri>=0&&t>h.tri))return;
  h.tri=t;h.actor=actor[t];h.part=part[t];h.d2=dd;h.point=q;h.b[0]=w0;h.b[1]=w1;h.b[2]=w2;h.u=vel(faces[t][0])*w0+vel(faces[t][1])*w1+vel(faces[t][2])*w2;V n=cross(ab,ac);h.normal=n*(1/std::sqrt(dot(n,n)));
 }
 void visit(int n,V p,int owner,Hit&h)const{
  const Node&q=nodes[n];if(q.box.lower2(p)>h.d2)return;
  if(q.left<0){for(int i=q.lo;i<q.hi;i++){int t=order[i];if(owner==0||actor[t]==owner)closest_triangle(t,p,h);}return;}
  int a=q.left,b=q.right;if(nodes[b].box.lower2(p)<nodes[a].box.lower2(p))std::swap(a,b);visit(a,p,owner,h);visit(b,p,owner,h);
 }
};
