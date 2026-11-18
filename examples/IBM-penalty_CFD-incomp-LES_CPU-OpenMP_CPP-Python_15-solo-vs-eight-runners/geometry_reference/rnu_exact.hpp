#pragma once
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
// NEW expansion predicates. A failed safe arithmetic-range certificate refuses;
// it never perturbs a point, rounds an event time, or substitutes another ray.
namespace RNUExact {
typedef std::vector<double> E;
static void sum2(double a,double b,double&s,double&e){s=a+b;if(!std::isfinite(s))throw std::runtime_error("PREDICATE_RANGE");double bv=s-a,av=s-bv,br=b-bv,ar=a-av;e=ar+br;}
static double product_checked(double a,double b){if(a!=0.&&b!=0.){int ea=0,eb=0;std::frexp(a,&ea);std::frexp(b,&eb);if(ea+eb<-1000||ea+eb>1000)throw std::runtime_error("PREDICATE_RANGE");}double r=a*b;if(!std::isfinite(r)||(r==0.&&a!=0.&&b!=0.))throw std::runtime_error("PREDICATE_RANGE");return r;}
static void prod2(double a,double b,double&p,double&e){p=product_checked(a,b);const double split=134217729.;double ac=product_checked(split,a),bc=product_checked(split,b),ahi=ac-(ac-a),alo=a-ahi,bhi=bc-(bc-b),blo=b-bhi;double err1=p-product_checked(ahi,bhi),err2=err1-product_checked(alo,bhi),err3=err2-product_checked(ahi,blo);e=product_checked(alo,blo)-err3;if(!std::isfinite(e))throw std::runtime_error("PREDICATE_RANGE");}
static E scalar(double x){if(!std::isfinite(x))throw std::runtime_error("PREDICATE_RANGE");return x==0.?E():E(1,x);}
static E grow(const E&a,double b){E o;double q=b;for(double x:a){double s,e;sum2(q,x,s,e);if(e)o.push_back(e);q=s;}if(q)o.push_back(q);return o;}
static E add(E a,const E&b){for(double x:b)a=grow(a,x);return a;}
static E neg(E a){for(double&x:a)x=-x;return a;}
static E sub(const E&a,const E&b){return add(a,neg(b));}
static E mul(const E&a,const E&b){E o;for(double x:a)for(double y:b){double p,e;prod2(x,y,p,e);if(e)o=grow(o,e);if(p)o=grow(o,p);}return o;}
static int sign(const E&a){for(size_t i=a.size();i>0;i--)if(a[i-1]!=0.)return a[i-1]>0.?1:-1;return 0;}
struct QV{E x[3];};
static QV point(const V&a){QV r;for(int i=0;i<3;i++)r.x[i]=scalar(a.x[i]);return r;}
static QV minus(const QV&a,const QV&b){QV r;for(int i=0;i<3;i++)r.x[i]=sub(a.x[i],b.x[i]);return r;}
static QV crossq(const QV&a,const QV&b){QV r;r.x[0]=sub(mul(a.x[1],b.x[2]),mul(a.x[2],b.x[1]));r.x[1]=sub(mul(a.x[2],b.x[0]),mul(a.x[0],b.x[2]));r.x[2]=sub(mul(a.x[0],b.x[1]),mul(a.x[1],b.x[0]));return r;}
static E dotq(const QV&a,const QV&b){return add(add(mul(a.x[0],b.x[0]),mul(a.x[1],b.x[1])),mul(a.x[2],b.x[2]));}
struct Ratio{E n,d;};
static Ratio parameter(V a,V b,V c,V p,V q){QV aa=point(a),e1=minus(point(b),aa),e2=minus(point(c),aa),dd=minus(point(q),point(p)),ss=minus(point(p),aa);Ratio r;r.d=dotq(e1,crossq(dd,e2));r.n=dotq(e2,crossq(ss,e1));if(sign(r.d)==0)throw std::runtime_error("PREDICATE_SINGULAR");return r;}
static bool equal(const Ratio&a,const Ratio&b){return sign(sub(mul(a.n,b.d),mul(b.n,a.d)))==0;}
}
