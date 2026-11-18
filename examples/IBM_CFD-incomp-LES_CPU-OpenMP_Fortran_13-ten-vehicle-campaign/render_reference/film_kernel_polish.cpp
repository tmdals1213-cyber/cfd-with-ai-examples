// NEW smooth scalar movie projection + variable-camera surface/Q renderer.
// Stored CFD pressure/velocity is read-only. Q surfaces are diagnostic geometry.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif
using V=std::array<double,3>;using F=std::array<int32_t,3>;using C=std::array<unsigned char,3>;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
template<class T>void rd(std::istream&f,T*p,size_t n){f.read(reinterpret_cast<char*>(p),sizeof(T)*n);need(bool(f),"Short input");}
template<class T>void wr(std::ostream&f,const T*p,size_t n){f.write(reinterpret_cast<const char*>(p),sizeof(T)*n);need(bool(f),"Output write");}
void fresh(std::ofstream&f,const std::string&name){std::ifstream old(name,std::ios::binary);need(!old.good(),"Output exists");old.close();f.open(name,std::ios::binary);need(bool(f),"Output open");}
double dot(const V&a,const V&b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
V cross(const V&a,const V&b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double symlog(double v){double a=std::abs(v),adj=.6/.9;return a<=5?v*adj:std::copysign(5*(adj+std::log10(a/5)),v);}
C color(const std::vector<C>&lut,double scalar,bool vort){
 double u=vort?(std::log10(std::max(scalar,10.))-1)/(std::log10(300.)-1):(symlog(scalar)+symlog(120))/(2*symlog(120));
 u=std::clamp(u,0.,1.)*4095;int i=std::min(4094,(int)u);double a=u-i;C c;
 for(int k=0;k<3;k++)c[k]=(unsigned char)std::lround((1-a)*lut[i][k]+a*lut[i+1][k]);return c;
}
void project(const char*cache,const char*batch){
 std::ifstream f(cache,std::ios::binary);int64_t magic;int32_t nn,np;rd(f,&magic,1);rd(f,&nn,1);rd(f,&np,1);
 need(magic==2026110210LL&&nn>0&&nn<=2500000&&np>0&&np<=1000000,"Sparse map header");
 std::vector<int32_t>indices(size_t(nn)*8);std::vector<double>w(size_t(nn)*8),nx(nn);
 rd(f,indices.data(),indices.size());rd(f,w.data(),w.size());rd(f,nx.data(),nx.size());need(f.peek()==EOF,"Map trailing bytes");f.close();
 for(int i=0;i<nn;i++){double sum=0;need(std::isfinite(nx[i])&&std::abs(nx[i])<=1.000000000001,"Display normal");
  for(int j=0;j<8;j++){size_t k=size_t(i)*8+j;need(indices[k]>=0&&indices[k]<np&&std::isfinite(w[k])&&w[k]>=0,"Sparse operand");sum+=w[k];}need(std::abs(sum-1)<1e-12,"Weight sum");}
 std::ifstream list(batch);int n;list>>n;need(bool(list)&&n>0&&n<=16,"Projection batch cap");std::string input,output;
 for(int j=0;j<n;j++){
  list>>input>>output;need(bool(list),"Projection filenames");std::ifstream frame(input,std::ios::binary);int32_t count;double t,p0;
  rd(frame,&magic,1);rd(frame,&count,1);rd(frame,&t,1);rd(frame,&p0,1);need(magic==2026110201LL&&count==np&&std::isfinite(t)&&std::isfinite(p0),"Pressure frame");
  std::vector<double>p(np);rd(frame,p.data(),p.size());need(frame.peek()==EOF,"Pressure trailing bytes");frame.close();
  for(double x:p)need(std::isfinite(x),"Nonfinite pressure");std::vector<float>values(nn);
  for(int i=0;i<nn;i++){double value=0;for(int k=0;k<8;k++)value+=w[size_t(i)*8+k]*p[indices[size_t(i)*8+k]];values[i]=(float)(-(value-p0)*925.9259259259259*nx[i]);need(std::isfinite(values[i]),"Nonfinite visual scalar");}
  std::ofstream out;fresh(out,output);magic=2026110209LL;wr(out,&magic,1);wr(out,&nn,1);wr(out,&t,1);wr(out,values.data(),values.size());out.close();need(bool(out),"Projection close");
 }
 std::cout<<"{\"projected_frames\":"<<n<<",\"CFD_calls\":0}\n";
}
struct Mesh{std::vector<V>v;std::vector<F>f;std::vector<float>s;};
Mesh body(const char*path){
 std::ifstream in(path,std::ios::binary);int64_t magic;int32_t nv,nf;rd(in,&magic,1);rd(in,&nv,1);rd(in,&nf,1);
 need(magic==2026110208LL&&nv>0&&nv<=2500000&&nf>0&&nf<=1500000,"Static body bounds");
 Mesh m;m.v.resize(nv);m.f.resize(nf);m.s.resize(nv);rd(in,m.v.data(),nv);rd(in,m.s.data(),nv);rd(in,m.f.data(),nf);need(in.peek()==EOF,"Body trailing bytes");return m;
}
Mesh wake(const char*path){
 Mesh m;if(std::string(path)=="-")return m;std::ifstream list(path);int n;list>>n;need(bool(list)&&n>=0&&n<=24,"Wake piece list");std::string name;
 for(int i=0;i<n;i++){
  list>>name;need(bool(list),"Wake piece path");std::ifstream f(name,std::ios::binary);int64_t magic;int32_t nv,nf;rd(f,&magic,1);rd(f,&nv,1);rd(f,&nf,1);
  need(magic==2026110207LL&&nv>0&&nv<=750000&&nf>0&&nf<=250000,"Indexed wake bounds");
  size_t offset=m.v.size();std::vector<std::array<float,4>>v(nv);std::vector<F>faces(nf);rd(f,v.data(),nv);rd(f,faces.data(),nf);need(f.peek()==EOF,"Wake trailing bytes");
  for(auto a:v){for(float x:a)need(std::isfinite(x),"Wake finite");m.v.push_back({a[0],a[1],a[2]});m.s.push_back(a[3]);}
  for(auto a:faces){for(int k=0;k<3;k++){need(a[k]>=0&&a[k]<nv,"Wake index");a[k]+=(int32_t)offset;}m.f.push_back(a);}
 }
 need(m.f.size()<=5000000&&m.v.size()<=5000000,"Whole wake bound");return m;
}
std::vector<C> lut(const char*path){std::ifstream f(path,std::ios::binary);std::vector<C>v(4096);rd(f,v.data(),v.size());need(f.peek()==EOF,"LUT bound");return v;}
void render(const Mesh&m,const Mesh&q,const std::vector<C>&jet,const std::vector<C>&magma,const std::string&dynamic,const std::string&output,int W,int H,double ppm,double az,double el,const V&center,int flow){
 need(W>0&&W<=2560&&H>0&&H<=1600&&W*H<=4096000&&ppm>0&&ppm<=500&&flow>=0&&flow<=1,"Camera bounds");
 std::vector<float>scalar=m.s;if(dynamic!="-"){std::ifstream f(dynamic,std::ios::binary);int64_t magic;int32_t n;double t;rd(f,&magic,1);rd(f,&n,1);rd(f,&t,1);need(magic==2026110209LL&&n==(int)m.v.size()&&std::isfinite(t),"Dynamic scalar header");rd(f,scalar.data(),n);need(f.peek()==EOF,"Dynamic scalar bytes");}
 double pi=std::acos(-1.);az*=pi/180;el*=pi/180;V eye={std::cos(el)*std::cos(az),std::sin(el),std::cos(el)*std::sin(az)},right=cross({0,1,0},eye);double rn=std::sqrt(dot(right,right));for(double&x:right)x/=rn;V up=cross(eye,right);
 auto screen=[&](const std::vector<V>&vertices){std::vector<V>s(vertices.size());for(size_t i=0;i<vertices.size();i++){V a;for(int k=0;k<3;k++){need(std::isfinite(vertices[i][k]),"Body finite");a[k]=vertices[i][k]-center[k];}s[i]={W*.5+ppm*dot(a,right),H*.55-ppm*dot(a,up),dot(a,eye)};}return s;};
 std::vector<double>depth(W*H,-std::numeric_limits<double>::infinity());std::vector<unsigned char>pixels(size_t(W)*H*3,0);
 std::vector<std::array<float,4>>accum(flow?size_t(W)*H:0,{0,0,0,0});
 auto draw=[&](const Mesh&mesh,const std::vector<float>&values,bool transparent){
  std::vector<V>s=screen(mesh.v);
  for(const auto&face:mesh.f){
   for(int k=0;k<3;k++)need(face[k]>=0&&face[k]<(int)s.size(),"Body index");
   const V&a=s[face[0]],&b=s[face[1]],&c=s[face[2]];
   double den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1]);if(std::abs(den)<1e-14)continue;
   int x0=std::max(0,(int)std::floor(std::min({a[0],b[0],c[0]}))),x1=std::min(W-1,(int)std::ceil(std::max({a[0],b[0],c[0]})));
   int y0=std::max(0,(int)std::floor(std::min({a[1],b[1],c[1]}))),y1=std::min(H-1,(int)std::ceil(std::max({a[1],b[1],c[1]})));
   for(int y=y0;y<=y1;y++)for(int x=x0;x<=x1;x++){
    double px=x+.5,py=y+.5,wa=((b[1]-c[1])*(px-c[0])+(c[0]-b[0])*(py-c[1]))/den,wb=((c[1]-a[1])*(px-c[0])+(a[0]-c[0])*(py-c[1]))/den,wc=1-wa-wb;
    if(wa<0||wb<0||wc<0)continue;double z=wa*a[2]+wb*b[2]+wc*c[2];size_t id=size_t(y)*W+x;if(z<=depth[id])continue;
    double value=wa*values[face[0]]+wb*values[face[1]]+wc*values[face[2]];need(std::isfinite(value),"Scalar finite");C rgb=color(transparent?magma:jet,value,transparent);
    if(!transparent){depth[id]=z;for(int k=0;k<3;k++)pixels[id*3+k]=rgb[k];}
    else{for(int k=0;k<3;k++)accum[id][k]+=rgb[k];accum[id][3]+=1;}
   }
  }
 };
 draw(m,scalar,false);if(flow){draw(q,q.s,true);
 double sigma=std::max(.75,ppm*.018);int radius=std::min(16,(int)std::ceil(3*sigma));std::vector<float>kernel(2*radius+1);double norm=0;
 for(int k=-radius;k<=radius;k++){kernel[k+radius]=(float)std::exp(-.5*k*k/(sigma*sigma));norm+=kernel[k+radius];}for(float&x:kernel)x=(float)(x/norm);
 std::vector<std::array<float,4>>tmp(accum.size(),{0,0,0,0});
 for(int y=0;y<H;y++)for(int x=0;x<W;x++)for(int k=-radius;k<=radius;k++){int xx=std::clamp(x+k,0,W-1);float w=kernel[k+radius];for(int c=0;c<4;c++)tmp[size_t(y)*W+x][c]+=w*accum[size_t(y)*W+xx][c];}
 std::fill(accum.begin(),accum.end(),std::array<float,4>{0,0,0,0});
 for(int y=0;y<H;y++)for(int x=0;x<W;x++)for(int k=-radius;k<=radius;k++){int yy=std::clamp(y+k,0,H-1);float w=kernel[k+radius];for(int c=0;c<4;c++)accum[size_t(y)*W+x][c]+=w*tmp[size_t(yy)*W+x][c];}
 for(size_t i=0;i<accum.size();i++)if(accum[i][3]>1e-6){double count=accum[i][3],a=std::min(.65,1-std::exp(-.055*count));for(int k=0;k<3;k++)pixels[i*3+k]=(unsigned char)std::lround((1-a)*pixels[i*3+k]+a*accum[i][k]/count);}}
 if(output=="-"){std::cout<<"P6\n"<<W<<" "<<H<<"\n255\n";wr(std::cout,pixels.data(),pixels.size());std::cout.flush();}
 else{std::ofstream out;fresh(out,output);out<<"P6\n"<<W<<" "<<H<<"\n255\n";wr(out,pixels.data(),pixels.size());out.close();need(bool(out),"PPM close");}
}
int main(int argc,char**argv){try{
 need(argc>1,"Mode");std::string mode=argv[1];
 if(mode=="project"){need(argc==4,"Project arguments");project(argv[2],argv[3]);}
 else if(mode=="render"){
#ifdef _WIN32
  _setmode(_fileno(stdout),_O_BINARY);
#endif
  need(argc==7,"Render arguments");Mesh m=body(argv[2]),q=wake(argv[3]);auto jet=lut(argv[4]),magma=lut(argv[5]);std::ifstream list(argv[6]);int n;list>>n;need(bool(list)&&n>0&&n<=12,"Render batch cap");
  for(int i=0;i<n;i++){std::string scalar,out;int W,H,flow;double ppm,az,el;V c;list>>scalar>>out>>W>>H>>ppm>>az>>el>>c[0]>>c[1]>>c[2]>>flow;need(bool(list),"Render row");render(m,q,jet,magma,scalar,out,W,H,ppm,az,el,c,flow);}
  std::cout<<"{\"rendered_frames\":"<<n<<",\"variable_camera\":true,\"source_geometry_modified\":false}\n";
 }else need(false,"Unknown fixed mode");return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
