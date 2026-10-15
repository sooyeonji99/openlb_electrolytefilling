#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <queue>
#include <random>
#include <stdexcept>
#include <vector>
namespace cell {
constexpr double width=4.3e-3,sourceWidth=1e-3,sourceHeight=.5e-3;
constexpr double porousHeight=4e-3/3, height=sourceHeight+porousHeight;
constexpr unsigned seed=20261014;
struct Circle {double x,y,r;};
struct Geometry {
  bool openOuter;
  double dx,pitch,targetPorosity,latticeHeight; int columns,rows;
  std::vector<Circle> circles;
  std::vector<std::vector<int>> bins;
  static constexpr double minDiameter=50e-6,maxDiameter=200e-6,minGap=10e-6;
  Geometry(double d,double porosity=.3,double =50e-6,bool open=true):openOuter(open),dx(d),pitch(maxDiameter+minGap),targetPorosity(porosity),latticeHeight(std::round(height/d)*d),columns(std::ceil(width/pitch)),rows(std::ceil(height/pitch)) {
    if(!std::isfinite(d)||d<=0||d>5e-6)throw std::runtime_error("dx must be in (0,5] um for mixed particles");
    std::mt19937 gen(seed);
    auto random=[&](){return double(gen())/double(gen.max());};
    const double wall=10e-6; // fixed physical packing margin across raster resolutions
    const double areaMargin=5e-6; // physical area reference; packing wall clearance stays 10 um
    const double area=(width-2*areaMargin)*(height-2*areaMargin)-sourceWidth*(sourceHeight-areaMargin);
    const double target=(1-porosity)*area;
    double occupied=0;
    while(occupied<target) {
      double r=(minDiameter+(maxDiameter-minDiameter)*random())/2;
      if(occupied+std::acos(-1.)*r*r>target) {
        r=std::sqrt((target-occupied)/std::acos(-1.));
        if(r<minDiameter/2)break;
      }
      circles.push_back({wall+r+(width-2*wall-2*r)*random(),wall+r+(height-2*wall-2*r)*random(),r});
      occupied+=std::acos(-1.)*r*r;
    }
    // Relax independent random centres into non-overlapping polydisperse disks.
    // Clearance is a physical distance, independent of requested raster resolution.
    double worst=1.;
    for(int iteration=0;iteration<30000 && worst>1e-10;++iteration) {
      rebuild();worst=0;
      for(size_t i=0;i<circles.size();++i) {
        auto& a=circles[i];int bx=std::clamp(int(a.x/pitch),0,columns-1),by=std::clamp(int(a.y/pitch),0,rows-1);
        for(int y=std::max(0,by-1);y<=std::min(rows-1,by+1);++y)
          for(int x=std::max(0,bx-1);x<=std::min(columns-1,bx+1);++x)
            for(int j:bins[y*columns+x])if(j>int(i)) {
              auto& b=circles[j];double vx=b.x-a.x,vy=b.y-a.y,length=std::hypot(vx,vy);
              double overlap=a.r+b.r+minGap-length;
              if(overlap>0) {
                worst=std::max(worst,overlap);
                if(length<1e-15){vx=1;vy=0;length=1;}
                double move=.45*overlap/length;a.x-=move*vx;a.y-=move*vy;b.x+=move*vx;b.y+=move*vy;
              }
            }
      }
      for(auto& c:circles) {
        c.x=std::clamp(c.x,wall+c.r,width-wall-c.r);c.y=std::clamp(c.y,wall+c.r,height-wall-c.r);
        double left=(width-sourceWidth)/2,right=(width+sourceWidth)/2;
        double nearX=std::clamp(c.x,left,right),nearY=std::clamp(c.y,0.,sourceHeight);
        double vx=c.x-nearX,vy=c.y-nearY,length=std::hypot(vx,vy),clearance=c.r+minGap;
        if(length<clearance) {
          worst=std::max(worst,clearance-length);
          if(length>1e-15){c.x+=vx*(clearance-length)/length;c.y+=vy*(clearance-length)/length;}
          else {
            double l=c.x-left,r=right-c.x,t=sourceHeight-c.y;
            if(l<=r&&l<=t)c.x=left-clearance;
            else if(r<=t)c.x=right+clearance;
            else c.y=sourceHeight+clearance;
          }
        }
      }
    }
    if(worst>1e-10)throw std::runtime_error("packing clearance did not converge");
    rebuild();
  }
  void rebuild() {
    bins.assign(columns*rows,{});
    for(size_t i=0;i<circles.size();++i){auto c=circles[i];int x=std::clamp(int(c.x/pitch),0,columns-1),y=std::clamp(int(c.y/pitch),0,rows-1);bins[y*columns+x].push_back(i);}
  }
  bool inSource(double x,double y)const {return std::abs(x-width/2)<=sourceWidth/2+dx*1e-7 && y>=-dx*1e-7 && y<=sourceHeight+dx*1e-7;}
  int material(double x,double y)const {
    x=std::round(x/dx)*dx;y=std::round(y/dx)*dx;
    if(x<0||x>width||y<0||y>latticeHeight+dx*1e-7)return 0;
    if(openOuter) {
      int ix=std::lround(x/dx),iy=std::lround(y/dx),nx=std::lround(width/dx),ny=std::lround(latticeHeight/dx);
      if(ix==0||ix==nx||iy==0||iy==ny)return 0; // exterior ghost ring
      bool side=ix==1||ix==nx-1, horizontal=iy==1||iy==ny-1;
      if(side&&horizontal)return 2; // four discrete corner closure nodes
      if(iy==1)return 2; // entire bottom is a wall; internal reservoir is unchanged
      if(ix==1)return 4;
      if(ix==nx-1)return 5;
      if(iy==ny-1)return 7;
    }
    // Closed outer mode retains the original two-node perimeter.
    if(x<1.51*dx||x>width-1.51*dx||y<1.51*dx||y>latticeHeight-1.51*dx)return 2;
    if(inSource(x,y))return 3;
    int i=int(x/pitch),j=int(y/pitch);
    for(int b=std::max(0,j-1);b<=std::min(rows-1,j+1);++b)
      for(int a=std::max(0,i-1);a<=std::min(columns-1,i+1);++a) {
        for(int index:bins[b*columns+a]) {
        auto c=circles[index];
        if(std::hypot(x-c.x,y-c.y)<=c.r)return 2;
        }
      }
    return 1;
  }
  void write(const std::string& path)const {
    std::filesystem::create_directories(path);
    int nx=std::lround(width/dx)+1,ny=std::lround(latticeHeight/dx)+1;
    std::vector<uint8_t> mask(nx*ny),seen(nx*ny);std::queue<int> q;
    long pore=0,solid=0,reservoir=0,connected=0;bool top=false;
    for(int j=0;j<ny;++j)for(int i=0;i<nx;++i){int k=j*nx+i;mask[k]=material(i*dx,j*dx);if(mask[k]==1)++pore;if(mask[k]==3){++reservoir;seen[k]=1;q.push(k);}if(mask[k]==2&&!inSource(i*dx,j*dx)&&i>1&&i<nx-2&&j>1&&j<ny-2)++solid;}
    while(!q.empty()){int k=q.front();q.pop();int i=k%nx,j=k/nx;if(mask[k]==1){++connected;if(j==ny-3)top=true;}for(int n:{i>0?k-1:-1,i+1<nx?k+1:-1,j>0?k-nx:-1,j+1<ny?k+nx:-1})if(n>=0&&!seen[n]&&(mask[n]==1||mask[n]==3)){seen[n]=1;q.push(n);}}
    double porosity=double(pore)/(pore+solid),fraction=double(connected)/pore;
    bool pass=reservoir>0&&std::abs(porosity-targetPorosity)<.005&&fraction>.99&&top;
    std::ofstream csv(path+"/geometry_validation.csv");csv<<std::setprecision(17)<<"status,nx,ny,dx_m,width_m,height_m,source_cells,pore_cells,porosity,source_connected_pore_fraction,reaches_top_interior\n"<<(pass?"PASS":"FAIL")<<','<<nx<<','<<ny<<','<<dx<<','<<width<<','<<latticeHeight<<','<<reservoir<<','<<pore<<','<<porosity<<','<<fraction<<','<<top<<'\n';
    std::ofstream vtk(path+"/geometry.vtk");vtk<<"# vtk DataFile Version 3.0\nCell horizontal material\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "<<nx<<' '<<ny<<" 1\nORIGIN 0 0 0\nSPACING "<<dx<<' '<<dx<<" 1\nPOINT_DATA "<<nx*ny<<"\nSCALARS material unsigned_char 1\nLOOKUP_TABLE default\n";for(auto m:mask)vtk<<int(m)<<'\n';
    std::ofstream particles(path+"/circles.csv");particles<<std::setprecision(17)<<"x_m,y_m,radius_m,removed_for_source\n";for(auto c:circles){double sx=std::max(std::abs(c.x-width/2)-sourceWidth/2,0.);double sy=std::max(c.y-sourceHeight,0.);particles<<c.x<<','<<c.y<<','<<c.r<<','<<(std::hypot(sx,sy)<=c.r)<<'\n';}
    double smallest=1e9,largest=0,spacing=1e9;
    std::vector<int> histogram(6);
    bool particlesValid=true;
    for(size_t i=0;i<circles.size();++i) {
      auto c=circles[i];double diameter=2*c.r;
      smallest=std::min(smallest,diameter);largest=std::max(largest,diameter);
      ++histogram[std::clamp(int((diameter-minDiameter)/25e-6),0,5)];
      double sx=std::max(std::abs(c.x-width/2)-sourceWidth/2,0.),sy=std::max(c.y-sourceHeight,0.);
      particlesValid &= diameter>=minDiameter-1e-12 && diameter<=maxDiameter+1e-12 && std::hypot(sx,sy)-c.r>=minGap-1e-10;
      for(size_t j=i+1;j<circles.size();++j){auto other=circles[j];spacing=std::min(spacing,std::hypot(c.x-other.x,c.y-other.y)-c.r-other.r);}
    }
    particlesValid &= spacing>=minGap-1e-10;
    std::ofstream stats(path+"/particle_validation.csv");
    stats<<std::setprecision(17)<<"status,particle_count,min_diameter_um,max_diameter_um,min_surface_gap_um,required_gap_um,seed\n"<<(particlesValid?"PASS":"FAIL")<<','<<circles.size()<<','<<smallest*1e6<<','<<largest*1e6<<','<<spacing*1e6<<','<<minGap*1e6<<','<<seed<<'\n';
    std::ofstream hist(path+"/diameter_histogram.csv");hist<<"diameter_min_um,diameter_max_um,count\n";for(int i=0;i<6;++i)hist<<50+25*i<<','<<75+25*i<<','<<histogram[i]<<'\n';
    if(!particlesValid)throw std::runtime_error("particle validation failed");
    if(!pass)throw std::runtime_error("geometry validation failed: see geometry_validation.csv");
  }
};
}
