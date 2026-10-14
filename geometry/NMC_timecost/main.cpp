#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct Options {
  double outerDiameter = 8.6e-3;
  double electrolyteDiameter = 1.785e-3;
  double dx = 2.e-6;
  double targetPorosity = 0.30;
  double particlePitch = 50.e-6;
  double tauLiquid = 0.8;
  double liquidNu = 1.e-6;
  double simulatedTime = 10.;
  // Calibrated from the dx=5 um OpenLB smoke run: 355468 KiB / 741321
  // allocated square-domain nodes = 491 bytes/node. Round upward.
  double bytesPerAllocatedNode = 500.;
  std::string output = "results/geometry_validation";
};

Options parse(int argc, char** argv) {
  Options o;
  for (int i=1; i<argc; ++i) {
    if (i+1 >= argc) throw std::runtime_error("missing value for " + std::string(argv[i]));
    const std::string key = argv[i];
    const std::string value = argv[++i];
    if (key == "--dx-um") o.dx = std::stod(value)*1.e-6;
    else if (key == "--outer-diameter-mm") o.outerDiameter = std::stod(value)*1.e-3;
    else if (key == "--electrolyte-diameter-mm") o.electrolyteDiameter = std::stod(value)*1.e-3;
    else if (key == "--porosity") o.targetPorosity = std::stod(value);
    else if (key == "--particle-pitch-um") o.particlePitch = std::stod(value)*1.e-6;
    else if (key == "--time-s") o.simulatedTime = std::stod(value);
    else if (key == "--output") o.output = value;
    else throw std::runtime_error("unknown option: " + key);
  }
  if (!(o.outerDiameter > 0 && o.electrolyteDiameter > 0 &&
        o.electrolyteDiameter < o.outerDiameter && o.dx > 0 &&
        o.targetPorosity > 0 && o.targetPorosity < 1 && o.particlePitch > 0)) {
    throw std::runtime_error("invalid geometry or numerical option");
  }
  return o;
}

struct Counts {
  std::uint64_t outside{}, electrolyte{}, pore{}, solid{};
};

void writeVtk(const fs::path& path, const std::vector<std::uint8_t>& material,
              int n, double dx) {
  std::ofstream f(path);
  f << "# vtk DataFile Version 3.0\nNMC quarter-sector validation\nASCII\n"
    << "DATASET STRUCTURED_POINTS\nDIMENSIONS " << n << ' ' << n << " 1\n"
    << "ORIGIN 0 0 0\nSPACING " << std::setprecision(17) << dx << ' ' << dx << " 1\n"
    << "POINT_DATA " << static_cast<std::uint64_t>(n)*n << "\n"
    << "SCALARS material unsigned_char 1\nLOOKUP_TABLE default\n";
  for (auto value : material) f << static_cast<int>(value) << '\n';
}

int main(int argc, char** argv) {
  const Options o = parse(argc, argv);
  const double radius = o.outerDiameter/2.;
  const double electrolyteRadius = o.electrolyteDiameter/2.;
  const int n = static_cast<int>(std::ceil(radius/o.dx))+1;
  const double particleRadius = o.particlePitch*std::sqrt((1.-o.targetPorosity)/M_PI);
  std::vector<std::uint8_t> material(static_cast<std::size_t>(n)*n, 0);
  Counts count;
  auto index = [n](int x, int y) { return static_cast<std::size_t>(y)*n+x; };

  for (int y=0; y<n; ++y) for (int x=0; x<n; ++x) {
    const double px=x*o.dx, py=y*o.dx, r=std::hypot(px,py);
    std::uint8_t value=0;
    if (r <= electrolyteRadius) value=1;
    else if (r <= radius) {
      const double cx=(std::floor(px/o.particlePitch)+.5)*o.particlePitch;
      const double cy=(std::floor(py/o.particlePitch)+.5)*o.particlePitch;
      value = std::hypot(px-cx,py-cy) <= particleRadius ? 3 : 2;
    }
    material[index(x,y)]=value;
    if (value==0) ++count.outside;
    else if (value==1) ++count.electrolyte;
    else if (value==2) ++count.pore;
    else ++count.solid;
  }

  std::vector<std::uint8_t> visited(material.size(),0);
  std::queue<std::pair<int,int>> queue;
  for (int y=0; y<n; ++y) for (int x=0; x<n; ++x) {
    if (material[index(x,y)] != 2) continue;
    const double r=std::hypot(x*o.dx,y*o.dx);
    if (r <= electrolyteRadius+1.5*o.dx) {
      visited[index(x,y)]=1; queue.push({x,y});
    }
  }
  std::uint64_t connected=0;
  bool reachesOuter=false;
  constexpr int d[4][2]={{1,0},{-1,0},{0,1},{0,-1}};
  while (!queue.empty()) {
    auto [x,y]=queue.front(); queue.pop(); ++connected;
    if (std::hypot(x*o.dx,y*o.dx) >= radius-1.5*o.dx) reachesOuter=true;
    for (auto& step:d) {
      int nx=x+step[0], ny=y+step[1];
      if (nx<0||ny<0||nx>=n||ny>=n) continue;
      auto k=index(nx,ny);
      if (!visited[k] && material[k]==2) { visited[k]=1; queue.push({nx,ny}); }
    }
  }

  fs::create_directories(o.output);
  writeVtk(fs::path(o.output)/"quarter_sector_geometry.vtk",material,n,o.dx);
  const double cellArea=o.dx*o.dx;
  const double analyticOuterArea=M_PI*radius*radius/4.;
  const double analyticElectrolyteArea=M_PI*electrolyteRadius*electrolyteRadius/4.;
  const auto annulusCells=count.pore+count.solid;
  const double porosity=annulusCells ? static_cast<double>(count.pore)/annulusCells : 0.;
  const double rasterOuterArea=(count.electrolyte+annulusCells)*cellArea;
  const double rasterElectrolyteArea=count.electrolyte*cellArea;

  std::ofstream summary(fs::path(o.output)/"geometry_validation.csv");
  summary << std::setprecision(17)
    << "outer_radius_m,electrolyte_radius_m,dx_m,nx,ny,outer_cells,electrolyte_cells,pore_cells,solid_cells,"
       "outer_area_relative_error,electrolyte_area_relative_error,annulus_porosity,target_porosity,"
       "connected_pore_fraction,reaches_outer_arc\n"
    << radius << ',' << electrolyteRadius << ',' << o.dx << ',' << n << ',' << n << ','
    << count.electrolyte+annulusCells << ',' << count.electrolyte << ',' << count.pore << ',' << count.solid << ','
    << std::abs(rasterOuterArea-analyticOuterArea)/analyticOuterArea << ','
    << std::abs(rasterElectrolyteArea-analyticElectrolyteArea)/analyticElectrolyteArea << ','
    << porosity << ',' << o.targetPorosity << ','
    << (count.pore ? static_cast<double>(connected)/count.pore : 0.) << ',' << (reachesOuter?1:0) << '\n';

  std::ofstream cost(fs::path(o.output)/"time_cost.csv");
  cost << "dx_um,quarter_sector_cells,allocated_square_nodes,dt_s,steps_1s,steps_requested,estimated_memory_GiB,runtime_at_1_MLUPS_h,runtime_at_10_MLUPS_h,runtime_at_50_MLUPS_h\n";
  for (double dxUm : {10.,5.,2.,1.,.5,.2,.1}) {
    const double dx=dxUm*1.e-6;
    const double cells=analyticOuterArea/(dx*dx);
    const double sideNodes=std::ceil(radius/dx)+1.;
    const double allocatedNodes=sideNodes*sideNodes;
    const double dt=(o.tauLiquid-.5)*dx*dx/(3.*o.liquidNu);
    const double steps=o.simulatedTime/dt;
    const double updates=cells*steps;
    cost << std::setprecision(17) << dxUm << ',' << cells << ',' << allocatedNodes << ',' << dt << ',' << 1./dt << ',' << steps << ','
         << allocatedNodes*o.bytesPerAllocatedNode/std::pow(1024.,3) << ','
         << updates/1.e6/3600. << ',' << updates/1.e7/3600. << ',' << updates/5.e7/3600. << '\n';
  }

  const bool pass = std::abs(rasterOuterArea-analyticOuterArea)/analyticOuterArea < .005 &&
                    std::abs(rasterElectrolyteArea-analyticElectrolyteArea)/analyticElectrolyteArea < .01 &&
                    std::abs(porosity-o.targetPorosity) < .03 && reachesOuter;
  std::ofstream result(fs::path(o.output)/"validation_result.csv");
  result << "status,geometry_area_ok,electrolyte_area_ok,porosity_ok,pore_connected\n"
         << (pass?"PASS":"FAIL") << ','
         << (std::abs(rasterOuterArea-analyticOuterArea)/analyticOuterArea < .005) << ','
         << (std::abs(rasterElectrolyteArea-analyticElectrolyteArea)/analyticElectrolyteArea < .01) << ','
         << (std::abs(porosity-o.targetPorosity) < .03) << ',' << reachesOuter << '\n';
  std::cout << "status=" << (pass?"PASS":"FAIL") << " grid=" << n << 'x' << n
            << " annulus_porosity=" << porosity << " connected_fraction="
            << (count.pore ? static_cast<double>(connected)/count.pore : 0.) << '\n';
  return pass ? 0 : 2;
}
