#include "gravity_excitation.h"
#include "olb2D.h"
#include "olb2D.hh"
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include "geometry.h"
#include "physical_range_coupling.h"
#include "phase_projection.h"
#include "phase_outlet.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;
using T = double;
using NSD = D2Q9<RHO,NABLARHO,LAPLACERHO,FORCE,EXTERNAL_FORCE,TAU_EFF,STATISTIC,SCALAR>;
using CHD = D2Q9<CONV_POPS,FORCE,SOURCE,VELOCITY,OLD_PHIU,STATISTIC,THETA,BOUNDARY>;
using Op = PhysicalRangeLiangPostProcessor;
using Maps = meta::map<NavierStokes,descriptors::VALUED_DESCRIPTOR<T,NSD>,
                       Component1,descriptors::VALUED_DESCRIPTOR<T,CHD>>;
using Coupling = SuperLatticeCoupling<Op,Maps>;

struct Options {
  bool sideWalls=true;
  bool openOuter=true, boundMaterial=true, projectPhi=true;
  std::string collision="crlb",sourceMode="legacy";
  T phaseTau=.5075,magic=3./16;
  T radius=cell::width, dx=5e-6;
  T porosity=.40, tau=.5075, nu=1e-6, rhoL=1000, rhoG=1.2;
  T nuGas=1.5e-5;
  T sigma=1e-5, theta=30, width=4;
  T reservoirPressure=2500, pressureRamp=0;
  T gravityY=-9.81;
  T duration=.12;
  GravityExcitation excitation;
  bool explicitSteps=false, explicitDuration=false;
  std::string phaseOutlet="backflow";
  T outletSpeedMax=.1, outletSpeedMin=1e-6;
  T anomalySpeed=.30, anomalyPhiTolerance=.10;
  long steps=1600000, outputEvery=10000, diagnosticEvery=1;
  std::string output="results/default_run";
};

Options parse(int argc,char**argv) {
  Options o;
  for (int i=1;i<argc;++i) {
    if (i+1>=argc) throw std::runtime_error("missing option value");
    std::string key=argv[i],value=argv[++i];
    if(key=="--dx-um")o.dx=std::stod(value)*1e-6;
    else if(key=="--reservoir-pressure-pa" || key=="--electrolyte-pressure-pa")o.reservoirPressure=std::stod(value);
    else if(key=="--outer-boundary") {if(value!="open"&&value!="closed"&&value!="top-open")throw std::runtime_error("outer-boundary must be open, closed or top-open");o.openOuter=value!="closed";o.sideWalls=value=="top-open";}
    else if(key=="--anomaly-speed-lattice")o.anomalySpeed=std::stod(value);
    else if(key=="--anomaly-phi-tolerance")o.anomalyPhiTolerance=std::stod(value);
    else if(key=="--phase-outlet-mode")o.phaseOutlet=value;
    else if(key=="--outlet-speed-max")o.outletSpeedMax=std::stod(value);
    else if(key=="--outlet-speed-min")o.outletSpeedMin=std::stod(value);
    else if(key=="--duration-s"){o.duration=std::stod(value);o.explicitDuration=true;}
    else if(key=="--gravity-scale")o.excitation.scale=std::stod(value);
    else if(key=="--excitation-frequency-hz"||key=="--excitation-max-frequency-hz")o.excitation.frequency=std::stod(value);
    else if(key=="--excitation-min-frequency-hz")o.excitation.minFrequency=std::stod(value);
    else if(key=="--excitation-sweep-period-s")o.excitation.sweepPeriod=std::stod(value);
    else if(key=="--excitation-on-s")o.excitation.onTime=std::stod(value);
    else if(key=="--excitation-rest-s")o.excitation.restTime=std::stod(value);
    else if(key=="--gravity-y")o.gravityY=std::stod(value);
    else if(key=="--pressure-ramp-s")o.pressureRamp=std::stod(value);
    else if(key=="--ns-collision")o.collision=value;
    else if(key=="--source-mode")o.sourceMode=value;
    else if(key=="--phase-tau")o.phaseTau=std::stod(value);
    else if(key=="--nu-liquid")o.nu=std::stod(value);
    else if(key=="--phase-projection") {if(value!="off"&&value!="conservative")throw std::runtime_error("phase-projection must be off or conservative");o.projectPhi=value=="conservative";}
    else if(key=="--material-interpolation") {if(value!="linear"&&value!="bounded")throw std::runtime_error("material-interpolation must be linear or bounded");o.boundMaterial=value=="bounded";}
    else if(key=="--rho-gas")o.rhoG=std::stod(value);
    else if(key=="--nu-gas")o.nuGas=std::stod(value);
    else if(key=="--contact-angle-deg")o.theta=std::stod(value);
    else if(key=="--sigma")o.sigma=std::stod(value);
    else if(key=="--tau")o.tau=std::stod(value);
    else if(key=="--steps"){o.steps=std::stol(value);o.explicitSteps=true;}
    else if(key=="--output-every")o.outputEvery=std::stol(value);
    else if(key=="--diagnostic-every")o.diagnosticEvery=std::stol(value);
    else if(key=="--output")o.output=value;
    else throw std::runtime_error("unknown option: "+key);
  }
  o.excitation.validate();
  if(o.phaseOutlet!="backflow"&&o.phaseOutlet!="convective"&&o.phaseOutlet!="legacy")throw std::runtime_error("phase-outlet-mode: backflow, convective, legacy");
  if(!std::isfinite(o.outletSpeedMax)||!std::isfinite(o.outletSpeedMin)||o.outletSpeedMin<0||o.outletSpeedMax<=0||o.outletSpeedMin>o.outletSpeedMax)throw std::runtime_error("invalid outlet convection speed limits");
  if(o.explicitSteps&&o.explicitDuration)throw std::runtime_error("use either --steps or --duration-s");
  if(!std::isfinite(o.duration)||o.duration<=0)throw std::runtime_error("invalid duration or excitation parameters");
  if(!std::isfinite(o.anomalySpeed)||o.anomalySpeed<=0||!std::isfinite(o.anomalyPhiTolerance)||o.anomalyPhiTolerance<0)throw std::runtime_error("invalid anomaly thresholds");
  if(!std::isfinite(o.gravityY))throw std::runtime_error("gravity-y must be finite (m/s^2)");
  if(!std::isfinite(o.rhoG)||o.rhoG<=0||o.rhoG>=o.rhoL)throw std::runtime_error("rho-gas must be finite and between zero and liquid density");
  if(!std::isfinite(o.theta)||o.theta<=0||o.theta>=180)throw std::runtime_error("contact-angle-deg must be finite and in (0,180)");
  if(!std::isfinite(o.nu)||o.nu<=0||!std::isfinite(o.nuGas)||o.nuGas<=0)throw std::runtime_error("viscosities must be finite and positive");
  if(!std::isfinite(o.tau)||o.tau<=.5)throw std::runtime_error("tau must exceed 0.5");
  if (!std::isfinite(o.dx) || o.dx<=0 || o.dx>5e-6 || !std::isfinite(o.reservoirPressure) || o.steps<0 || o.outputEvery<=0 || o.diagnosticEvery<=0) throw std::runtime_error("invalid options: dx must be in (0,5] um; intervals positive");
  if(!std::isfinite(o.pressureRamp)||o.pressureRamp<0)throw std::runtime_error("pressure-ramp-s must be finite and nonnegative");
  if(o.collision!="bgk"&&o.collision!="trt"&&o.collision!="itrt"&&o.collision!="rlb"&&o.collision!="crlb")throw std::runtime_error("ns-collision: bgk, trt, itrt, rlb or crlb");
  if(o.sourceMode!="legacy"&&o.sourceMode!="equilibrium"&&o.sourceMode!="finite"&&o.sourceMode!="decoupled")throw std::runtime_error("source-mode: legacy, equilibrium, finite or decoupled");
  if(o.phaseTau==0)o.phaseTau=o.tau;
  if(!std::isfinite(o.phaseTau)||o.phaseTau<=.5||!std::isfinite(o.sigma)||o.sigma<0)throw std::runtime_error("invalid phase tau or sigma");
  if(o.sourceMode=="finite"&&o.reservoirPressure!=0)throw std::runtime_error("finite reservoir comparison requires pressure-pa 0");
  return o;
}

class SourcePhase : public AnalyticalF2D<T,T> {
  T _w;
public:
  SourcePhase(T w):AnalyticalF2D<T,T>(1),_w(w){}
  bool operator()(T out[],const T p[]) override {
    const T distance=std::max(std::abs(p[0]-cell::width/2)-cell::sourceWidth/2,p[1]-cell::sourceHeight);
    out[0]=.5*(1-std::tanh(2*distance/_w)); return true;
  }
};

class Phi : public SuperF2D<T,T> {
  SuperLattice<T,CHD>& _lattice;
  SuperGeometry<T,2>& _geometry;
  bool _raw;
public:
  Phi(SuperLattice<T,CHD>& lattice, SuperGeometry<T,2>& geometry, bool raw=false)
    :SuperF2D<T,T>(lattice,1),_lattice(lattice),_geometry(geometry),_raw(raw){this->getName()=raw?"phi_raw_populations":"phi_fluid_only";}
  bool operator()(T out[],const int p[]) override {
    auto& load=_lattice.getLoadBalancer();
    if(!load.isLocal(p[0]))return false;
    const int material=_geometry.getBlockGeometry(load.loc(p[0])).get(p[1],p[2]);
    if(material!=1 && material!=3 && !(material>=4&&material<=7)) {
      out[0]=std::numeric_limits<T>::quiet_NaN();
      return true;
    }
    auto cell=_lattice.getBlock(load.loc(p[0])).get(p[1],p[2]);
    out[0]=_raw?cell.computeRho():cell.template getFieldComponent<STATISTIC>(0);
    return true;
  }
};

struct Sample { T mean{},minimum{1e9},maximum{-1e9},cells{},rhoMin{1e9},tauMin{1e9},outside{},overshootL1{}; };
Sample sample(SuperGeometry<T,2>& geometry,SuperLattice<T,CHD>& ch,SuperLattice<T,NSD>& ns) {
  Sample s;
  for (int iC=0; iC<ch.getLoadBalancer().size(); ++iC) {
    auto& g=geometry.getBlockGeometry(iC); auto& lattice=ch.getBlock(iC);
    g.forSpatialLocations([&](LatticeR<2> r){
      int material=g.get(r);
      if(lattice.isPadding(r)||(material!=1&&material!=3&&!(material>=4&&material<=7)))return;
      T phi=lattice.get(r).template getFieldComponent<STATISTIC>(0);
      if(!std::isfinite(phi))s.mean=std::numeric_limits<T>::quiet_NaN();
      if(phi < -1e-12 || phi > 1.+1e-12)++s.outside;
      s.overshootL1+=std::abs(phi-std::clamp(phi,T(0),T(1)));
      T rho=ns.getBlock(iC).get(r).template getFieldComponent<RHO>(0);
      T tau=ns.getBlock(iC).get(r).template getFieldComponent<TAU_EFF>(0);
      if(!std::isfinite(rho)||!std::isfinite(tau))s.mean=std::numeric_limits<T>::quiet_NaN();
      s.rhoMin=std::min(s.rhoMin,rho);s.tauMin=std::min(s.tauMin,tau);
      s.minimum=std::min(s.minimum,phi);s.maximum=std::max(s.maximum,phi);
      if(material==1){s.mean+=phi;++s.cells;}
    });
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(s.outside,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.overshootL1,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.rhoMin,MPI_MIN);
  singleton::mpi().reduceAndBcast(s.tauMin,MPI_MIN);
  singleton::mpi().reduceAndBcast(s.mean,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.cells,MPI_SUM);
  singleton::mpi().reduceAndBcast(s.minimum,MPI_MIN);
  singleton::mpi().reduceAndBcast(s.maximum,MPI_MAX);
#endif
  if(s.cells)s.mean/=s.cells;
  return s;
}

PhaseProjectionStats projectPhaseLattice(SuperGeometry<T,2>& geometry,SuperLattice<T,CHD>& ch, T tolerance, const std::function<void(T,T,bool)>& onAnomaly, bool syncOutlets) {
  struct Location {int block;LatticeR<2> r;};
  std::vector<Location> locations;std::vector<double> values;
  locations.reserve(125000);values.reserve(125000);
  for(int b=0;b<ch.getLoadBalancer().size();++b){auto& g=geometry.getBlockGeometry(b);auto& lattice=ch.getBlock(b);
    g.forSpatialLocations([&](LatticeR<2> r){int m=g.get(r);if(lattice.isPadding(r)||(m!=1&&m!=4&&m!=5&&m!=7))return;
      locations.push_back({b,r});values.push_back(lattice.get(r).computeRho());});
  }
  T rawMin=std::numeric_limits<T>::infinity(),rawMax=-rawMin;
  int invalid=0;
  for(T v:values){if(!std::isfinite(v))invalid=1;else {rawMin=std::min(rawMin,v);rawMax=std::max(rawMax,v);}}
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(rawMin,MPI_MIN);
  singleton::mpi().reduceAndBcast(rawMax,MPI_MAX);
  singleton::mpi().reduceAndBcast(invalid,MPI_MAX);
#endif
  if(invalid||rawMin < -tolerance||rawMax > 1+tolerance)onAnomaly(rawMin,rawMax,invalid);
  auto stats=projectPhase(values,[](double v){
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(v,MPI_SUM);
#endif
    return v;
  });
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(stats.rawMin,MPI_MIN);singleton::mpi().reduceAndBcast(stats.rawMax,MPI_MAX);
#endif
  for(size_t i=0;i<values.size();++i){auto cell=ch.getBlock(locations[i].block).get(locations[i].r);T delta=values[i]-cell.computeRho();
    // Preserve the first population moment: sum_i w_i*c_i = 0.
    for(int q=0;q<CHD::q;++q)cell[q]+=descriptors::t<T,CHD>(q)*delta;
    if(syncOutlets){
      const int m=geometry.getBlockGeometry(locations[i].block).get(locations[i].r);
      if(m==4||m==5||m==7){auto saved=cell.template getField<CONV_POPS>();for(int q=0;q<CHD::q;++q)saved[q]=cell[q];cell.template setField<CONV_POPS>(saved);}
    }
  }
  return stats;
}

int main(int argc,char**argv) {
  initialize(&argc,&argv);
  if (argc==2 && (std::string(argv[1])=="--help" || std::string(argv[1])=="-h")) {
    std::cout << R"(Usage: ./run.sh [options] or ./main [options]
  --outer-boundary MODE         open, closed, top-open (left/right/bottom walls)
  --electrolyte-pressure-pa PA  Electrolyte source pressure (default 2500 Pa)
  --anomaly-speed-lattice VALUE  Immediate VTK speed threshold (default 0.30)
  --anomaly-phi-tolerance VALUE  Raw phi tolerance outside [0,1] (default 0.10)
  --gravity-y ACCEL             Signed gravity acceleration in m/s^2 (default -9.81)
  --phase-outlet-mode MODE      backflow (default), convective, legacy
  --outlet-speed-max VALUE      Convective update speed cap (default .1 lattice)
  --outlet-speed-min VALUE      Convective update speed floor (default 1e-6)
  --duration-s SEC              Requested physical duration (default 0.12 s)
  --steps N                     Update count instead of --duration-s
  --gravity-scale M            Peak excitation acceleration / 9.81 (default 1; 0 = no excitation)
  --excitation-min-frequency-hz HZ  Minimum instantaneous frequency (default 10 Hz)
  --excitation-max-frequency-hz HZ  Maximum instantaneous frequency (default 100 Hz)
  --excitation-sweep-period-s SEC   Sinusoidal frequency round trip (default 0.12 s)
  --excitation-frequency-hz HZ      Alias for maximum frequency
  --excitation-on-s SEC          Excitation window (default 0.005 s)
  --excitation-rest-s SEC        Rest window, gravity/pressure maintained (default 0.001 s)
  --pressure-ramp-s SEC         Linear pressure ramp (default 0 s = constant from start)
  --output-every N              VTK interval (default 10000)
  --diagnostic-every N          CSV interval (default 1)
  --output DIRECTORY            Output directory (run.sh creates a unique default)
  --dx-um UM                    Grid spacing (default 5)
  --tau VALUE                   Liquid relaxation time (default 0.5075)
  --phase-tau VALUE             Phase relaxation time (default 0.5075)
  --nu-liquid VALUE             Liquid kinematic viscosity (default 1e-6 m^2/s)
  --nu-gas VALUE                Gas kinematic viscosity (default 1.5e-5 m^2/s)
  --rho-gas VALUE               Gas density (default 1.2 kg/m^3)
  --contact-angle-deg VALUE     Contact angle (default 30)
  --sigma VALUE                 Surface tension (default 1e-5 N/m)
  --outer-boundary open|closed
  --ns-collision bgk|trt|itrt|rlb|crlb
  --source-mode legacy|equilibrium|finite|decoupled
  --phase-projection off|conservative
  --material-interpolation linear|bounded
Alias: --reservoir-pressure-pa is identical to --electrolyte-pressure-pa.
--steps controls update count, not physical dt. At defaults dt=6.25e-8 s.
Carrier phase continues through rest; only excitation acceleration is gated.
Physical dt=(tau-0.5)*dx^2/(3*nu-liquid); changing tau/dx changes stability.
)";
    return 0;
  }
  try { Options o=parse(argc,argv);
  std::filesystem::create_directories(o.output);
  singleton::directories().setOutputDir(o.output+"/");
  cell::Geometry layout(o.dx,o.porosity,50e-6,o.openOuter,o.sideWalls);
  if(singleton::mpi().isMainProcessor()) layout.write(o.output);
  IndicatorCuboid2D<T> domain({cell::width,layout.latticeHeight},{0,0});
  CuboidDecomposition2D<T> cuboids(domain, o.dx, singleton::mpi().getSize());
  HeuristicLoadBalancer<T> load(cuboids);
  SuperGeometry<T,2> geometry(cuboids,load,2); geometry.rename(0,1);
  for(int iC=0;iC<load.size();++iC) {
    auto& block=geometry.getBlockGeometry(iC);
    block.forSpatialLocations([&](LatticeR<2> r){
      const auto p=block.getPhysR(r);
      block.set(r,layout.material(p[0],p[1]));
    });
  }
  geometry.getStatistics().getStatisticsStatus()=true;geometry.communicate();geometry.checkForErrors();geometry.print();

  SuperLattice<T,NSD> ns(geometry); SuperLattice<T,CHD> ch(geometry);
  Coupling coupling(Op{},NavierStokes{},ns,Component1{},ch);
  const int resolution=static_cast<int>(std::lround(o.radius/o.dx));
  ns.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,NSD>>(
    resolution,o.tau,1.,o.radius,o.nu,o.rhoL);
  ch.setUnitConverter(ns.getUnitConverter()); const auto& converter=ns.getUnitConverter();
  if(!o.explicitSteps) {
    const T ratio=o.duration/converter.getPhysDeltaT();
    const T count=std::ceil(ratio-32*std::numeric_limits<T>::epsilon()*std::max(T(1),ratio));
    if(!std::isfinite(count)||count>T(std::numeric_limits<int>::max()))throw std::runtime_error("requested duration exceeds VTK step index limit");
    o.steps=static_cast<long>(count);
  }
  o.excitation.endTime=converter.getPhysTime(o.steps);
  const T sigmaLattice=o.sigma/converter.getConversionFactorSurfaceTension();
  const T rhoRatio=o.rhoG/o.rhoL;
  const T tauGasValue=.5+(o.nuGas/o.nu)*(o.tau-.5);
  auto pore=geometry.getMaterialIndicator(1),wall=geometry.getMaterialIndicator(2);
  auto reservoir=geometry.getMaterialIndicator(3),fluid=geometry.getMaterialIndicator({1,3,4,5,6,7});
  auto all=geometry.getMaterialIndicator({0,1,2,3,4,5,6,7});
  if(o.collision=="itrt")dynamics::set<MultiPhaseIncompressibleInterfaceTRTdynamics<T,NSD>>(ns,fluid);
  else if(o.collision=="crlb")dynamics::set<CorrectedMultiPhaseIncompressibleRLBdynamics<T,NSD>>(ns,fluid);
  else if(o.collision=="rlb")dynamics::set<MultiPhaseIncompressibleRLBdynamics<T,NSD>>(ns,fluid);
  else if(o.collision=="trt")dynamics::set<MultiPhaseIncompressibleTRTdynamics<T,NSD>>(ns,fluid);
  else dynamics::set<MultiPhaseIncompressibleBGKdynamics<T,NSD>>(ns,fluid);
  dynamics::set<ConservativePhaseFieldBGKdynamics<T,CHD>>(ch,fluid);
  boundary::set<boundary::BounceBackIncompressible>(ns,wall);
  boundary::set<boundary::PhaseFieldCurvedWall>(ch,wall);
  if(o.openOuter) {
    auto outlets=geometry.getMaterialIndicator({4,5,6,7});
    boundary::set<boundary::IncompressibleZouHePressure>(ns,outlets);
    if(o.phaseOutlet=="backflow") {
      dynamics::set<NoCollideDynamicsExternalVelocity>(ch,outlets);
      ch.addPostProcessor<stage::PostCollide>(geometry.getMaterialIndicator(4),meta::id<BackflowPhaseA<-1,0>>());
      ch.addPostProcessor<stage::PostStream>(geometry.getMaterialIndicator(4),meta::id<BackflowPhaseB<-1,0>>());
      ch.addPostProcessor<stage::PostCollide>(geometry.getMaterialIndicator(5),meta::id<BackflowPhaseA<1,0>>());
      ch.addPostProcessor<stage::PostStream>(geometry.getMaterialIndicator(5),meta::id<BackflowPhaseB<1,0>>());
      ch.addPostProcessor<stage::PostCollide>(geometry.getMaterialIndicator(7),meta::id<BackflowPhaseA<0,1>>());
      ch.addPostProcessor<stage::PostStream>(geometry.getMaterialIndicator(7),meta::id<BackflowPhaseB<0,1>>());
    } else setConvectivePhaseFieldBoundary<T,CHD>(ch,outlets);
    ch.setParameter<MAX_VELOCITY>(o.phaseOutlet=="legacy"?0.:o.outletSpeedMin);
    ch.addPostProcessor<stage::InitOutlet>(geometry.getMaterialIndicator(4),meta::id<SetOutletCells2D<-1,0>>());
    ch.addPostProcessor<stage::InitOutlet>(geometry.getMaterialIndicator(5),meta::id<SetOutletCells2D<1,0>>());
    ch.addPostProcessor<stage::InitOutlet>(geometry.getMaterialIndicator(7),meta::id<SetOutletCells2D<0,1>>());
  }
  ch.addPostProcessor<stage::PreCoupling>(fluid,meta::id<RhoStatistics>());
  if(o.sourceMode=="decoupled")coupling.restrictTo(geometry.getMaterialIndicator({1,4,5,6,7}));
  else coupling.restrictTo(fluid); coupling.setParameter<Op::SIGMA>(sigmaLattice);
  coupling.setParameter<Op::CORRECTED>(o.collision=="crlb"?1.:0.);
  coupling.setParameter<Op::BOUND_MATERIAL>(o.boundMaterial?1.:0.);
  coupling.setParameter<Op::W>(o.width); coupling.setParameter<Op::TAUS>({tauGasValue,o.tau});
  coupling.setParameter<Op::RHOS>({rhoRatio,T(1)}); coupling.setParameter<Op::SWITCH>(1.);
  ns.setParameter<OMEGA>(1/o.tau);ch.setParameter<OMEGA>(1/o.phaseTau);
  if(o.collision=="trt")ns.setParameter<collision::TRT::MAGIC>(o.magic);
  if(o.collision=="itrt") {
    ns.setParameter<collision::ITRT::TAU_MINUS>(1.5);
    ns.setParameter<collision::ITRT::MAXNABLARHO>((1-rhoRatio)/o.width);
  }
  ch.setParameter<INTERFACE_WIDTH>(o.width);ch.setParameter<THETA>(o.theta*M_PI/180);
  ch.setParameter<EPSILON>(1.5*o.width);
  auto& communicator=ch.getCommunicator(stage::PreCoupling());
  communicator.requestOverlap(1);communicator.requestField<STATISTIC>();communicator.exchangeRequests();

  struct OutletNode {int block;LatticeR<2> r;int nx,ny;};
  std::vector<OutletNode> outletNodes;
  if(o.openOuter)for(int b=0;b<load.size();++b){auto& g=geometry.getBlockGeometry(b);g.forSpatialLocations([&](LatticeR<2> r){if(ns.getBlock(b).isPadding(r))return;int m=g.get(r);if(m==4||m==5||m==7)outletNodes.push_back({b,r,m==4?-1:(m==5?1:0),m==7?1:0});});}
  struct OutletStats {T maxSpeed{},maxOutflow{},backflow{},outflow{},invalid{},minDenom{1e9};};
  auto scanOutlets=[&](){OutletStats out;
    for(const auto& n:outletNodes){
      auto cell=ns.getBlock(n.block).get(n.r);auto inner=cell.neighbor({-n.nx,-n.ny});
      T u[2]{},v[2]{};cell.computeU(u);inner.computeU(v);
      const T un=n.nx*v[0]+n.ny*v[1];
      if(!std::isfinite(u[0])||!std::isfinite(u[1])||!std::isfinite(v[0])||!std::isfinite(v[1]))++out.invalid;
      else {out.maxSpeed=std::max({out.maxSpeed,std::hypot(u[0],u[1]),std::hypot(v[0],v[1])});out.maxOutflow=std::max(out.maxOutflow,un);out.backflow+=un < -1e-12;out.outflow+=un>0;}
      const T denom=cell.template getFieldComponent<RHO>(0)+cell.template getFieldComponent<LAPLACERHO>(0)/4;
      if(!std::isfinite(denom))++out.invalid;else out.minDenom=std::min(out.minDenom,denom);
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(out.maxSpeed,MPI_MAX);singleton::mpi().reduceAndBcast(out.maxOutflow,MPI_MAX);
    singleton::mpi().reduceAndBcast(out.backflow,MPI_SUM);singleton::mpi().reduceAndBcast(out.outflow,MPI_SUM);
    singleton::mpi().reduceAndBcast(out.invalid,MPI_SUM);singleton::mpi().reduceAndBcast(out.minDenom,MPI_MIN);
#endif
    return out;
  };
  auto updateOutletSpeed=[&](const OutletStats& out){
    const T speed=o.phaseOutlet=="legacy"?T(0):std::clamp(out.maxOutflow,o.outletSpeedMin,o.outletSpeedMax);
    ch.setParameter<MAX_VELOCITY>(speed);return speed;
  };
  SourcePhase phase(o.width*o.dx);
  AnalyticalConst2D<T,T> zero(0),one(1),two(2),velocity(Vector<T,2>{0,0});
  // EXTERNAL_FORCE is acceleration: coupling adds it after dividing capillary force by density.
  const T gravityYLattice=o.gravityY*converter.getPhysDeltaT()*converter.getPhysDeltaT()/converter.getPhysDeltaX();
  AnalyticalConst2D<T,T> gravity(Vector<T,2>{0,gravityYLattice});
  fields::set<EXTERNAL_FORCE>(ns,fluid,gravity);
  constexpr T gReference=GravityExcitation::referenceG;
  const T accelerationAmplitude=o.excitation.peak();
  auto inertialAt=[&](T time){return o.excitation.at(time);};
  auto accelerationAt=[&](T time){return o.gravityY+inertialAt(time);};
  // Gate frequency-modulated excitation; retain integrated phase, gravity and pressure through rest.
  auto applyAcceleration=[&](T time){
    if(accelerationAmplitude==0)return;
    const T latticeAcceleration=accelerationAt(time)*converter.getPhysDeltaT()*converter.getPhysDeltaT()/converter.getPhysDeltaX();
    AnalyticalConst2D<T,T> acceleration(Vector<T,2>{0,latticeAcceleration});
    fields::set<EXTERNAL_FORCE>(ns,fluid,acceleration);
  };
  auto pressureAt=[&](T time) {return o.reservoirPressure*(o.pressureRamp>0?std::clamp(time/o.pressureRamp,T(0),T(1)):T(1));};
  auto applyPressure=[&](T time) {
    AnalyticalConst2D<T,T> p(pressureAt(time)/converter.getConversionFactorPressure());
    ns.iniEquilibrium(reservoir,p,velocity);
  };
  AnalyticalConst2D<T,T> rhoGas(rhoRatio),rhoDifference(1-rhoRatio),tauGas(tauGasValue),tauDifference(o.tau-tauGasValue);
  AnalyticalIdentity2D<T,T> rho(rhoGas+rhoDifference*phase),tau(tauGas+tauDifference*phase);
  fields::set<RHO>(ns,all,rho);fields::set<TAU_EFF>(ns,all,tau);
  fields::set<SCALAR>(ns,all,one);fields::set<SCALAR>(ns,wall,two);
  fields::set<BOUNDARY>(ch,fluid,one);fields::set<BOUNDARY>(ch,wall,two);
  AnalyticalConst2D<T,T> wallTheta(o.theta*M_PI/180);fields::set<THETA>(ch,wall,wallTheta);
  fields::set<OLD_PHIU>(ch,all,velocity);
  momenta::setOrderParameter(ch,all,phase);momenta::setIncompressiblePressure(ns,all,zero);
  momenta::setVelocity(ns,all,velocity);ns.iniEquilibrium(all,zero,velocity);ch.iniEquilibrium(all,phase,velocity);
  ns.initialize();ch.initialize();ch.iniEquilibrium(all,phase,velocity);ns.iniEquilibrium(all,zero,velocity);
  auto applySource=[&](T time) {
    if(o.sourceMode=="finite")return;
    applyPressure(time);
    if(o.sourceMode=="equilibrium"||o.sourceMode=="decoupled") {
      ch.iniEquilibrium(reservoir,one,velocity);
      fields::set<OLD_PHIU>(ch,reservoir,velocity);
      fields::set<VELOCITY>(ch,reservoir,velocity);
      fields::set<FORCE>(ch,reservoir,velocity);
      if(o.sourceMode=="decoupled") {
        AnalyticalConst2D<T,T> liquidTau(o.tau);
        fields::set<RHO>(ns,reservoir,one);
        fields::set<TAU_EFF>(ns,reservoir,liquidTau);
        fields::set<FORCE>(ns,reservoir,velocity);
        fields::set<NABLARHO>(ns,reservoir,velocity);
        fields::set<LAPLACERHO>(ns,reservoir,zero);
      }
    } else momenta::setOrderParameter(ch,reservoir,one);
  };
  momenta::setOrderParameter(ch,reservoir,one); // identical initial liquid bath in all ablations
  if(o.sourceMode!="finite")applySource(0);
  if(o.openOuter) {
    momenta::setIncompressiblePressure(ns,geometry.getMaterialIndicator({4,5,6,7}),zero);
    ch.executePostProcessors(stage::InitOutlet());
  }
  applyAcceleration(0); // initialize the prescribed acceleration before the first VTK output
  communicator.communicate();ns.setProcessingContext(ProcessingContext::Simulation);ch.setProcessingContext(ProcessingContext::Simulation);
  ch.executePostProcessors(stage::PreCoupling());
  communicator.communicate();

  if(singleton::mpi().isMainProcessor()) {
    std::ofstream cfg(o.output+"/run_config.csv");
    cfg<<std::setprecision(17)<<"key,value\nmodel_variant,regularized_phase_gravity_duty_5to1\nphase_outlet_mode,"<<o.phaseOutlet<<"\noutlet_speed_max,"<<o.outletSpeedMax<<"\noutlet_speed_min,"<<o.outletSpeedMin<<"\nphase_laplacian_normalization,abs_L_plus_epsilon\nrequested_duration_s,"<<(o.explicitSteps?converter.getPhysTime(o.steps):o.duration)<<"\npressure_duration_s,inf\npressure_schedule,constant\npressure_after_duration_pa,"<<o.reservoirPressure<<"\nexcitation_frequency_min_hz,"<<o.excitation.minFrequency<<"\nexcitation_frequency_max_hz,"<<o.excitation.frequency<<"\nexcitation_sweep_period_s,"<<o.excitation.sweepPeriod<<"\nexcitation_frequency_law,sinusoidal_min_max_min\nexcitation_mode,sinusoidal_frequency_sweep_on_rest\nexcitation_on_s,"<<o.excitation.onTime<<"\nexcitation_rest_s,"<<o.excitation.restTime<<"\nexcitation_set_period_s,"<<o.excitation.period()<<"\nexcitation_duty_fraction,"<<o.excitation.onTime/o.excitation.period()<<"\nexcitation_active_time_s,"<<o.excitation.activeTime(o.excitation.endTime)<<"\nexcitation_requested_sets,"<<o.excitation.endTime/o.excitation.period()<<"\nexcitation_control,acceleration_g\ngravity_scale,"<<o.excitation.scale<<"\nexcitation_g,"<<o.excitation.scale<<"\ng_reference_m_s2,"<<gReference<<"\nexcitation_acceleration_amplitude_m_s2,"<<accelerationAmplitude<<"\nexcitation_frame,rigid_porous_frame\nexcitation_phase_rad,0\nexcitation_start_s,0\nexcitation_duration_s,"<<(o.explicitSteps?converter.getPhysTime(o.steps):o.duration)<<"\nexcitation_end_s,"<<(o.explicitSteps?converter.getPhysTime(o.steps):o.duration)<<"\nexcitation_phase_policy,continuous_through_rest\nexcitation_acceleration_model,gated_sinusoidal_frequency_modulation\nanomaly_speed_lattice,"<<o.anomalySpeed<<"\nanomaly_phi_tolerance,"<<o.anomalyPhiTolerance<<"\ngravity_y_m_s2,"<<o.gravityY<<"\ngravity_y_lattice,"<<gravityYLattice<<"\nphase_projection,"<<(o.projectPhi?"conservative":"off")<<"\nmaterial_interpolation,"<<(o.boundMaterial?"bounded":"linear")<<"\nns_collision,"<<o.collision<<"\nsource_mode,"<<o.sourceMode<<"\nphase_tau,"<<o.phaseTau<<"\ntrt_magic,"<<o.magic<<"\nitrt_tau_minus,1.5\nitrt_max_rho_gradient,"<<(1-rhoRatio)/o.width<<"\nrequested_steps,"<<o.steps<<"\nouter_boundary,"<<(o.openOuter?"open":"closed")<<"\nside_boundary,"<<(o.sideWalls||!o.openOuter?"wall":"open")<<"\ntop_boundary,"<<(o.openOuter?"open":"wall")<<"\nbottom_boundary,wall\noutlet_pressure_pa,0\nmax_speed_lattice_limit,disabled\nwidth_m,"<<cell::width<<"\nheight_m,"<<layout.latticeHeight<<"\nrequested_height_m,"<<cell::height<<"\nsource_width_m,"<<cell::sourceWidth<<"\nsource_height_m,"<<cell::sourceHeight<<"\ngeometry_version,random_mixed_disks_wide_pore_v1\nparticle_min_diameter_m,"<<cell::Geometry::minDiameter<<"\nparticle_max_diameter_m,"<<cell::Geometry::maxDiameter<<"\nparticle_min_gap_m,"<<cell::Geometry::minGap<<"\ntarget_porosity,"<<o.porosity<<"\nseed,"<<cell::seed<<"\ndx_m,"<<o.dx<<"\ndt_s,"<<converter.getPhysDeltaT()<<"\nrho_liquid_kg_m3,"<<o.rhoL<<"\nrho_gas_kg_m3,"<<o.rhoG<<"\nnu_liquid_m2_s,"<<o.nu<<"\ntau_liquid,"<<o.tau<<"\ntau_gas,"<<tauGasValue<<"\nnu_gas_m2_s,"<<o.nuGas<<"\nsigma_N_m,"<<o.sigma<<"\ncontact_angle_deg,"<<o.theta<<"\npressure_ramp_s,"<<o.pressureRamp<<"\npressure_ramp_shape,linear\npressure_pa,"<<o.reservoirPressure<<"\npressure_conversion_Pa,"<<converter.getConversionFactorPressure()<<"\npressure_lattice,"<<o.reservoirPressure/converter.getConversionFactorPressure()<<"\n";
  }
  const bool isMain=singleton::mpi().isMainProcessor();
  std::ofstream csv,forcing,outletCsv;
  if(isMain) {
    outletCsv.open(o.output+"/outlet_diagnostics.csv");outletCsv<<std::setprecision(17)<<"step,time_s,convection_speed_lattice,max_outlet_and_inner_speed_lattice,max_outflow_normal_lattice,inflow_nodes,outflow_nodes,nonfinite_fields,min_corrected_density_lattice\n";
    forcing.open(o.output+"/forcing_history.csv");
    forcing<<std::setprecision(17)<<"step,time_s,source_pressure_pa,frame_acceleration_y_m_s2,inertial_acceleration_y_m_s2,effective_acceleration_y_m_s2,excitation_frequency_hz,phase_rad,gravity_scale,excitation_active\n";
    csv.open(o.output+"/multiphase_diagnostics.csv");
    csv<<std::setprecision(17)<<"step,time_s,reservoir_pressure_pa,pore_mean_phi,min_phi,max_phi,max_speed_lattice,finite,min_rho_lattice,min_tau_eff,out_of_range_phi_cells,phi_overshoot_l1,raw_min_phi,raw_max_phi,phase_projection_l1,phase_projection_mass_residual,phase_projection_cells\n";
  }
  SuperLatticePhysIncPressure2D<T,NSD> pressure(ns,converter);pressure.getName()="pressure_phys_Pa";
  Phi phi(ch,geometry);SuperLatticePhysVelocity2D<T,NSD> physVelocity(ns,converter);physVelocity.getName()="velocity_phys_m_s";
  SuperGeometryF<T,2> material(geometry);material.getName()="material";
  SuperLatticeField2D<T,NSD,EXTERNAL_FORCE> externalAcceleration(ns);externalAcceleration.getName()="external_acceleration_lattice";
  SuperLatticeField2D<T,NSD,RHO> rhoField(ns);rhoField.getName()="rho_lattice";
  SuperLatticeField2D<T,NSD,LAPLACERHO> laplaceRhoField(ns);laplaceRhoField.getName()="laplace_rho_lattice";
  SuperLatticeField2D<T,NSD,FORCE> forceField(ns);forceField.getName()="total_acceleration_lattice";
  SuperLatticeField2D<T,CHD,FORCE> phaseForceField(ch);phaseForceField.getName()="phase_force_lattice";
  SuperVTMwriter2D<T> writer("cell_horizontal");writer.addFunctor(pressure);writer.addFunctor(material);writer.addFunctor(phi);writer.addFunctor(physVelocity);writer.addFunctor(externalAcceleration);writer.addFunctor(rhoField);writer.addFunctor(laplaceRhoField);writer.addFunctor(forceField);writer.addFunctor(phaseForceField);writer.createMasterFile();
  Phi rawPhi(ch,geometry,true);
  SuperVTMwriter2D<T> rawWriter("anomaly_raw");
  rawWriter.addFunctor(rawPhi);rawWriter.addFunctor(material);rawWriter.createMasterFile();
  // Raw snapshots deliberately contain only current phase and geometry: NS fields
  // and STATISTIC have not yet been updated to this step at the projection stage.
  std::ofstream events;
  if(isMain){events.open(o.output+"/anomaly_events.csv");events<<"step,time_s,stage,reason,value_min,value_max\n";}
  auto event=[&](long step,const char* stage,const char* reason,T low,T high){
    if(isMain){events<<std::setprecision(17)<<step<<','<<converter.getPhysTime(step)<<','<<stage<<','<<reason<<','<<low<<','<<high<<'\n';events.flush();}
  };
  long lastRawStep=-1;
  auto writeRaw=[&](long step){if(lastRawStep!=step){rawWriter.write(static_cast<int>(step));lastRawStep=step;}};
  long completed=0;
  bool finite=true;
  PhaseProjectionStats projection;projection.rawMin=0;projection.rawMax=1;
  for(long step=0;step<=o.steps;++step){
    try {
    if(step){applySource(converter.getPhysTime(step-1));ns.collideAndStream();ch.collideAndStream();applySource(converter.getPhysTime(step));if(o.projectPhi)projection=projectPhaseLattice(geometry,ch,o.anomalyPhiTolerance,[&](T low,T high,bool invalid){
      event(step,"pre_projection",invalid?"nonfinite_phi":"raw_phi_range",low,high);
      writeRaw(step);
    },o.phaseOutlet!="legacy");communicator.communicate();ch.executePostProcessors(stage::PreCoupling());communicator.communicate();applyAcceleration(converter.getPhysTime(step));coupling.apply();}
    } catch(const std::exception& e) {
      event(step,"exception","simulation_exception",0,0);
      if(isMain)std::cerr<<"Exception at step "<<step<<": "<<e.what()<<std::endl;
      // Raw writer reads populations directly, including NaN/Inf, without needing
      // a successful projection or coupling update.
      writeRaw(step);
      event(step,"after_stream_before_projection","ns_state_on_exception",0,0);writer.write(static_cast<int>(step));
      throw;
    }
    completed=step;
    if(step%o.diagnosticEvery==0||step%o.outputEvery==0||step==o.steps){auto s=sample(geometry,ch,ns);T u=ns.getStatistics().getMaxU();bool ok=std::isfinite(s.mean)&&std::isfinite(s.minimum)&&std::isfinite(s.maximum)&&std::isfinite(u);ok=ok&&s.minimum>=-.1&&s.maximum<=1.1&&s.rhoMin>0&&s.tauMin>.5;finite&=ok;if(isMain)csv<<step<<','<<converter.getPhysTime(step)<<','<<pressureAt(converter.getPhysTime(step))<<','<<s.mean<<','<<s.minimum<<','<<s.maximum<<','<<u<<','<<ok<<','<<s.rhoMin<<','<<s.tauMin<<','<<s.outside<<','<<s.overshootL1<<','<<(o.projectPhi?projection.rawMin:s.minimum)<<','<<(o.projectPhi?projection.rawMax:s.maximum)<<','<<projection.l1<<','<<(projection.massAfter-projection.massBefore)<<','<<projection.changed<<'\n';csv.flush();}
    if(isMain&&(step%o.diagnosticEvery==0||step%o.outputEvery==0||step==o.steps)){
      const T time=converter.getPhysTime(step), inertial=inertialAt(time);
      forcing<<step<<','<<time<<','<<pressureAt(time)<<','<<-inertial<<','<<inertial<<','<<accelerationAt(time)<<','<<o.excitation.frequencyAt(time)<<','<<o.excitation.phase(time)<<','<<o.excitation.scale<<','<<o.excitation.active(time)<<'\n';forcing.flush();
    }
    const auto outletState=scanOutlets();const T convectiveSpeed=updateOutletSpeed(outletState);
    if(isMain&&(step%o.diagnosticEvery==0||step%o.outputEvery==0||step==o.steps)){
      outletCsv<<step<<','<<converter.getPhysTime(step)<<','<<convectiveSpeed<<','<<outletState.maxSpeed<<','<<outletState.maxOutflow<<','<<outletState.backflow<<','<<outletState.outflow<<','<<outletState.invalid<<','<<outletState.minDenom<<'\n';outletCsv.flush();
    }
    const T eventU=ns.getStatistics().getMaxU();
    const bool outletEvent=outletState.invalid>0||outletState.maxSpeed>.1;
    if(outletEvent)event(step,"post_coupling","outlet_velocity_or_nonfinite",outletState.minDenom,outletState.maxSpeed);
    if(outletState.invalid>0)finite=false;
    const bool speedEvent=!std::isfinite(eventU)||eventU>o.anomalySpeed;
    if(speedEvent)event(step,"post_coupling",std::isfinite(eventU)?"speed_threshold":"nonfinite_speed",eventU,eventU);
    if(step%o.outputEvery==0||step==o.steps||speedEvent||outletEvent)writer.write(static_cast<int>(step));
    if(!finite) {
      if(step%o.outputEvery!=0 && step!=o.steps && !speedEvent && !outletEvent)writer.write(static_cast<int>(step));
      break;
    }
  }
  auto final=sample(geometry,ch,ns);
  const bool pass=finite&&final.minimum>=-.1&&final.maximum<=1.1&&std::isfinite(ns.getStatistics().getMaxU());
  if(isMain) {
    std::ofstream result(o.output+"/multiphase_validation.csv");
    result<<"status,dx_m,dt_s,steps,physical_time_s,reservoir_pressure_pa,pore_mean_phi,min_phi,max_phi,max_speed_lattice\n"
          <<(pass?"PASS":"FAIL")<<','<<o.dx<<','<<converter.getPhysDeltaT()<<','<<completed<<','<<converter.getPhysTime(completed)<<','<<pressureAt(converter.getPhysTime(completed))<<','<<final.mean<<','<<final.minimum<<','<<final.maximum<<','<<ns.getStatistics().getMaxU()<<'\n';
    std::cout<<"multiphase_status="<<(pass?"PASS":"FAIL")<<" mean_phi="<<final.mean<<" phi_range=["<<final.minimum<<','<<final.maximum<<"] max_u="<<ns.getStatistics().getMaxU()<<'\n';
  }
  return pass?0:2;
  } catch(const std::exception& e) {std::cerr<<"NMC_longterm_case: "<<e.what()<<std::endl;return 1;}
}
